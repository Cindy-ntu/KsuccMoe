"""Verify released KsuccMoE datasets without downloading or modifying anything."""
from pathlib import Path
import csv
import gzip
import hashlib
import io

ROOT = Path(__file__).resolve().parents[1]


def check(condition, message):
    if not condition:
        raise ValueError(message)


def rows(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle))


def sequences(path):
    result = {}
    accession = None
    with gzip.open(path, 'rt', encoding='utf-8-sig') as handle:
        for line in handle:
            if line.startswith('>'):
                accession = line[1:].split('|')[0].split()[0]
                check(accession not in result, f'Duplicate FASTA identifier: {accession}')
                result[accession] = ''
            elif line.strip():
                result[accession] += line.strip()
    return result


def identity(row):
    return row['accession'], int(row['position']), int(row['label'])


def main():
    manifest = (ROOT / 'SHA256SUMS').read_text().splitlines()
    for line in manifest:
        expected, rel = line.split('  ', 1)
        check(hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() == expected,
              f'Checksum mismatch: {rel}')
    expected_paths = {line.split('  ', 1)[1] for line in manifest}
    check(expected_paths == {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*')
                              if p.is_file() and p.name != 'SHA256SUMS' and '.git' not in p.relative_to(ROOT).parts
                              and '__pycache__' not in p.parts}, 'Release inventory mismatch')
    tables = {}
    for item in rows(ROOT / 'data/dataset_summary.csv'):
        name, split = item['dataset'], item['split']
        table = rows(ROOT / f'data/{name}/{split}.csv')
        seqs = sequences(ROOT / f'data/{name}/full_sequences.fasta.gz')
        check(len(table) == int(item['sites']), f'{name}/{split}: site count')
        check(sum(int(r['label']) for r in table) == int(item['positive']), f'{name}/{split}: positives')
        check(sum(int(r['label']) == 0 for r in table) == int(item['negative']), f'{name}/{split}: negatives')
        check(len({r['accession'] for r in table}) == int(item['proteins']), f'{name}/{split}: protein count')
        check(len({identity(r) for r in table}) == len(table), f'{name}/{split}: duplicate sites')
        check(len({(r['accession'], r['position']) for r in table}) == len(table), f'{name}/{split}: conflicting labels')
        for i, row in enumerate(table):
            check(int(row['sample_index']) == i, 'Sample order mismatch')
            seq, pos = seqs[row['accession']], int(row['position'])
            check(int(row['label']) in (0, 1), 'Invalid label')
            check(1 <= pos <= len(seq) and seq[pos-1] == 'K', 'Site is not lysine')
            window = ''.join(seq[j] if 0 <= j < len(seq) else '-' for j in range(pos-17, pos+16))
            check(row['window'] == window, f'Window mismatch at {row["accession"]}:{pos}')
            if 'original_window31' in row:
                original = ''.join(seq[j] if 0 <= j < len(seq) else 'X' for j in range(pos-16, pos+15))
                check(row['original_window31'] == original, 'Original mouse window mismatch')
        tables[name, split] = table
    ktrain, ind1, ind2 = (tables['succinylation', split] for split in ('train', 'ind1', 'ind2'))
    check(not {r['accession'] for r in ktrain} & {r['accession'] for r in ind2}, 'Succinylation protein overlap')
    check({identity(r) for r in ind1} <= {identity(r) for r in ind2}, 'Ind1 is not a subset of Ind2')
    mtrain, mtest = (tables['mouse_lactylation', split] for split in ('train', 'test'))
    check(len({r['accession'] for r in mtrain} & {r['accession'] for r in mtest}) == 318, 'Mouse protein overlap count')
    check(not {identity(r) for r in mtrain} & {identity(r) for r in mtest}, 'Mouse site overlap')
    try:
        import numpy as np
    except ImportError:
        np = None
    for name in ('succinylation', 'mouse_lactylation'):
        table = tables[name, 'train']
        membership = rows(ROOT / f'data/{name}/splits/fold_membership.csv.gz')
        groups = {}
        for row in membership:
            key = int(row['fold']), row['role']
            groups.setdefault(key, []).append(int(row['sample_index']))
        if np is not None:
            archive = np.load(ROOT / f'data/{name}/splits/shared_folds.npz', allow_pickle=False)
            convert = {int(v): i for i, v in enumerate(archive['original_index'])}
            check([identity(r) for r in table] == [(str(a), int(p), int(y)) for a, p, y in
                  zip(archive['uniprot'], archive['position'], archive['y'])], 'Protocol identity mismatch')
        for fold in range(1, 11):
            role_sets = {role: set(groups[fold, role]) for role in
                         ('outer_train', 'outer_validation', 'selection_train', 'monitor')}
            for role, indices in role_sets.items():
                check(len(groups[fold, role]) == len(indices), 'Duplicate partition indices')
                if np is not None:
                    check(groups[fold, role] == [convert[int(v)] for v in archive[f'{role}_{fold:02d}']], 'Original partition mismatch')
            for a, b in (('outer_train', 'outer_validation'), ('selection_train', 'monitor')):
                check(role_sets[a].isdisjoint(role_sets[b]), 'Partition overlap')
                check({table[i]['accession'] for i in role_sets[a]}.isdisjoint(
                      {table[i]['accession'] for i in role_sets[b]}), 'Within-fold protein overlap')
            check(role_sets['outer_train'] | role_sets['outer_validation'] == set(range(len(table))), 'Outer coverage mismatch')
            check(role_sets['selection_train'] | role_sets['monitor'] == role_sets['outer_train'], 'Monitor coverage mismatch')
        if np is not None:
            archive.close()
    print(f'PASS: {len(manifest)} file checksums; five dataset splits; '
          f'{sum(map(len, tables.values()))} sites; sequence windows; saved grouped partitions.'
          + (' Original NPZ files verified.' if np is not None else ' Install NumPy to also verify original NPZ arrays.'))


if __name__ == '__main__':
    main()
