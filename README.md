# KsuccMoE datasets

Sequence datasets and saved data partitions used in **KsuccMoE: Biology-Aware Dynamic Routing of Multimodal Experts for Protein Lysine Succinylation Site Prediction**.

## Dataset overview

| Modification | Split | Positive sites | Negative sites | Total sites | Proteins |
|---|---|---:|---:|---:|---:|
| Lysine succinylation | Train | 4,743 | 4,526 | 9,269 | 2,191 |
| Lysine succinylation | Ind1 | 253 | 242 | 495 | 123 |
| Lysine succinylation | Ind2 | 253 | 2,971 | 3,224 | 123 |
| Mouse lysine lactylation | Train | 1,419 | 1,419 | 2,838 | 492 |
| Mouse lysine lactylation | Test | 355 | 355 | 710 | 331 |

These are the retained datasets used for the reported experiments. Labels and saved partitions have not been changed for this release.

## Files

```text
data/
  dataset_summary.csv
  provenance.json
  validation.json
  succinylation/
    train.csv
    ind1.csv
    ind2.csv
    full_sequences.fasta.gz
    sequence_sources.csv.gz
    splits/
      shared_folds.npz
      shared_folds.json
      fold_membership.csv.gz
  mouse_lactylation/
    train.csv
    test.csv
    full_sequences.fasta.gz
    source_windows/
      train_31.csv.gz
      test_31.csv.gz
    splits/
      shared_folds.npz
      shared_folds.json
      fold_membership.csv.gz
tools/
  verify_datasets.py
SHA256SUMS
```

### Site tables

Each row contains a protein `accession`, a **1-based** `position` in the archived full-length sequence, a binary `label` (1: annotated positive; 0: benchmark negative), and a 33-residue `window` centered on lysine. A hyphen (`-`) pads sequence termini. `sample_index` is **0-based** within the corresponding CSV. Negative labels follow the annotations and sampling used by the source benchmarks.

The succinylation training table also includes `outer_fold` (1–10) and `source_row_index`, the row index in the normalized, pre-deduplication input table. Mouse tables additionally retain the original 31-residue window, source, gene and protein name. Original mouse input tables are supplied in `source_windows/`.

### Full-length sequences

FASTA identifiers begin with the accession, followed by `|` metadata. The succinylation archive contains 2,314 recovered sequences from the original preparation pool, including proteins excluded during final filtering; the mouse archive contains 505 sequences. Every released site window has been checked against its archived sequence. The succinylation `sequence_sources.csv.gz` records sequence retrieval sources and hashes.

### Saved partitions and indexing

Training uses ten protein-grouped outer folds with grouped monitoring partitions inside each outer training fold. `fold_membership.csv.gz` exposes the exact saved partition membership using the released training CSV's **0-based `sample_index`**. Roles are `outer_train`, `outer_validation`, `selection_train` and `monitor`; selection and monitoring together partition the outer training subset.

The original `shared_folds.npz` files are also included without modification. Their `y`, `uniprot`, `position` and `outer_fold` arrays align with the released training table. **The NPZ partition-index arrays use the original input-row index space**, not necessarily the released CSV row index space. To convert an NPZ partition index into `sample_index`, look it up in the NPZ `original_index` array. For succinylation, this mapping is required because deduplication and protein filtering removed rows. For mouse lactylation, `original_index` is the identity mapping. Use the provided membership CSVs when working directly with the released tables.

### Relationships between splits

The succinylation training set shares no protein accessions with either test set. Ind1 is a subset of Ind2: they share the same 123 proteins and all 253 positive sites, with different numbers of negative sites.

The mouse experiment retains the supplied window-clustered split. Its training and test sets share 318 protein accessions, but no accession–position sites. This is a held-out-site evaluation. The original windows contain 31 residues with `X` padding; the model inputs were re-extracted from verified full-length sequences as 33-residue windows with `-` padding.

## Sources and citation

The succinylation training records and Ind1 originate from **CBiLSuccSite**; Ind2 originates from **LMSuccSite**. Protein identifiers, sequence versions and site coordinates were standardized, duplicate sites were removed, and training proteins present in the test sets were excluded.

- Tran et al. *Integrating CNN and Bi-LSTM for protein succinylation sites prediction based on Natural Language Processing technique.* Computers in Biology and Medicine (2025). [doi:10.1016/j.compbiomed.2025.109664](https://doi.org/10.1016/j.compbiomed.2025.109664).
- Pokharel et al. *Improving protein succinylation sites prediction using embeddings from protein language model.* Scientific Reports (2022). [doi:10.1038/s41598-022-21366-2](https://doi.org/10.1038/s41598-022-21366-2).

Mouse lactylation positives originate from UniProt/Swiss-Prot annotations and the published cardiac lactylome records below. The supplied dataset construction, balancing and split were retained.

- Wu et al. *Characterization of protein lactylation in healthy and ischemic mouse hearts.* Frontiers in Cardiovascular Medicine (2025). [doi:10.3389/fcvm.2025.1644886](https://doi.org/10.3389/fcvm.2025.1644886).
- [UniProt](https://www.uniprot.org/).

Please cite the relevant original sources when using these data. Original third-party data remain subject to their source terms; this release does not assign a new license to those data.

## Download and verify

Download the repository using **Code → Download ZIP**, or clone it:

```bash
git clone https://github.com/Cindy-ntu/KsuccMoe.git
cd KsuccMoe
python tools/verify_datasets.py
```

The verifier uses the Python standard library for checksums, counts, site identities and sequence windows. If NumPy is installed, it also checks the original saved fold files against the released membership tables. `SHA256SUMS` covers every release file except itself.

This release contains the sequence datasets and partition definitions. Pretrained model weights, derived embedding tensors, structure caches and training checkpoints are not included.

