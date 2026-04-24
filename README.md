# PHACE for concatenated multi‑protein MSAs with artificial‑gap masking (PHACE-MultiProteins)

This repository contains a **modified PHACE** pipeline that can run on a **single concatenated MSA** composed of multiple protein MSAs (protein “blocks”), while **discarding false signal introduced by artificial gap blocks** (i.e., blocks filled with gaps because a species lacks the ortholog for that protein).

In the original PHACE workflow, gaps (`-`) in the alignment can be interpreted as real evolutionary events. In a concatenated multi‑protein MSA, however, entire protein blocks can be gaps purely due to *missing orthologs* (not biological indels). This version prevents those artificial gaps from contributing to tolerance scoring and coevolution scoring.

---

## What’s new compared to original PHACE

### 1) Metadata‑driven masking (two required files)
You must provide:

1. **Ortholog presence table** (`ortholog_selection_table.tsv`)
   - **Rows:** species names (must match tree tip labels)
   - **Columns:** protein identifiers (must match boundary file protein names)
   - **Cells:** non‑empty = ortholog present; empty/NA = ortholog missing (artificial gap block)

2. **Protein boundary map** (`proteinBoundriesForMSA.csv` / `.txt`)
   - Maps concatenated MSA coordinates to protein blocks
   - Must include (at minimum): `protein_name` (or `protein_id`), `start`, `end`
   - Coordinates are **1‑based, inclusive**, in the concatenated alignment.

### 2) Artificial gaps are treated as missing data (upstream)
When a species lacks the ortholog for the protein block containing a position, this version encodes that position as **`?`** (missing) in **MSA1** and **MSA2** encodings *for that species only* (dynamic, per position / per pair).

### 3) Missing‑data short‑circuit in change matrices
In `coev_diff_MSA1.R` and `coev_diff_MSA2.R`, if either parent or child ancestral state is `?`, the branch is forced to **“no change”** for that position (diff=0, change label `"--"`, etc.). This prevents spurious change labels like `C?`, `A?`, etc.

### 4) Pairwise scoring masks leaf **and internal** branches
For a position pair `(i1, i2)` mapping to proteins `(P1, P2)`:

- Species missing **P1 or P2** are masked (they contribute zero).
- Additionally, **any internal branch whose descendant subtree contains zero species that have both orthologs** is masked (branch weight set to 0 and per‑branch diffs zeroed).

### 5) Strict MSA ↔ tree consistency
Whenever a script reads both an MSA and a tree, the MSA is reordered to match `tree$tip.label` and the pipeline errors if the sets differ.

---

## Directory layout expected by the scripts

The scripts follow the same folder conventions as original PHACE. Typical layout:

```
PHACE_Codes/
  load_metadata.R
  ToleranceScore.R
  MSA1.R
  MSA2.R
  Part1_MSA1.R
  Part1_MSA2.R
  coev_diff_MSA1.R
  coev_diff_MSA2.R
  GetTotalChangeMatrix.R
  PHACE_parallel.R
  merge_PHACE_scores.py
  

Data/
  vals_MSA1.txt
  vals_MSA2.txt

ToleranceScores/
MSA1/
MSA2/
AncestralStates/
Part1_AC/
Part1_Gap/
totalChanges/
PHACE_scores/<id>/
```

(If a folder does not exist, create it before running the corresponding step.)

---

## Inputs you need

At minimum:

- **Concatenated AA alignment** (FASTA)
- **Tree** (Newick / IQ‑TREE `.treefile`) 
- **ASR output** (IQ‑TREE `.state`) 
- **Ortholog presence table** (`.tsv`)
- **Protein boundaries** (`.csv`/`.txt`)

> Partitioned tree search and ASR are optional, but recommended for concatenated multi-protein alignments. If partitioned ASR outputs are provided, the pipeline recognizes the partition structure and maps sites accordingly

For ASR of MSA1 and MSA2 you will also use:
- `Data/vals_MSA1.txt`
- `Data/vals_MSA2.txt`

---

## Step‑by‑step workflow

Below, `<id>` is your analysis identifier (used as filename prefix), and paths are examples.

### Step 0 — (Optional) Create folders
```
mkdir -p ToleranceScores MSA1 MSA2 AncestralStates Part1_AC Part1_Gap totalChanges PHACE_scores/<id>
```

### Step 1 — Tolerance scores (AA MSA + tree + ASR .state)

**Command**
```bash
Rscript PHACE_Codes/ToleranceScore.R \
  <id> \
  <concatenated_AA.fasta> \
  <tree.treefile> \
  <AA_ASR.state> \
  <boundaries.csv> \
  <ortholog_selection_table.tsv>
```

**Output**
- `ToleranceScores/<id>.csv`

> Note: This step is the only one that needs the AA ASR `.state` file.

---

### Step 2 — Build MSA1 encoding (C/A/-/?) and run ASR

**Build MSA1**
```bash
Rscript PHACE_Codes/MSA1.R \
  <id> \
  <concatenated_AA.fasta> \
  <boundaries.csv> \
  <ortholog_selection_table.tsv>
```

**Output**
- `MSA1/<id>_MSA1.fasta`

**Run IQ‑TREE2 ASR on MSA1**
```bash
iqtree2 -s MSA1/<id>_MSA1.fasta -te <tree.treefile> -blfix \
  -m Data/vals_MSA1.txt -asr --prefix AncestralStates/<id>_MSA1 --safe
```

This produces:
- `AncestralStates/<id>_MSA1.state`
- `AncestralStates/<id>_MSA1.treefile` (and other IQ‑TREE outputs)

---

### Step 3 — Build MSA2 encoding (C/G/?) and run ASR

**Build MSA2**
```bash
Rscript PHACE_Codes/MSA2.R \
  <id> \
  <concatenated_AA.fasta> \
  <boundaries.csv> \
  <ortholog_selection_table.tsv>
```

**Output**
- `MSA2/<id>_MSA2.fasta`

**Run IQ‑TREE2 ASR on MSA2**
```bash
iqtree2 -s MSA2/<id>_MSA2.fasta -te <tree.treefile> -blfix \
  -m Data/vals_MSA2.txt -asr --prefix AncestralStates/<id>_MSA2 --safe
```

---

### Step 4 — Build per‑position change matrices

**MSA1 (amino‑acid change mapping)**
```bash
Rscript PHACE_Codes/Part1_MSA1.R <id>
```

**MSA2 (gap change mapping)**
```bash
Rscript PHACE_Codes/Part1_MSA2.R <id>
```

Outputs go to:
- `Part1_AC/`
- `Part1_Gap/`

---

### Step 5 — Merge MSA1 + MSA2 into TotalChange matrix

```bash
Rscript PHACE_Codes/GetTotalChangeMatrix.R <id>
```

Output:
- `totalChanges/<id>_TotalChange.RData`

---

### Step 6 — Compute PHACE coevolution scores (parallel / SLURM array)

`PHACE_parallel.R` is designed to be run as a job array. Arguments:

1. `<id>`
2. `<array_task_id>` (1‑based)
3. `<num_jobs>` (total array size)
4. `<ortholog_selection_table.tsv>`
5. `<boundaries.csv>`
6. `<concatenated_AA.fasta>` (original AA fasta; used for position→protein mapping)

Example (single task):
```bash
Rscript PHACE_Codes/PHACE_parallel.R \
  <id> \
  1 \
  100 \
  <ortholog_selection_table.tsv> \
  <boundaries.csv> \
  <concatenated_AA.fasta>
```

Outputs:
- `PHACE_scores/<id>/<id>_PHACE_part<task>.RData`

---

### Step 7 — Merge array parts into final table (Python)

This repo includes `merge_PHACE_scores.py` which merges `PHACE_part*.RData` outputs into a single file and writes:
- `<id>_PHACE_internalBranchEffectRemoved_scores.feather`
- `<id>_PHACE_internalBranchEffectRemoved_scores.npz`

Usage:
```bash
python merge_PHACE_scores.py <id> <num_parts>
```

> Dependencies: `numpy`, `pandas`, `pyreadr`.  
> If you see a `NameError: sys is not defined`, add `import sys` at the top of `merge_PHACE_scores.py`.

---

## File format requirements (strict)

### Species naming
- MSA sequence names must match tree tip labels **exactly**.
- Species names in `ortholog_selection_table.tsv` must match those as well.
- This version will error out if the sets differ.

### Boundary ↔ ortholog table consistency
- Boundary `protein_name` values must match the ortholog table column names.
- If a protein in boundaries is not present as a column in the ortholog table, masking for that protein cannot be applied.

---

## Troubleshooting

### “Protein X not found in ortholog table columns”
Your boundary file contains a protein name not present in `ortholog_selection_table.tsv` header. Fix by making names identical in both files.

### “species labels do not match tree tip labels”
Your FASTA rownames and tree tip labels differ. Fix naming (or reorder/rename sequences) so they match exactly.

### Internal node labels missing
Internal‑branch masking requires that `totalChanges` branch labels can be mapped to tree nodes. In typical runs, `mat_info[,2]` contains tip labels and internal node labels. If your tree lacks internal node labels, you must either label internal nodes or change the mapping strategy.

## Citation
If you use PHACE in published research, cite the original PHACE paper / repository, and describe these modifications (multi‑protein concatenation and artificial‑gap masking) in Methods.
