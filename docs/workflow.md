## 1. Adding new MOTU table(s)

1.  Create a folder under `data/raw/` named exactly like your experiment ID\
    (e.g. `data/raw/EXP003/`).

2.  Drop the **four CSV files** and name them\
    `EXP003_reads.csv`, `EXP003_motu.csv`, `EXP003_samples.csv`, `EXP003_stations.csv`.

3.  In R (can read multiple tables at a time):

\`\`\`r source("src/data_preprocessing/read_meta.R") meta \<- read_meta("EXP003")

## 2. Pooled and rarefied MOTU tables

PCR/extraction replicates are pooled so there is one column per sample. The tables are rarefied to 10,000 reads (default) and only MOTUs remaining in rarefied data retained for further processing. Samples with \<30,000 reads are discarded.

Each run of `rarefy_motus()` writes one CSV per experiment:`data/processed/<experiment>_rarefied_motu.csv`

### File structure

| column | description | example |
|-------------------|------------------------------------|-------------------|
| **id** | MOTU identifier (`reads$id` / `motu$id`) | `MOTU_00051` |
| *others* | One column per **sample_name** (after pooling PCR replicates) | `fjord_A`, `fjord_B` |

-   **id** is always the **first** column\
-   All remaining columns are named exactly after the `sample_name` values found in\
    the corresponding `samples.csv` file for that experiment.\
-   Cell values are read counts after rarefaction to the chosen depth (default = 30 000).

### Example

| Id         | fjord_A | fjord_B | fjord_C |
|------------|--------:|--------:|--------:|
| MOTU_00001 |   1 234 |     567 |      89 |
| MOTU_00002 |     456 |     321 |      78 |
| …          |       … |       … |       … |

## 4. Assigning MINs to new MOTU tables

After rarefaction each experiment has a file\
`data/processed/<experiment>_motu_list.csv` that holds the sequences present in that run.\
`assign_min.py` compares every sequence against the master reference, assigns an existing **MIN** (97 % identity, `vsearch`), or creates a new MIN when no match is found.

### Prerequisites

| item | notes |
|---------------------------------|---------------------------------------|
| **vsearch** | Installed in `/usr/local/bin` (or adjust the script). |
| `pandas` | Listed in `requirements.txt`; `pip install -r requirements.txt`. |
| Folder layout | `data/processed/` for inputs/outputs, `data/reference/` for the growing DB. |

### One-line usage

\`\`\`bash \# from the repo root python src/MIN_clustering.py

### Workflow

| step | action |
|-----------|-------------------------------------------------------------|
| 1 | Scan `data/processed/` for `*_motu_list.csv`. |
| 2 | For every sequence:<br> • write a temp FASTA<br> • run `vsearch --usearch_global` vs `data/reference/OTU_reference.fasta` at 0.97 ID<br> • parse best hit (if any). |
| 3 | **No hit?**<br> • Create a new MIN `MZP_#########` (9-digit counter).<br> • Append the sequence to `OTU_reference.fasta`.<br> • Append a row to `OTU_database_updated.csv` with timestamp and source experiment. |
| 4 | Write `<experiment>_assigned.csv` (two columns: `id`, `MIN`) to `data/processed/`. |
| 5 | Append the source filename to `assigned_motu_files.txt` so it is not re-processed. |
| 6 | Save the updated reference CSV (`OTU_database_updated.csv`). |

### Input and output file

| path | description |
|---------------------------------------|----------------------------------------------|
| `data/processed/<exp>_motu_list.csv` | **Input** table from rarefaction step (`id`, `sequence`). |
| `data/reference/OTU_reference.fasta` | FASTA of all sequences that already have a MIN. Auto-created on first run. |
| `data/reference/OTU_database_updated.csv` | Tabular version of the reference (MIN, sequence, timestamp, source). |
| `data/processed/<exp>_assigned.csv` | **Output** – MOTU ids mapped to MINs (`id`, `MIN`). |
| `data/processed/assigned_motu_files.txt` | Log of files already processed. |

### Example

| id         | MIN           |
|------------|---------------|
| MOTU_00001 | MZP_000000123 |
| MOTU_00002 | MZP_000000045 |
| …          | …             |

### Re-running

You can safely rerun the command at any time; previously processed MOTU tables will be skipped. If you delete a line from assigned_motu_files.txt, that table will be re-processed and the reference DB will only grow (never overwritten).
