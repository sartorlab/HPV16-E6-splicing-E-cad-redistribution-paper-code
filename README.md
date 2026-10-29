# HPV16 E6 splicing: figure code and data

Code and data for **Figure 7B–F, Figure S7A–F and Figure S8A–E**, using the 230-patient dataset and UM-SCC-38 RNA-seq data. Original analyses: **Bailey F. Garb**. Revision analyses and code preparation: **Shaomiao Xia**. Responsible group author: **Maureen A. Sartor**.

## Find a panel

Each entry script names its input tables, analysis scripts and method. Related panels share the same scientific functions. Count matrices and analysis result tables are uncompressed TSV or CSV files.

| Panel | Analysis | Entry script |
|---|---|---|
| 7B | All-HPV overall and recurrence-free survival | [Figure7B.py](panels/Figure7B.py) |
| 7C | IMU/KRT score comparison | [Figure7C.py](panels/Figure7C.py) |
| 7D | HPV integration score comparison | [Figure7D.py](panels/Figure7D.py) |
| 7E | HPV16 overall and recurrence-free survival | [Figure7E.py](panels/Figure7E.py) |
| 7F | Selected 16-association network | [Figure7F.py](panels/Figure7F.py) |
| S7A | E6*I versus FL volcano | [FigureS7A.py](panels/FigureS7A.py) |
| S7B | E6*I versus FLm volcano | [FigureS7B.py](panels/FigureS7B.py) |
| S7C | E6*I versus PCD volcano | [FigureS7C.py](panels/FigureS7C.py) |
| S7D | FL versus PCD volcano | [FigureS7D.py](panels/FigureS7D.py) |
| S7E | FLm versus PCD volcano | [FigureS7E.py](panels/FigureS7E.py) |
| S7F | FL versus FLm volcano | [FigureS7F.py](panels/FigureS7F.py) |
| S8A | E6*I versus FL pathway enrichment | [FigureS8A.py](panels/FigureS8A.py) |
| S8B | E6*I versus FLm pathway enrichment | [FigureS8B.py](panels/FigureS8B.py) |
| S8C | E6*I versus PCD pathway enrichment | [FigureS8C.py](panels/FigureS8C.py) |
| S8D | FL versus PCD pathway enrichment | [FigureS8D.py](panels/FigureS8D.py) |
| S8E | FLm versus PCD pathway enrichment | [FigureS8E.py](panels/FigureS8E.py) |

[panels.tsv](panels.tsv) gives the input/output paths for every panel. Module instructions are in [figure7](figure7/README.md), [cell_line](cell_line/README.md) and [patient_scores](patient_scores/README.md).

## Install

Use Python 3.11.14 and R 4.4.3 on Linux. Python package versions are in `requirements.txt`; R package versions are in `renv.lock`. Run the following commands from the repository root:

```bash
conda env create --file environment.yml
conda activate e6splice-panels
python -m pip install --requirement requirements.txt
Rscript --vanilla restore_environment.R .
```

The environment file supplies the system tools for PDF rendering. Plotting uses DejaVu Sans through R's `systemfonts` package.

## Generate figures

Generate a single panel, or all 16:

```bash
python panels/Figure7B.py --output build/Figure7B
python panels/FigureS8A.py --output build/FigureS8A
python reproduce.py figures --output build/all_figures
```

Each output directory must be new or empty. PDFs appear under `figure7/figures/` or `cell_line/figures/` within that directory, alongside the numerical tables used to draw them. The Figure 7 models are fitted from the supplied patient inputs; S7/S8 use the supplied complete DESeq2/ORA result tables.

To recalculate patient scores and cell-line differential expression and enrichment from count tables, supply the complete gene-set files described in the [reference instructions](cell_line/resources/gene_sets/README.md):

```bash
python reproduce.py all --output build/from_counts \
  --gene-set-dir /path/to/exact_gene_sets
```

This route passes newly calculated patient scores into Figure 7 and uses the supplied cell-fraction tables. Separate `patient-scores` and `cell-line-de` modes are also available. Use `--rscript` or `--r-library` to select an existing installation.

## Citation and license

See [CITATION.cff](CITATION.cff) for authorship and citation details. The [MIT license](LICENSE) covers the code. Patient data and third-party references retain their source terms.
