# UM-SCC-38 RNA-seq: Figures S7 and S8

Code and data for six volcano panels and five overrepresentation analysis (ORA)
panels. Differential expression analysis: Bailey F. Garb. ORA and figure workflow:
Shaomiao Xia.

| Comparison | Volcano | GO-BP / Hallmark / KEGG ORA |
|---|---|---|
| E6*I vs FL | S7A | S8A |
| E6*I vs FLm | S7B | S8B |
| E6*I vs PCD | S7C | S8C |
| FL vs PCD | S7D | S8D |
| FLm vs PCD | S7E | S8E |
| FL vs FLm | S7F | — |

Run from the repository root after setting up the documented environment:

```sh
# Draw one panel or all eleven from the supplied result tables.
python3 cell_line/run.py --mode figures --panel S7A --output /tmp/S7A
python3 cell_line/run.py --mode figures --panel S8A --output /tmp/S8A
python3 cell_line/run.py --mode figures --output /tmp/cell_figures

# Calculate all six DE comparisons from the count matrix.
python3 cell_line/run.py --mode de --output /tmp/cell_de

# Calculate DE and ORA, then draw all eleven panels.
python3 cell_line/run.py --mode all --output /tmp/cell_all \
  --gene-set-dir /path/to/frozen_gene_sets
```

`analysis` performs DE and ORA without drawing figures. `--panel` is available
only in `figures` mode. Both `analysis` and `all` require the complete reference
files described in the [gene-set instructions](resources/gene_sets/README.md).
Use `--rscript /path/to/Rscript` if needed. Each command requires a new or empty
output directory. PDFs appear in `figures/`, plotted values in
`figures/source_data/`, and calculated tables in `results/de/` and `results/ora/`.

All supplied data are readable TSV files. `inputs/gene_counts.tsv` contains
60,675 HTSeq-count genes × 12 samples; `config/sample_metadata.tsv` defines four
conditions with three independent biological replicates each. `data/de/` contains
six complete DE result tables and `data/ora/` contains all 40,180 tested ORA rows.
`resources/` supplies the annotation, Ensembl-to-Entrez mappings and GO-BP/Hallmark
definitions. Figure drawing uses the supplied DE/ORA tables; DE recalculation
starts from the count matrix. Full ORA recalculation additionally requires the
exact KEGG reference snapshot.

Plotting uses DejaVu Sans. Calculation and plotting steps are annotated in the
numbered R scripts.
