# UM-SCC-38 RNA-seq: Figures S7 and S8

Code and data for six volcano panels and five overrepresentation analysis (ORA)
panels. The DESeq2 analysis retains Bailey F. Garb's model and contrasts; Shaomiao
Xia prepared the manuscript ORA and figure workflow.

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
  --gene-set-dir /path/to/authorized/frozen_gene_sets
```

`analysis` performs DE and ORA without drawing figures. `--panel` is available
only in `figures` mode. Use `--rscript /path/to/Rscript` if needed. Each command
requires a new output directory; PDFs appear in `figures/`, and the exact plotted
values appear in `figures/source_data/`.

All supplied data are readable TSV files. `inputs/gene_counts.tsv` contains
60,675 HTSeq-count genes × 12 samples; `config/sample_metadata.tsv` defines four
conditions with three independent biological replicates each. `data/de/` contains
six complete DE result tables and `data/ora/` contains all 40,180 tested ORA rows.
`resources/` supplies the annotation, Ensembl-to-Entrez mappings and available
pathway definitions. The numbered R scripts document their inputs, calculations,
assumptions and outputs.

The `~cell_line` DESeq2 model retains genes with ≥10 counts in ≥3 samples
(17,725 genes). `results()` uses `alpha=0.1` for independent filtering; reported
DEGs require `padj < 0.05` and `abs(log2FoldChange) > 1`. Comparisons to PCD use
the PCD-first fitted contrast followed by LFC/Wald sign reversal; displayed
positive values mean greater expression in the stated numerator.

ORA selects one unambiguous Entrez representative per gene by greatest
all-sample baseMean, with Ensembl identifier breaking ties. For each collection,
the background contains represented genes with finite DESeq2 adjusted P values
and an annotation in that collection. Pathways with 15–500 background genes enter
separate up/down hypergeometric tests and BH correction (30 families). Zero-hit
pathways remain in the correction. S8 shows up to ten significant pathways per
direction/database, ranked by adjusted P, raw P and pathway identifier.

Figure drawing and DE recalculation use supplied data. Full ORA recalculation
requires the exact authorized KEGG snapshot described in
[reference instructions](resources/gene_sets/README.md). Raw sequencing reads
are not included; their accession has not yet been assigned. Plotting uses
installed DejaVu Sans; font substitutions can change label positions.
