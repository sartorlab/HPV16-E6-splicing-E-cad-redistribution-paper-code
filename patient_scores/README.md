This module supplies the eleven expression scores and seven cell fractions
used by Figure 7B–F. The analysis methods and gene lists are Bailey Garb's;
the portable implementation was organized by Shaomiao Xia.

Run with Python 3 and R containing `sva` and `edgeR`:

```bash
python3 patient_scores/run.py --mode all --output results/patient_scores
```

The output directory must be new or empty. Use `--rscript /path/to/Rscript`
when R is not on your PATH. The result `patient_scores_n230.csv` can be passed
to `figure7/run.py --patient-scores`. Its `samples` column uses the same pseudonymous
patient identifiers as the Figure 7 data. All input and output matrices are
uncompressed CSV files.

The four `data/counts_*.csv` files contain the final 230 patients' count-source
subsets. The analysis joins these sources by gene symbol, preserves repeated
symbols during the join, removes rows with at least 55 zero counts, and applies
ComBat_seq using eight sequencing batches without biological covariates. edgeR
then retains genes with CPM ≥1 in at least three patients, recomputes library
sizes, applies upperquartile normalization and calculates logCPM. The first
duplicate symbol is retained after normalization. The three intermediate count
and logCPM tables are written beside the final scores.

**Each gene is ranked across all 230 patients.** Rank-sum signatures combine
these ranks per patient, then standardize the sums across patients. The E6
percent-FL influence score instead uses `(sum of ranks − mean of ranks) / SD of
ranks` within each patient, with direction-specific minimum ranks. The scripts
explain tied ranks and missing-gene handling. Keep the full gene background and
original patient order throughout normalization and scoring.

`references/gene_lists/` contains the fixed E6 percent-FL, EMT, partial EMT,
T-cell, B-cell, keratinization, mesenchymal, cell-adhesion, respiratory-chain and
two KEGG pathway lists. The cell-adhesion membership is GO:0098609 (GOALL).
Its annotation and the Ensembl-to-symbol map use
[org.Hs.eg.db 3.17.0](https://bioconductor.org/packages/3.17/data/annotation/src/contrib/org.Hs.eg.db_3.17.0.tar.gz)
(Artistic-2.0). Applying the supplied E6 signature is separate from deriving it.
Third-party gene lists retain their original provenance and licensing.

Seven cell fractions are merged from saved results for the original 210
patients and **Job42_Adjusted** results for 20 added patients. Job42 records
describe a Puram HNSCC single-cell reference, 100 permutations, B-mode batch
correction and no quantile normalization. This package reproduces the saved
merge; it cannot rerun CIBERSORTx because the signature matrix, original210
settings and Job42 adjustment recipe are unavailable. Raw-read alignment and
count generation are outside this count-table analysis.
