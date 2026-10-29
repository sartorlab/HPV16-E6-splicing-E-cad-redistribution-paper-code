# Patient scores for Figure 7B–F

Calculate eleven expression scores for the 230 patients and combine them with
seven supplied cell fractions. Analysis methods and gene lists: Bailey F. Garb.
Code preparation: Shaomiao Xia.

Run from the repository root after installing the [documented environment](../README.md#install):

```bash
python3 patient_scores/run.py --mode all --output results/patient_scores
```

The output directory must be new or empty. Use `--rscript /path/to/Rscript`
to select an R installation. Pass the resulting scores to Figure 7 with:

```bash
python3 figure7/run.py --output results/figure7_from_counts \
  --patient-scores results/patient_scores/patient_scores_n230.csv
```

Inputs are the four `data/counts_*.csv` matrices, `data/cohort_n230.csv`,
`data/cell_fractions_original_n210.csv`, `data/cell_fractions_added_n20.csv`,
and the annotation and gene lists in `references/`.

The workflow joins counts, applies ComBat-seq batch correction, calculates
normalized logCPM and expression scores, and merges the supplied cell fractions.
Output CSV files are:

- `uncorrected_joined_counts_n230.csv`
- `batch_corrected_counts_n230.csv`
- `log_cpm_n230.csv`
- `patient_scores_n230.csv`

The final table contains one row per patient, with the same `samples` identifiers
used by the Figure 7 inputs. Calculation steps are annotated in `scripts/`.

The annotation and Ensembl-to-symbol mapping use
[org.Hs.eg.db 3.17.0](https://bioconductor.org/packages/3.17/data/annotation/src/contrib/org.Hs.eg.db_3.17.0.tar.gz)
(Artistic-2.0). Third-party gene lists retain their source attribution and licenses.
