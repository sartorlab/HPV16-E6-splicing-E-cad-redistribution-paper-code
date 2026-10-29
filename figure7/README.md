# Figure 7B-F

Bailey F. Garb developed the patient scoring, regression, and survival analyses.
Shaomiao Xia adapted their execution and figure rendering for this repository.

| Panel | Analysis | Input table |
|---|---|---|
| 7B | Overall and recurrence-free survival, all HPV-positive patients | `data/survival_n230.csv` |
| 7C | IMU versus KRT E6 score | `data/subtype_integration_n230.csv` |
| 7D | HPV integration-negative versus positive E6 score | `data/subtype_integration_n230.csv` |
| 7E | Overall and recurrence-free survival, HPV16 patients | `data/survival_n230.csv` |
| 7F | Sixteen batch-adjusted associations | `data/network_n230.csv` |

From the repository root, generate all five panels or one panel:

```bash
python3 figure7/run.py --output results/figure7
python3 figure7/run.py --panel 7B --output results/figure7b
```

Use `--rscript /path/to/Rscript` if needed. Output must be new or empty. Each
selected PDF is written to `figures/`, with its numerical source results in
`results/`. `--mode analysis` writes only statistics; `--mode figures` and the
default `--mode all` calculate statistics and draw panels.

To use regenerated patient scores, run the following from the repository root:

```bash
python3 patient_scores/run.py --output results/patient_scores
python3 figure7/run.py --output results/figure7_from_counts \
  --patient-scores results/patient_scores/patient_scores_n230.csv
```

The score table supplies eleven expression scores and seven cell fractions,
matched to the clinical tables by patient ID.

The three input tables share 230 pseudonymous `samples` IDs. Survival times are
in months, age in years, and events are 1 (event) or 0 (censored). Subtype is
TRUE (IMU) or FALSE (KRT); integration is `neg` or `pos`. Missing values are empty
CSV fields. Cell fractions use a 0-1 scale; expression scores retain their
original units. `data/network_candidates.csv` lists the network candidates.
Calculation steps are annotated in the R analysis scripts and Python plotting
functions.
