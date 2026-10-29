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
default `--mode all` calculate statistics and draw panels. `--patient-scores PATH`
uses the upstream regenerated eleven scores and seven cell fractions, matched
by patient ID, while preserving the endpoint and annotation snapshots.

The three input tables share 230 pseudonymous `samples` IDs. Survival times are
in months, age in years, and events are 1 (event) or 0 (censored). Subtype is
TRUE (IMU) or FALSE (KRT); integration is `neg` or `pos`. Missing values are empty
CSV fields. Cell fractions use a 0-1 scale; expression scores retain their
original units. `data/network_candidates.csv` names the 18 specified network
candidates and their definitions. Detailed calculation steps are commented in
the R analysis scripts and Python plotting functions.

For survival, the high/low split is the population mean **before** endpoint
exclusions. Follow-up is administratively censored at 60 months. The adjusted
Cox model uses continuous E6 score, age, pooled AJCC stage, and five clinical
cohorts. For C/D, the 204 patients with observed annotations are compared using
a two-sided Student test. The original variance setting was not recovered, so
this equal-variance choice is stated explicitly in the panels. For F, the 18 specified candidates are fitted
with **eight sequencing batches** as adjustment variables. Integration uses a
logistic model; other candidates use linear models. Nominal p < 0.05 selects
16 edges (10 positive, 6 negative). Color shows coefficient sign; line width
increases with -log10(p).

The regenerated figures use model-derived labels. In the October manuscript,
HR50 labels for B-OS/B-RFS/E-OS/E-RFS are 0.78/0.74/0.78/0.69; the model-derived
values round to 0.76/0.73/0.79/0.65. At-risk values differ for B-RFS high at
48 months (8 versus calculated 6), and E-RFS low at 48/60 months (27/13 versus
19/5). Cox p values agree. The F methods text says cohort adjustment, whereas
the supplied analysis uses eight sequencing batches. These discrepancies should
be reconciled in the manuscript; the code does not reproduce manually entered
labels that conflict with its calculations.
