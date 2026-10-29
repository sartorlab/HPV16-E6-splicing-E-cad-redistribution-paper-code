"""Draw Figure 7B-E from the R-generated numerical results.

The survival panels show Kaplan-Meier steps, confidence bands, censor marks,
log-rank tests, continuous-score adjusted Cox effects, and numbers at risk.
C/D show the observed patient scores with a violin density, deterministic point
jitter, median, and explicitly labelled Student t test. Bailey's statistical
choices and survival colors are retained; Shaomiao Xia adapted the rendering.

Usage: python render_survival_and_groups.py OUTPUT [all|7B|7C|7D|7E]
Each selected panel is written as one vector PDF in OUTPUT/figures.
"""

from pathlib import Path
import csv
import hashlib
import os
import sys

os.environ.setdefault("MPLCONFIGDIR", "/tmp/e6_publication_matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(sys.argv[1]).resolve()
OUT = ROOT / "figures"
COLORS = {"E6FL Low": "#185FA5", "E6FL High": "#993C1D"}
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "pdf.fonttype": 42,
    "svg.fonttype": "none", "axes.spines.top": False, "axes.spines.right": False,
})


def read(name):
    """Load numerical source results without rounding their stored values."""
    with (ROOT / "results" / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def save(fig, name):
    """Write a vector PDF with author credit in its document metadata."""
    fig.savefig(OUT / (name + ".pdf"), metadata={
        "Title": name,
        "Author": "Bailey (statistical analysis); Shaomiao (publication adaptation)",
    })
    plt.close(fig)


def pvalue(value):
    """Format p values for legible labels without changing analysis precision."""
    value = float(value)
    return f"{value:.3g}" if value < .001 else f"{value:.4f}"


def survival_panel(panel):
    """Draw OS and RFS curves and risk tables for one HPV-defined population."""
    stats = [row for row in read("Figure7BE_survival_statistics.csv") if row["panel"] == panel]
    steps = [row for row in read("Figure7BE_KM_steps.csv") if row["panel"] == panel]
    risks = [row for row in read("Figure7BE_risk_tables.csv") if row["panel"] == panel]
    fig = plt.figure(figsize=(180 / 25.4, 99 / 25.4))
    grid = fig.add_gridspec(2, 2, height_ratios=[3.55, 1], hspace=.48, wspace=.35)
    fig.text(.015, .98, panel[-1], fontweight="bold", fontsize=12, va="top")
    for index, endpoint in enumerate(("OS", "RFS")):
        statistic = next(row for row in stats if row["endpoint"] == endpoint)
        ax = fig.add_subplot(grid[0, index])
        risk_ax = fig.add_subplot(grid[1, index])
        for group in ("E6FL Low", "E6FL High"):
            rows = [row for row in steps if row["endpoint"] == endpoint and row["group"] == group]
            x = np.r_[0, [float(row["month"]) for row in rows]]
            y = np.r_[1, [float(row["survival"]) for row in rows]]
            lower = np.r_[1, [float(row["lower"]) if row["lower"] else np.nan for row in rows]]
            upper = np.r_[1, [float(row["upper"]) if row["upper"] else np.nan for row in rows]]

            # Keep the final step constant through the administrative horizon;
            # the curve itself is determined entirely by the R survival fit.
            if x[-1] < 60:
                x = np.r_[x, 60]
                y = np.r_[y, y[-1]]
                lower = np.r_[lower, lower[-1]]
                upper = np.r_[upper, upper[-1]]
            ax.step(x, y, where="post", color=COLORS[group], lw=1.25, label=group)
            ax.fill_between(x, lower, upper, step="post", color=COLORS[group], alpha=.12, lw=0)
            censored = [row for row in rows if int(row["n_censor"]) > 0]
            ax.plot(
                [float(row["month"]) for row in censored],
                [float(row["survival"]) for row in censored],
                ls="none", marker="+", markersize=4, markeredgewidth=.7, color=COLORS[group],
            )
        title = "Overall survival" if endpoint == "OS" else "Recurrence-free survival"
        if panel == "Figure7E":
            title += " (HPV16+)"
        ax.set(
            title=title, xlim=(0, 60), ylim=(0, 1.04), xticks=range(0, 61, 12),
            xlabel="Time (months)", ylabel="Survival probability",
        )
        ax.tick_params(labelsize=7)
        text = (
            f"Log-rank p = {pvalue(statistic['logrank_p'])}\n"
            f"Adjusted Cox p = {pvalue(statistic['cox_p'])}\n"
            f"HR per 50 points = {float(statistic['HR_per_50_points']):.2f} "
            f"({float(statistic['HR50_ci_low']):.2f}–{float(statistic['HR50_ci_high']):.2f})\n"
            f"KM n = {statistic['km_n']}; Cox n = {statistic['cox_n']}"
        )
        ax.text(.03, .04, text, transform=ax.transAxes, va="bottom", fontsize=6.7, linespacing=1.35)
        ax.legend(frameon=False, loc="lower right", fontsize=6.4, handlelength=1.2)
        risk_ax.set(xlim=(-17, 62), ylim=(-.55, 2.05))
        risk_ax.axis("off")
        risk_ax.text(-17, 1.85, "Number at risk", fontsize=7)
        for index, group in enumerate(("E6FL Low", "E6FL High")):
            y = 1 - index
            risk_ax.text(-17, y, group.replace("E6FL ", "").capitalize(), color=COLORS[group], va="center", fontsize=7)
            for row in risks:
                if row["endpoint"] == endpoint and row["group"] == group:
                    risk_ax.text(float(row["month"]), y, row["n_risk"], ha="center", va="center", fontsize=7)
        for month in range(0, 61, 12):
            risk_ax.text(month, -.55, str(month), ha="center", va="center", fontsize=6.5, color="#4b5965")
    fig.subplots_adjust(left=.08, right=.985, top=.91, bottom=.08)
    save(fig, panel)


def comparison_panel(panel):
    """Draw a molecular-group comparison using the explicitly chosen Student test."""
    data = read("Figure7CD_plot_data.csv")
    tests = read("Figure7CD_ttests.csv")
    if panel == "Figure7C":
        groups, labels = ("IMU", "KRT"), ("IMU", "KRT")
    else:
        groups = ("Negative", "Positive")
        labels = ("HPV integration\nnegative", "HPV integration\npositive")
    fig, ax = plt.subplots(figsize=(83 / 25.4, 84 / 25.4))
    test = next(row for row in tests if row["panel"] == panel)
    maximum = 0
    for index, group in enumerate(groups):
        rows = [row for row in data if row["panel"] == panel and row["group"] == group]
        values = np.array([float(row["E6FL_E6ALL_score"]) for row in rows])
        violin = ax.violinplot([values], positions=[index], widths=.65, showextrema=False)
        for body in violin["bodies"]:
            body.set_facecolor(("#43a047", "#ef8a29")[index])
            body.set_edgecolor("black")
            body.set_linewidth(.8)
            body.set_alpha(.8)

        # A stable ID hash supplies horizontal jitter only; score values and
        # statistical calculations are never randomized or altered.
        jitter = np.array([
            (int(hashlib.sha256(row["sample_id"].encode()).hexdigest()[:8], 16) / 0xffffffff - .5) * .35
            for row in rows
        ])
        ax.scatter(index + jitter, values, s=10, alpha=.65, color="#525b66", linewidths=0)
        ax.plot([index - .23, index + .23], [np.median(values)] * 2, color="#111111", lw=1.5)
        maximum = max(maximum, max(values))
    ax.plot([0, 0, 1, 1], [maximum + 5, maximum + 12, maximum + 12, maximum + 5], lw=.8, color="black")
    ax.text(.5, maximum + 16, "p = " + pvalue(test["p_value"]), ha="center", va="bottom", fontsize=8)
    ax.set(
        xticks=[0, 1],
        xticklabels=[labels[0] + f"\n(n = {test['n1']})", labels[1] + f"\n(n = {test['n2']})"],
        xlim=(-.55, 1.55), ylim=(0, max(800, maximum + 75)), ylabel="E6FL:E6ALL influence score",
    )
    ax.text(-.18, 1.06, panel[-1], transform=ax.transAxes, fontweight="bold", fontsize=12)
    fig.text(.54, .025, "Two-sided Student unpaired t-test", ha="center", fontsize=6.5)
    fig.subplots_adjust(left=.24, right=.97, top=.89, bottom=.24)
    save(fig, panel)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    selection = sys.argv[2] if len(sys.argv) > 2 else "all"
    for letter in "BCDE":
        if selection not in ("all", "7" + letter):
            continue
        panel = "Figure7" + letter
        if letter in "BE":
            survival_panel(panel)
        else:
            comparison_panel(panel)
    print("Rendered the selected Figure 7B-E PDF panels.")
