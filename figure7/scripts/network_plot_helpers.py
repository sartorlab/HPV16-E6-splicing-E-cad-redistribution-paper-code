"""Drawing helpers for the Bailey-style Figure 7F association network.

The companion render_network.py supplies labels and fixed node coordinates.
This module reads the fitted candidate models from the selected output directory,
encodes coefficient sign as color and -log10(p) as line width, and writes the
network PDF plus its plotted edge values. It does not fit statistical models.
"""

from pathlib import Path
import csv
import math
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MODULE = Path(__file__).resolve().parents[1]
ROOT = Path(sys.argv[1]).resolve()
SOURCE = ROOT / "results/Figure7F_candidate_models.csv"
OUT = ROOT / "figures"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "pdf.fonttype": 42,
    "mathtext.fontset": "custom", "mathtext.rm": "DejaVu Sans",
    "mathtext.it": "DejaVu Sans:italic", "mathtext.bf": "DejaVu Sans:bold",
    "mathtext.cal": "DejaVu Sans",
})
with SOURCE.open(newline="") as handle:
    ALL = list(csv.DictReader(handle))
DATA = {row["Variable"]: row for row in ALL if float(row["P_value"]) < .05}
RED, BLUE, INK = "#ff0000", "#0000ff", "#242B32"
HUB = r"$\mathrm{E6FL}\!:\!\mathrm{E6}_{\mathrm{ALL}}$" + "\nscore"
LABELS = {}


def color(variable):
    """Encode a positive model coefficient as red and a negative one as blue."""
    return RED if float(DATA[variable]["Beta"]) > 0 else BLUE


def width(variable):
    """Map the selected edge's -log10(p) linearly onto 0.4-2.3 points."""
    return .4 + 1.9 * (EVIDENCE[variable] - LO) / (HI - LO)


def pwidth(p):
    """Apply the same width scale to the legend's minimum and maximum p values."""
    return .4 + 1.9 * (-math.log10(p) - LO) / (HI - LO)


def canvas(height):
    """Create a 180 mm wide panel using millimetres as drawing coordinates."""
    fig = plt.figure(figsize=(180 / 25.4, height / 25.4), dpi=180)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, 180), ylim=(0, height))
    ax.axis("off")
    return fig, ax


def text(ax, x, y, label, size=9.8, **kwargs):
    """Place a node or legend label with the common font and spacing settings."""
    defaults = dict(ha="left", va="center", color=INK, fontsize=size, linespacing=1.04, zorder=5)
    defaults.update(kwargs)
    return ax.text(x, y, label, **defaults)


def edge_row(variable, x, y):
    """Export the exact fitted values and visual encodings used for one edge."""
    row = DATA[variable]
    return {
        "Variable": variable, "Label": LABELS[variable].replace("\n", " "),
        "N": row["N"], "Beta": row["Beta"], "P_value": row["P_value"],
        "Direction": "Positive" if float(row["Beta"]) > 0 else "Negative",
        "Color": color(variable), "Line_width_pt": width(variable), "x_mm": x, "y_mm": y,
    }


def save(fig, edges):
    """Write the publication panel and its compact scientific source-data table."""
    fig.savefig(OUT / "Figure7F.pdf", facecolor="white", metadata={
        "Title": "Figure7F", "Author": "Bailey F. Garb (analysis); Shaomiao Xia (adaptation)",
        "Subject": "16 nominal p<0.05 E6FL:E6ALL associations from 18 fixed display candidates",
    })
    plt.close(fig)
    with (ROOT / "results/Figure7F_edges.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(edges[0]))
        writer.writeheader()
        writer.writerows(edges)
