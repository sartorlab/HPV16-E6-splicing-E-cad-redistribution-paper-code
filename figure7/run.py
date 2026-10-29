#!/usr/bin/env python3
"""Generate Figure 7B-F, or one selected panel, from patient-level source data.

Bailey F. Garb developed the statistical analyses. Shaomiao Xia adapted their
execution and rendering for this repository. The three bundled data tables
retain the endpoint and molecular-annotation snapshots used in those analyses.
An optional upstream score table replaces the numeric expression scores and
cell fractions by patient ID without changing clinical annotations.

Examples, run from the repository root:
    python3 figure7/run.py --output results/figure7
    python3 figure7/run.py --panel 7B --output results/figure7b
    python3 figure7/run.py --mode analysis --output results/statistics

Both ``all`` and ``figures`` recompute the statistics needed for their panels;
``analysis`` writes only numerical source results. Output directories must be
new or empty so a run cannot silently overwrite another set of figures.
"""

import argparse
import csv
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
PATIENT_FILES = (
    "survival_n230.csv",
    "subtype_integration_n230.csv",
    "network_n230.csv",
)


def read_table(path):
    """Read a CSV while preserving numeric strings and original column names."""
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames, list(reader)


def prepare_inputs(score_file, output):
    """Join regenerated scores to the frozen clinical tables by patient ID.

    All 230 IDs and all eleven expression scores plus seven cell fractions
    must be present. Only those numeric columns are replaced; event indicators,
    follow-up times, subgroup labels and adjustment variables remain unchanged.
    The resulting input files are saved alongside the derived statistics so
    readers can inspect the exact data used by their run.
    """
    if score_file is None:
        return ROOT / "data"

    fields, generated = read_table(score_file)
    if not fields or "samples" not in fields:
        raise ValueError("--patient-scores requires a samples column")
    scores = {row["samples"]: row for row in generated}
    _, network = read_table(ROOT / "data/network_n230.csv")
    required = set(network[0]) - {"samples", "batch", "HPV_integration"}
    if (
        len(generated) != 230
        or len(scores) != 230
        or set(scores) != {row["samples"] for row in network}
    ):
        raise ValueError("--patient-scores requires the same 230 unique patient IDs")
    missing = required - set(fields)
    if missing:
        raise ValueError("Missing patient-score columns: " + ", ".join(sorted(missing)))
    if any(
        not math.isfinite(float(row[column]))
        for row in generated
        for column in required
    ):
        raise ValueError("Patient-score columns must contain finite numeric values")

    inputs = output / "inputs_used"
    inputs.mkdir()
    for filename in PATIENT_FILES:
        columns, rows = read_table(ROOT / "data" / filename)
        for row in rows:
            for column in required & set(columns):
                row[column] = scores[row["samples"]][column]
        with (inputs / filename).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            writer.writerows(rows)
    return inputs


def main():
    """Resolve inputs, run the requested models, then write the selected PDFs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("figures", "analysis", "all"), default="all")
    parser.add_argument("--panel", choices=tuple("7" + letter for letter in "BCDEF"))
    parser.add_argument("--output", type=Path, required=True, help="New or empty directory")
    parser.add_argument("--rscript", default=os.environ.get("RSCRIPT", "Rscript"))
    parser.add_argument(
        "--patient-scores",
        type=Path,
        help="Regenerated eleven scores and seven cell fractions, keyed by samples",
    )
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        parser.error("--output must be a new or empty directory")
    executable = shutil.which(args.rscript)
    if executable is None:
        parser.error("Rscript executable not found: " + args.rscript)
    rscript = str(Path(executable).resolve())

    output.mkdir(parents=True, exist_ok=True)
    (output / "results").mkdir()
    if args.mode != "analysis":
        (output / "figures").mkdir()
    environment = os.environ.copy()
    # Direct module calls also use the restored library; explicit R_LIBS_USER wins.
    library_pointer = ROOT.parent / ".R_library_path"
    if "R_LIBS_USER" not in environment and library_pointer.is_file():
        library_text = library_pointer.read_text().strip()
        library = (ROOT.parent / Path(library_text).expanduser()).resolve()
        if not library_text or not library.is_dir():
            parser.error("The R library recorded in .R_library_path does not exist.")
        environment["R_LIBS_USER"] = str(library)
    environment.update(PYTHONDONTWRITEBYTECODE="1", TZ="UTC")
    environment.setdefault("MPLCONFIGDIR", str(output / ".matplotlib"))
    inputs = prepare_inputs(args.patient_scores, output)
    selection = args.panel or "all"

    analysis_scripts = []
    if selection != "7F":
        analysis_scripts.append("analyze_survival_and_groups.R")
    if selection in ("all", "7F"):
        analysis_scripts.append("analyze_network.R")
    for script in analysis_scripts:
        command = [
            rscript, "--vanilla", str(ROOT / "scripts" / script),
            str(ROOT), str(inputs), str(output / "results"), selection,
        ]
        subprocess.run(command, check=True, env=environment)

    if args.mode in ("figures", "all"):
        render_scripts = []
        if selection != "7F":
            render_scripts.append("render_survival_and_groups.py")
        if selection in ("all", "7F"):
            render_scripts.append("render_network.py")
        for script in render_scripts:
            subprocess.run(
                [sys.executable, "-B", str(ROOT / "scripts" / script), str(output), selection],
                check=True,
                env=environment,
            )
    print("Figure 7 output:", output)


if __name__ == "__main__":
    main()
