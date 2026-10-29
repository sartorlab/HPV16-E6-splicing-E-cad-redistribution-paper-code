#!/usr/bin/env python3
"""Run the analyses for the October 27 manuscript's 16 included panels.

The ``figures`` route fits the Figure 7 models using the supplied patient scores
and plots S7/S8 using the supplied differential-expression/enrichment tables.
The ``all`` route additionally recalculates scores, DESeq2 and ORA from counts;
the exact external KEGG reference is required for that route.

Individual scripts in panels/ call this entry point with one panel selected.
Scientific calculations live in figure7/, patient_scores/ and cell_line/.
"""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
PANELS = (
    ["7" + letter for letter in "BCDEF"]
    + ["S7" + letter for letter in "ABCDEF"]
    + ["S8" + letter for letter in "ABCDE"]
)


def main(argv=None, panel=None):
    """Parse a workflow or single-panel request and run the needed modules.

    ``panel`` is set by the panel-specific entry scripts. ``argv`` permits calling
    the same command-line interface from Python. Generated files go to a new
    output folder, leaving source data intact.
    """
    parser = argparse.ArgumentParser(description=f"Generate Figure {panel}." if panel else __doc__)
    if panel is None:
        parser.add_argument(
            "mode", choices=["figures", "all", "patient-scores", "cell-line-de"],
            help="Plot supplied results, repeat all analyses, or run one upstream stage",
        )
        parser.add_argument("--panel", choices=PANELS, help="Generate one panel (figures mode)")
    else:
        parser.set_defaults(mode="figures", panel=panel)
    parser.add_argument("--output", type=Path, required=True, help="A new or empty output folder")
    parser.add_argument("--rscript", default=os.environ.get("RSCRIPT", "Rscript"))
    parser.add_argument("--r-library", type=Path, help="Optional R library restored from renv.lock")
    parser.add_argument("--gene-set-dir", type=Path, help="Exact external references for ORA recalculation")
    args = parser.parse_args(argv)

    # Fail before analysis when a requested input or interpreter is unavailable.
    if args.panel and args.mode != "figures":
        parser.error("--panel is used with figures mode; all recalculates the shared analyses.")
    executable = shutil.which(args.rscript)
    if executable is None:
        parser.error("Rscript was not found; install the documented environment or use --rscript.")
    rscript = str(Path(executable).resolve())
    if args.mode == "all" and args.gene_set_dir is None:
        parser.error("all requires --gene-set-dir; see cell_line/resources/gene_sets/README.md.")
    if args.gene_set_dir is not None and not args.gene_set_dir.is_dir():
        parser.error("--gene-set-dir must be an existing directory.")
    output = args.output.resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        parser.error("Output must be new or empty to preserve previous results.")

    # Apply the restored R environment to every child process. Absolute paths
    # allow the repository to be called from any working directory.
    env = os.environ.copy()
    library = args.r_library
    library_pointer = ROOT / ".R_library_path"
    if library is None and library_pointer.is_file():
        library = Path(library_pointer.read_text().strip())
    if library is not None:
        if not library.is_dir():
            parser.error("The selected R library does not exist.")
        env["R_LIBS_USER"] = str(library.resolve())
    env.update(PYTHON=sys.executable, PYTHONDONTWRITEBYTECODE="1", TZ="UTC")
    output.mkdir(parents=True, exist_ok=True)
    env.setdefault("MPLCONFIGDIR", str(output / ".matplotlib_cache"))

    def run(module, mode, extra=()):
        """Run a module in its own output subfolder and stop on a failed stage."""
        command = [
            sys.executable, str(ROOT / module / "run.py"),
            "--mode", mode, "--output", str(output / module),
            "--rscript", rscript, *map(str, extra),
        ]
        print(f"Running {module}: {mode}", flush=True)
        subprocess.run(command, cwd=output, env=env, check=True)

    # Figure 7 must consume newly calculated scores in a full run. Figure-only
    # runs instead use the supplied score tables.
    patient_scores = None
    if args.mode in ("all", "patient-scores"):
        run("patient_scores", "all")
        patient_scores = output / "patient_scores" / "patient_scores_n230.csv"
        if not patient_scores.is_file():
            raise RuntimeError("The patient-score stage did not produce its output table.")

    if args.mode in ("figures", "all"):
        if args.panel is None or args.panel.startswith("7"):
            extra = ["--patient-scores", patient_scores] if patient_scores else []
            if args.panel:
                extra += ["--panel", args.panel]
            run("figure7", "all", extra)
        if args.panel is None or args.panel.startswith("S"):
            extra = ["--gene-set-dir", args.gene_set_dir.resolve()] if args.mode == "all" else []
            if args.panel:
                extra += ["--panel", args.panel]
            run("cell_line", "all" if args.mode == "all" else "figures", extra)
    elif args.mode == "cell-line-de":
        run("cell_line", "de")

    # Require exactly one PDF per requested panel to catch incomplete runs.
    figures = sorted(output.rglob("*.pdf"))
    if args.mode in ("figures", "all"):
        selected = [args.panel] if args.panel else PANELS
        wanted = {"Figure" + name + ".pdf" for name in selected}
        if len(figures) != len(wanted) or {path.name for path in figures} != wanted:
            raise RuntimeError("The requested panel PDFs were not all generated; inspect the output.")
    print(f"Completed: {len(figures)} panel PDFs. Output: {output}")


if __name__ == "__main__":
    main()
