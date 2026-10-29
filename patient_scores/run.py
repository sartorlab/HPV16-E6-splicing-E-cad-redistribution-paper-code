#!/usr/bin/env python3
"""Generate the 230-patient score and cell-fraction inputs for Figure 7B-F.

The analysis follows Bailey Garb's count processing and scoring methods.
Python runs the count join and R analysis, then combines the eleven calculated
scores with seven saved CIBERSORTx cell fractions using pseudonymous patient IDs.
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
SCORE_COLUMNS = [
    'e6%FL_influence_scores', 'mes_score', 'tcell_score', 'cell_adhesion',
    'go_resp_electron_transport', 'p_emt_zt', 'kegg_ox_phos_score',
    'kerat_score', 'emt_score', 'bcell_score', 'kegg_focal_adhesion',
]
FRACTION_COLUMNS = [
    'Endothelial', 'Malignant', 'T cells CD4', 'Macrophage',
    'Fibroblast', 'T cells CD8', 'Dendritic',
]


def read_rows(path):
    """Read CSV rows as strings, preserving cell-fraction decimal values."""
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle))


def index_patients(rows, required_columns):
    """Require unique patient IDs and finite values before joining tables."""
    indexed = {}
    for row in rows:
        sample = row['samples']
        if not sample or sample in indexed:
            raise ValueError('Patient identifiers must be present and unique.')
        for column in required_columns:
            if not math.isfinite(float(row[column])):
                raise ValueError(f'Non-finite patient value in {column}.')
        indexed[sample] = row
    return indexed


def assemble_fractions(output):
    """Join saved original210 + adjusted20 fractions to the calculated scores.

    The added patients use Job42_Adjusted. This operation preserves the saved
    deconvolution results; it does not rerun the CIBERSORTx service or recreate
    the adjustment. Values remain strings so CSV writing does not round them.
    """
    old = read_rows(ROOT/'data/cell_fractions_original_n210.csv')
    added = read_rows(ROOT/'data/cell_fractions_added_n20.csv')
    cohort = read_rows(ROOT/'data/cohort_n230.csv')
    if len(old) != 210 or len(added) != 20:
        raise ValueError('Expected 210 original and 20 added patient fractions.')

    added_ids = {row['samples'] for row in added}
    cohort_added = {
        row['samples'] for row in cohort if row['cohort_version'] == 'added20'
    }
    if added_ids != cohort_added:
        raise ValueError('Added patient fractions differ from cohort metadata.')

    fractions = index_patients(old+added, FRACTION_COLUMNS)
    score_file = output/'patient_scores_n230.csv'
    scores = index_patients(read_rows(score_file), SCORE_COLUMNS)
    cohort_ids = {row['samples'] for row in cohort}
    if len(scores) != 230 or set(scores) != set(fractions) or set(scores) != cohort_ids:
        raise ValueError('Scores, cell fractions and cohort must contain the same 230 patients.')

    # Sorting happens only after all across-patient ranks have been calculated.
    # The expression matrix's original patient order is preserved during scoring.
    with score_file.open('w', newline='') as handle:
        writer = csv.DictWriter(
            handle, fieldnames=['samples']+SCORE_COLUMNS+FRACTION_COLUMNS
        )
        writer.writeheader()
        for sample in sorted(scores):
            writer.writerow(dict(scores[sample], **fractions[sample]))


def main():
    """Run each dependency in order and write all results to a new directory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['all'], default='all')
    parser.add_argument('--output', type=Path, required=True,
                        help='New or empty directory for generated CSV tables.')
    parser.add_argument('--rscript', default=os.environ.get('RSCRIPT', 'Rscript'))
    args = parser.parse_args()

    output = args.output.resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        parser.error('--output must be a new or empty directory.')
    rscript = shutil.which(args.rscript)
    if rscript is None:
        parser.error('Rscript not found; supply --rscript /path/to/Rscript.')
    rscript = str(Path(rscript).resolve())
    output.mkdir(parents=True, exist_ok=True)

    # Pass the module location explicitly: R's --file argument can encode spaces.
    env = os.environ.copy()
    # Direct module calls also use the restored library; explicit R_LIBS_USER wins.
    library_pointer = ROOT.parent / '.R_library_path'
    if 'R_LIBS_USER' not in env and library_pointer.is_file():
        library_text = library_pointer.read_text().strip()
        library = (ROOT.parent / Path(library_text).expanduser()).resolve()
        if not library_text or not library.is_dir():
            parser.error('The R library recorded in .R_library_path does not exist.')
        env['R_LIBS_USER'] = str(library)
    env['E6_PATIENT_MODULE_ROOT'] = str(ROOT)
    stages = [
        '00_prepare_counts.py',
        '00_batch_correct.R',
        '01_normalize_logcpm.R',
        '02_calculate_scores.R',
    ]
    for name in stages:
        command = [sys.executable] if name.endswith('.py') else [rscript, '--vanilla']
        subprocess.run(
            command+[str(ROOT/'scripts'/name), str(output)], env=env, check=True
        )
    assemble_fractions(output)
    print('Created patient_scores_n230.csv: 230 patients, 11 scores, 7 cell fractions.')


if __name__ == '__main__':
    main()
