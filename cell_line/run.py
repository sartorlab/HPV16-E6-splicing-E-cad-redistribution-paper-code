#!/usr/bin/env python3
"""Run the UM-SCC-38 analyses or draw selected manuscript panels.

Inputs are ordinary TSV files in this module. Figure mode reads supplied DE/ORA
results; analysis modes begin with the 12-sample integer count matrix. Outputs
always go to a new directory, preserving the supplied inputs.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

MODULE = Path(__file__).resolve().parent
PANELS = ['S7' + letter for letter in 'ABCDEF'] + ['S8' + letter for letter in 'ABCDE']
VOLCANO_INPUTS = dict(zip(
    ['S7' + letter for letter in 'ABCDEF'],
    ['01_E6IvsFL', '02_E6IvsFLm', '03_E6IvsPCD',
     '05_FLvsPCD', '06_FLmvsPCD', '04_FLvsFLm'],
))


def copy_inputs(source, destination):
    """Copy readable TSV inputs or scripts into the isolated run directory."""
    shutil.copytree(source, destination, dirs_exist_ok=True)


def check_counts(output):
    """Require nonnegative integer counts in the declared 12-sample order."""
    with (output / 'inputs/gene_counts.tsv').open() as stream:
        rows = csv.reader(stream, delimiter='\t')
        header = next(rows)
        number_of_genes = 0
        for row in rows:
            if len(row) != 13 or any(int(value) < 0 for value in row[1:]):
                raise ValueError('Counts must have 12 nonnegative integer columns.')
            number_of_genes += 1
    with (output / 'config/sample_metadata.tsv').open() as stream:
        samples = list(csv.DictReader(stream, delimiter='\t'))
    if number_of_genes != 60675 or len(samples) != 12:
        raise ValueError('This study input must contain 60,675 genes and 12 samples.')
    if header[1:] != [sample['count_column'] for sample in samples]:
        raise ValueError('Count columns do not match sample_metadata.tsv.')


def run_script(rscript, script, output, environment):
    """Run one documented R stage and retain its console output on failure."""
    command = [rscript, '--vanilla', str(output / 'scripts' / script), str(output)]
    print('+ ' + ' '.join(command), flush=True)
    logfile = output / 'logs' / (Path(script).stem + '.log')
    with logfile.open('w') as stream:
        result = subprocess.run(
            command, stdout=stream, stderr=subprocess.STDOUT,
            env=environment, cwd=output,
        )
    if result.returncode:
        print(logfile.read_text()[-10000:], file=sys.stderr)
        raise SystemExit(result.returncode)


def main():
    """Prepare the requested workflow's inputs and execute its R stages."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['figures', 'de', 'analysis', 'all'], default='figures')
    parser.add_argument('--panel', choices=PANELS, help='One panel; only with --mode figures.')
    parser.add_argument('--output', required=True, type=Path, help='New or empty output directory.')
    parser.add_argument('--rscript', default=os.environ.get('RSCRIPT', 'Rscript'))
    parser.add_argument('--gene-set-dir', type=Path,
                        help='Exact complete reference directory; required by analysis/all.')
    args = parser.parse_args()
    if args.panel and args.mode != 'figures':
        parser.error('--panel is available only with --mode figures.')
    rscript = shutil.which(args.rscript)
    if not rscript:
        parser.error('Rscript executable not found: ' + args.rscript)
    rscript = str(Path(rscript).resolve())

    # Reference versions define the analysis. Newer KEGG definitions cannot
    # silently substitute for the snapshot used in the manuscript.
    if args.mode in ['analysis', 'all'] and args.gene_set_dir is None:
        parser.error('analysis/all requires --gene-set-dir; see resources/gene_sets/README.md.')
    references = []
    if args.gene_set_dir is not None:
        checksums = json.loads((MODULE / 'resources/gene_sets/frozen_complete_checksums.json').read_text())
        for name, checksum in checksums.items():
            source = args.gene_set_dir / name
            if not source.is_file() or hashlib.sha256(source.read_bytes()).hexdigest() != checksum:
                parser.error('Missing or changed reference snapshot: ' + str(source))
            references.append(source)

    output = args.output.resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        parser.error('Output must be new or empty.')
    output.mkdir(parents=True, exist_ok=True)
    for folder in ['config', 'scripts']:
        copy_inputs(MODULE / folder, output / folder)
    for folder in ['results/de', 'results/ora', 'results/mapping', 'figures/source_data', 'logs']:
        (output / folder).mkdir(parents=True, exist_ok=True)

    if args.mode == 'figures':
        # A volcano requires its own complete DE table; an ORA panel selects its
        # comparison from the shared complete ORA results.
        if args.panel is None:
            copy_inputs(MODULE / 'data', output / 'results')
        elif args.panel.startswith('S7'):
            filename = 'res_' + VOLCANO_INPUTS[args.panel] + '.tsv'
            shutil.copy2(MODULE / 'data/de' / filename, output / 'results/de' / filename)
        else:
            shutil.copy2(MODULE / 'data/ora/all_ORA_results.tsv', output / 'results/ora/all_ORA_results.tsv')
    else:
        copy_inputs(MODULE / 'inputs', output / 'inputs')
        copy_inputs(MODULE / 'resources', output / 'resources')
        for source in references:
            shutil.copy2(source, output / 'resources/gene_sets' / source.name)
        check_counts(output)

    environment = os.environ.copy()
    # Direct module calls also use the restored library; explicit R_LIBS_USER wins.
    library_pointer = MODULE.parent / '.R_library_path'
    if 'R_LIBS_USER' not in environment and library_pointer.is_file():
        library_text = library_pointer.read_text().strip()
        library = (MODULE.parent / Path(library_text).expanduser()).resolve()
        if not library_text or not library.is_dir():
            parser.error('The R library recorded in .R_library_path does not exist.')
        environment['R_LIBS_USER'] = str(library)
    environment['PYTHON'] = sys.executable
    environment['TZ'] = 'UTC'
    environment['PYTHONDONTWRITEBYTECODE'] = '1'
    environment['MPLCONFIGDIR'] = str(output / 'logs/matplotlib')
    environment['FIGURE_PANEL'] = args.panel or ''
    scripts = []
    if args.mode != 'figures':
        scripts.append('01_deseq2.R')
    if args.mode in ['analysis', 'all']:
        scripts.extend(['02_identifier_mapping.R', '03_ora.R'])
    if args.mode in ['figures', 'all']:
        if args.panel is None or args.panel.startswith('S7'):
            scripts.append('04_volcano_plots.R')
        if args.panel is None or args.panel.startswith('S8'):
            scripts.append('05_ora_plots.R')
    for script in scripts:
        run_script(rscript, script, output, environment)
    print('Finished: ' + str(output))


if __name__ == '__main__':
    main()
