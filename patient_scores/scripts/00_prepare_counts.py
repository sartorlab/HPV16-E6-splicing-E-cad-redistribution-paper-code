#!/usr/bin/env python3
"""Join Bailey Garb's four patient count sources before batch correction.

The HVC table determines gene-row order. The join preserves repeated symbols
and every matching row combination, as in the original inner joins. TCGA
Ensembl IDs are mapped to symbols using the supplied org.Hs.eg.db 3.17 map.
Patient columns remain in source order: HVC, FFPE-family, FF, then TCGA.
"""

import argparse
import csv
import itertools
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
GROUPS = ['HVC', 'FFPE_family', 'FF', 'TCGA']


def read_count_sources():
    """Return source rows without collapsing duplicate gene symbols."""
    with (ROOT/'references/ensembl_to_symbol_orgHs3.17.tsv').open() as handle:
        gene_map = {
            row['ensembl_id']: row['symbol']
            for row in csv.DictReader(handle, delimiter='\t')
        }

    tables = {}
    sample_order = []
    left_rows = []
    for group in GROUPS:
        data = {}
        with (ROOT/'data'/f'counts_{group}.csv').open(newline='') as handle:
            reader = csv.reader(handle)
            header = next(reader)
            sample_order.extend(header[1:])
            for row in reader:
                if len(row) != len(header):
                    raise ValueError(f'Count-row width differs from header in {group}.')
                gene = row[0]
                if group == 'TCGA':
                    # Remove an Ensembl version suffix before the frozen mapping.
                    gene = gene_map.get(re.sub(r'\.\d+$', '', gene))
                if not gene:
                    continue
                values = row[1:]
                data.setdefault(gene, []).append(values)
                if group == 'HVC':
                    left_rows.append((gene, values))
        tables[group] = data
    return tables, sample_order, left_rows


def main():
    """Write the complete joined count matrix needed for normalization."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    tables, sample_order, left_rows = read_count_sources()
    with (ROOT/'data/cohort_n230.csv').open() as handle:
        cohort = {row['samples'] for row in csv.DictReader(handle)}
    if len(sample_order) != 230 or len(set(sample_order)) != 230 or set(sample_order) != cohort:
        raise ValueError('Count sources must contain the same 230 unique cohort patients.')

    joined_rows = 0
    destination = args.output/'uncorrected_joined_counts_n230.csv'
    with destination.open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['gene']+sample_order)
        for gene, values in left_rows:
            # Empty matches remove a gene from the inner join. Multiple matches
            # produce all combinations; do not silently collapse duplicate genes.
            matching_rows = (tables[group].get(gene, []) for group in GROUPS[1:])
            for right_rows in itertools.product(*matching_rows):
                counts = values+list(itertools.chain.from_iterable(right_rows))
                writer.writerow([gene]+counts)
                joined_rows += 1
    if joined_rows == 0:
        raise ValueError('No gene symbols are shared by all four count sources.')
    print(f'Joined {joined_rows} gene rows across 230 patients.')


if __name__ == '__main__':
    main()
