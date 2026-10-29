#!/usr/bin/env python3
"""Assemble exact pathway definitions from supplied MSigDB and authorized KEGG data.

Inputs: bundled GO-BP/Hallmark TSVs and the four KEGG source files listed in
kegg_raw_checksums.json. Outputs: complete membership and metadata TSVs in a new
reference directory. Matching checksums enforce the reference version used in
this study. No statistical analysis or figure rendering occurs in this step.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import urllib.request

MODULE = Path(__file__).resolve().parents[1]
RESOURCES = MODULE / 'resources/gene_sets'
URLS = {
    'kegg_human_pathways.txt': 'https://rest.kegg.jp/list/pathway/hsa',
    'kegg_human_membership.txt': 'https://rest.kegg.jp/link/hsa/pathway',
    'kegg_pathway_info.txt': 'https://rest.kegg.jp/info/pathway',
    'kegg_human_info.txt': 'https://rest.kegg.jp/info/hsa',
}


def read_table(name):
    """Read a bundled TSV while retaining identifier strings and source row order."""
    with (RESOURCES / name).open() as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        return reader.fieldnames, list(reader)


def require_snapshot(directory, checksums):
    """Reject missing or changed references, including a newer KEGG release."""
    for name, checksum in checksums.items():
        source = directory / name
        if not source.is_file() or hashlib.sha256(source.read_bytes()).hexdigest() != checksum:
            raise SystemExit(
                'Missing or changed reference snapshot: ' + name +
                '. Obtain the exact authorized snapshot used for this paper.'
            )


def main():
    """Merge the two reference sources without changing their gene memberships."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-dir', required=True, type=Path,
                        help='Directory containing the four authorized KEGG source files.')
    parser.add_argument('--output', required=True, type=Path,
                        help='New directory for complete reference definitions.')
    parser.add_argument('--download', action='store_true',
                        help='Attempt academic REST download; changed definitions are rejected.')
    args = parser.parse_args()

    # Downloads are optional because a live service may no longer supply the
    # exact reference version. Existing source files are never overwritten.
    if args.download:
        args.raw_dir.mkdir(parents=True, exist_ok=True)
        for name, url in URLS.items():
            destination = args.raw_dir / name
            if not destination.exists():
                with urllib.request.urlopen(url, timeout=120) as response:
                    destination.write_bytes(response.read())
    require_snapshot(args.raw_dir, json.loads((RESOURCES / 'kegg_raw_checksums.json').read_text()))
    if args.output.exists() and (not args.output.is_dir() or any(args.output.iterdir())):
        raise SystemExit('Output must be new or empty.')
    args.output.mkdir(parents=True, exist_ok=True)

    fields, members = read_table('gene_set_membership.tsv')
    metadata_fields, metadata = read_table('gene_set_metadata.tsv')
    with (args.raw_dir / 'kegg_human_pathways.txt').open() as stream:
        pathways = {
            row[0].removeprefix('path:'): row[1]
            for row in csv.reader(stream, delimiter='\t')
        }
    with (args.raw_dir / 'kegg_human_membership.txt').open() as stream:
        for pathway, gene in csv.reader(stream, delimiter='\t'):
            pathway = pathway.removeprefix('path:')
            if pathway in pathways:
                members.append({
                    'collection': 'KEGG_PATHWAY', 'set_id': pathway,
                    'entrez_id': gene.removeprefix('hsa:'),
                })
    release = '; '.join(
        line.strip() for line in (args.raw_dir / 'kegg_pathway_info.txt').read_text().splitlines()
        if re.search('Release|Last update', line, re.I)
    )
    for pathway, name in pathways.items():
        metadata.append({
            'collection': 'KEGG_PATHWAY', 'set_id': pathway,
            'set_name': re.sub(r' - Homo sapiens \(human\)$', '', name),
            'source_id': pathway, 'original_source': 'KEGG human PATHWAY',
            'description': name, 'source_url': 'https://www.kegg.jp/pathway/' + pathway,
            'version': release,
        })

    # Preserve MSigDB source collation; order KEGG's numeric identifiers explicitly.
    # Quoting and line endings are fixed so exact reference checksums are portable.
    tables = [
        ('gene_set_membership.tsv', fields, members, ['collection', 'set_id', 'entrez_id']),
        ('gene_set_metadata.tsv', metadata_fields, metadata, ['collection', 'set_id']),
    ]
    for name, columns, rows, keys in tables:
        msigdb_rows = [row for row in rows if row['collection'] != 'KEGG_PATHWAY']
        kegg_rows = sorted(
            [row for row in rows if row['collection'] == 'KEGG_PATHWAY'],
            key=lambda row: tuple(row[key] for key in keys),
        )
        with (args.output / name).open('w') as stream:
            writer = csv.DictWriter(
                stream, columns, delimiter='\t', quoting=csv.QUOTE_ALL, lineterminator='\n',
            )
            writer.writeheader()
            writer.writerows(msigdb_rows + kegg_rows)
    require_snapshot(args.output, json.loads((RESOURCES / 'frozen_complete_checksums.json').read_text()))
    print('Complete reference definitions: ' + str(args.output))


if __name__ == '__main__':
    main()
