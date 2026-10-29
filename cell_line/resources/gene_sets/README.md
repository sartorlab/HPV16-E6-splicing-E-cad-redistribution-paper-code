# Gene-set references

The bundled GO-BP and Hallmark definitions are human MSigDB 2026.1.Hs, obtained
through msigdbr 26.1.0. Cite MSigDB and the source collections; these distributions
use [CC BY 4.0](https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp).
Version, acquisition date and download source are in
[`reference_sources.tsv`](reference_sources.tsv).

All S8 panels can be drawn from the supplied `cell_line/data/ora/` results.
Full ORA recalculation uses the human KEGG PATHWAY REST snapshot last updated
2026-10-15, subject to [KEGG terms](https://www.kegg.jp/kegg/legal.html).

Supply `gene_set_membership.tsv` and `gene_set_metadata.tsv` matching
[`frozen_complete_checksums.json`](frozen_complete_checksums.json), or assemble
them from these four KEGG source files matching
[`kegg_raw_checksums.json`](kegg_raw_checksums.json):

- `kegg_human_pathways.txt`
- `kegg_human_membership.txt`
- `kegg_pathway_info.txt`
- `kegg_human_info.txt`

Run from the repository root:

```sh
python3 cell_line/scripts/prepare_gene_sets.py \
  --raw-dir /path/to/kegg_snapshot --output /tmp/frozen_gene_sets
python3 cell_line/run.py --mode all --output /tmp/cell_all \
  --gene-set-dir /tmp/frozen_gene_sets
```

The preparation script combines KEGG with the bundled GO-BP/Hallmark definitions
and checks the resulting files against the complete-reference checksums.
Both preparation and analysis require the exact file checksums above.
