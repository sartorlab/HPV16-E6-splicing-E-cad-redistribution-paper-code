The bundled GO-BP and Hallmark definitions are human MSigDB 2026.1.Hs, obtained
through msigdbr 26.1.0. Cite MSigDB and the source collections; these distributions
use [CC BY 4.0](https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp).
Version, acquisition date and download source are in `reference_sources.tsv`.

KEGG definitions came from the human PATHWAY REST snapshot last updated
2026-10-15. Direct redistribution is unresolved, so the raw snapshot is excluded.
The complete calculated ORA results are deposited under `data/ora/` and allow all
S8 panels to be redrawn. See [KEGG terms](https://www.kegg.jp/kegg/legal.html).

For full ORA recalculation, supply `gene_set_membership.tsv` and
`gene_set_metadata.tsv` matching `frozen_complete_checksums.json`, or reconstruct
them from the four authorized frozen raw files named in `kegg_raw_checksums.json`:

```sh
python3 cell_line/scripts/prepare_gene_sets.py \
  --raw-dir /path/to/authorized/kegg_snapshot --output /tmp/frozen_gene_sets
python3 cell_line/run.py --mode all --output /tmp/cell_all \
  --gene-set-dir /tmp/frozen_gene_sets
```

The optional `--download` attempts the current academic REST service. Current
definitions may have changed: the script rejects any checksum mismatch rather
than silently changing the analysis. Exact-snapshot access remains unresolved
for public users; a successful author-side rerun does not resolve that access.
