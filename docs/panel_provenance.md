# Panel provenance contract

The chr21 delta-to-real decision has one authoritative input panel:

```
results/framework/tables/mapping_<tool>_Hsapiens_chr21_30x.json
results/framework/tables/variant_<tool>_Hsapiens_chr21_30x.json
results/framework/tables/assembly_flye_<tool>_Hsapiens_chr21_30x.json
```

These are produced by `scripts/titan_chr21_30x_panel.sh`. The real anchor is:

```
results/framework/tables/mapping_real_ont30_Hsapiens_chr21.json
results/framework/tables/variant_real_ont30_Hsapiens_chr21.json
results/framework/tables/assembly_flye_real_ont30_Hsapiens_chr21.json
```

The public reviewer package ships these authoritative JSON files in
`docs/reproducibility/results_bundle/framework/tables/`. `reproduce.sh`
materialises that bundle into `results/` automatically for `--audit-panel`,
`--framework-figures`, and `--delta-to-real`.

The broader framework tables also contain chr21 `_alignment.json`,
`_fair_v2.json`, and `_fair_v3.json` files. Those belong to the main
cross-tool panel and are **not** the coverage-matched 30x decision panel.
They must never be substituted for the `*_30x.json` files merely because an
exported CSV is newer or easier to consume.

Required checks before a decision or manuscript update:

1. Run `scripts/audit_panel_provenance.py --project-dir <project>`.
2. Require `PANEL_PROVENANCE_OK`.
3. Require the decision output to say `top-2 in **2 of 2**`.
4. Require exported `mapping.csv`, `variant.csv`, and `assembly.csv` chr21
   values to equal the 30x JSON values; the exporter now aborts rather than
   silently falling back.
5. Never edit the decision conclusion to match a figure whose source files
   are not the 30x files. Fix the figure/export source selection instead.
