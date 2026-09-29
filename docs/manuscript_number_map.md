# Manuscript number map

This file maps the manuscript figures, tables and Additional files to the
repository artifacts and the command that should be used to inspect or rebuild
each item. It is a provenance index, not a replacement for the frozen files.

## Figure map

| Manuscript item | Primary artifact | Generator / source | Reviewer command |
|---|---|---|---|
| Fig. 1 | `docs/gcerrhmm_main_figures_refined_20260927_v2/figure1_overview_innovation.*` | `scripts/build_gcerrhmm_main_figures.py` | `bash reproduce.sh --framework-figures` |
| Fig. 2 | background figure set in `docs/background_figures/` | `scripts/build_background_figures.py` | `bash reproduce.sh --framework-figures` |
| Fig. 3 | `docs/figures_supplementary/figS_dup_ambiguity.*` | supplementary figure builder | inspect frozen artifact / `--smoke` |
| Fig. 4 | `figure1_overview_innovation.*` | main-figure builder | `--framework-figures` |
| Fig. 5 | `figure2_reads_cross_tool.*` | main-figure builder | `--framework-figures` |
| Fig. 6 | `figure3_alignment.*` | main-figure builder | `--framework-figures` |
| Fig. 7 | `figure4_variant_calling.*` | main-figure builder | `--framework-figures` |
| Fig. 8 | `figure5_assembly_phasing.*` | main-figure builder | `--framework-figures` |
| Fig. 9 | `figure6_ablation.*` | main-figure builder | `--framework-figures` |

The frozen figure archives are checked by:

```bash
bash reproduce.sh --smoke
```

## Table map

| Manuscript item | Primary artifact | Source / command |
|---|---|---|
| Table 1 | species source metadata | `config/species_config_v2.json`, `config/data_sources.json` |
| Table 2 | alignment metrics | `results/framework/stats/mapping.csv`; `--framework-figures` |
| Table 3 | chr21 30x delta-to-real panel | `scripts/delta_to_real_decision.py`; `--audit-panel` |
| Table 4 | assembly metrics | `results/framework/stats/assembly.csv`; `--framework-figures` |
| Table 5 | GC-bin ablation | `results/stats/gc_bins_cross_species_with_1bin_seeds_*.csv` |
| Table 6 | depth/uncertainty analysis | coverage and depth outputs under `results/stats/` |
| Table 7 | six-species stratified panel | `results/stats/gc_fidelity_summary.csv` (when present in full tree) |
| Table 8 | simulator relationship table | `docs/final_cross_tool_comparison.md` and manuscript Table 8 |

## Additional files

| Additional file | Expected artifact | Status |
|---|---|---|
| Additional file 1 | `Additional_file_1_species_panel.xlsx` | to be generated/verified |
| Additional file 2 | `Additional_file_2_level1_sub_scores.xlsx` | to be generated/verified |
| Additional file 3 | `Additional_file_3_gc_bin_ablation.xlsx` | to be generated/verified |
| Additional file 4 | `Additional_file_4_delta_to_real.xlsx` | to be generated/verified |
| Additional file 5 | `Additional_file_5_instrument_audit.pdf` | to be generated/verified |
| Additional file 6 | `Additional_file_6_uncertainty_layers.pdf` | to be generated/verified |
| Additional file 7 | `Additional_file_7_hifi_crossplatform_note.pdf` | to be generated/verified |

## Model parameters and seeds

- MAPQ: chr21 uses `0`; the remaining paper panel uses `20`.
- Truth-variant seed: `42`.
- GC-bin ablation seeds: `11, 22, 33`.
- Coverage-replicate seeds: `101, 202, 303`.
- Quick demo seed: `42`.

The machine-readable seed map is `config/experiment_seeds.json`.
