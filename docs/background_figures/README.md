# Background figure provenance

The W2 Background section uses three generated figure pairs (PNG and PDF).
All values are either taken from frozen project tables or shown as a
qualitative literature synthesis; none of the background panels is used as a
performance ranking.

## B1. GC context

* Panel a: `docs/paper_figures/table1_species_panel.csv`.
* Panel b: `docs/paper_innovation/stats/gc_error_curve.csv`.
* Purpose: show that the six-genome panel spans different GC contexts and
  that real and simulated error curves can differ in shape.

## B2. Simulator landscape

* Panel a: publication years and primary modelling emphasis for NanoSim [1],
  NPBSS [3], PBSim3 [2], Meta-NanoSim [9], Squigulator [10], Icarust [11],
  Platinum Pedigree [13], context-aware simulation [14] and NanoSimFormer
  [15].
* Panel b: qualitative feature coverage based on the primary descriptions in
  those papers. `yes`, `partial` and `-` are design labels, not benchmark
  scores.

## B3. Context layers

* Panel a: `docs/paper_innovation/stats/homopolymer_deletion.csv`, real
  E. coli ONT reads only. The log-scaled y axis shows the context dependence
  of deletion rate across reference homopolymer lengths.
* Panel b: evidence-layer synthesis using the 2026 ONT simulator benchmark
  [16] together with the signal, benchmark and downstream papers cited in
  the W2 reference list. The 2026 benchmark is explicitly a preprint.

## Build

From the code root, run:

```bash
nice -n 19 python3 scripts/build_background_figures.py \
  --project-dir /path/to/errhmm_project \
  --output-dir /path/to/output
```

The builder writes deterministic PNG and PDF pairs. The W2 manuscript builder
expects them under `background_figures/` in the package directory.
The handoff package also carries the three input tables under
`background_data/` so the panels can be checked without the full result tree.
