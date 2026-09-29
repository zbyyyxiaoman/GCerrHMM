# GC-aware simulator algorithm

## Purpose

The implementation learns an error-state model from an aligned long-read
BAM and applies it to a target reference. The GC-aware model and the
one-bin control are trained from exactly the same reads; the only
intentional difference is whether transition probabilities are conditioned
on local GC content.

## State model

The default state space is:

`M`, `S`, `I`, `D1`, `D2`, `D3`, `D4+`

where:

- `M` is a reference-matching aligned base;
- `S` is a substitution;
- `I` is an insertion;
- `D1` to `D3` are short deletions;
- `D4+` is a deletion of four or more bases.

For a configured number of GC bins `K`, each 100-bp reference window is
assigned to one bin from local GC content. The training pass estimates
`P(state_{i+1} | state_i, GC_bin)` and a nucleotide substitution matrix
from primary alignments passing the MAPQ threshold.

The one-bin control sets `K = 1`, which removes the GC conditioning while
keeps the rest of the state space and training data unchanged.

## Simulation

For each simulated read:

1. sample a start position and read length from the requested platform or
   the supplied real-read profile;
2. walk along the target reference;
3. look up the GC bin at the current reference position;
4. sample the next error state from the corresponding transition row;
5. emit match, substitution, insertion, or deletion bases until the read
   length target is reached;
6. emit FASTQ sequence and quality values.

The target reference can be an ordinary reference or a truth-variant
reference. The latter is used for downstream variant-calling evaluation so
that all simulators see the same underlying genome.

## Reproducible improvement check

`scripts/reproduce_gc_improvement.py` trains both models from one BAM,
validates that the GC-aware model contains different transition matrices
across bins, simulates both routes plus the empirical baseline, and
evaluates them against the same real FASTQ.

The output is a C-minus-B delta table. The gate records positive or
negative direction for composite, GC, and k-mer scores. A negative result
is preserved; it is not converted into a positive claim.

## Reuse outside this repository

The algorithm is not tied to the six species panel. Supply:

- `--bam`: an aligned long-read BAM;
- `--reference`: the reference used for the BAM alignment;
- `--simulation-reference`: the sequence to simulate from;
- `--real-fastq`: the real reads used for the matched evaluation profile;
- `--species`: a label for output provenance.

The model is written as plain JSON and can be loaded by
`src/generate_simulated.py` without retraining.

## Parallel execution

`generate_simulated.py --threads N` now uses process-level parallelism for
CPU-bound Python simulation. The parent samples the read layout once and
derives a deterministic seed for every read. Workers generate fixed-size
read chunks, and results are written in read-index order.

This gives two properties that matter for a real simulator:

1. changing the worker count does not change the generated FASTQ;
2. process-level parallelism avoids the Python GIL and can use multiple
   cores without changing the model or the output distribution.

`scripts/benchmark_parallel_simulation.py` records wall time, reads/s,
speedup versus one worker, and read-count/base-length equivalence fields.
The worker ceiling is capped at 16 inside the generator, matching the
shared-server discipline.

## HiFi quality-value audit

The public E. coli HiFi FASTQ (`ERR12723508`) was checked directly:
raw ASCII quality characters have mean `113.28`; subtracting Phred+33
gives mean `80.28`, and `pysam.FastxRecord.get_quality_array()` gives
mean `80.08` on the same sample. Therefore the reported mean QV
`83.45` is not an off-by-33 script bug. The run itself contains many
`~` (Phred 93) bases, so its absolute QV distribution is
encoding-limited and should not be treated as a typical PacBio HiFi
quality benchmark. HiFi claims should rely on length, GC, and k-mer
results unless a quality-complete source BAM is used.
