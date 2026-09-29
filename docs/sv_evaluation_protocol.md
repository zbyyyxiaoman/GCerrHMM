# Structural-variant evaluation protocol

This note fixes the scoring rules used for the R4 structural-variant (SV)
layer. It exists because the deletion sub-score reaches F1 = 1.0000 for two
tools, and a perfect score invites the question of whether the matching rule
is too permissive. The rule below is deliberately explicit.

## 1. Truth set

The public truth panels (GIAB v4.2.1 for human chr21, SGD-derived panels for
yeast, NCBI-derived panel for *E. coli*) contain **no variants >= 50 bp**, so
SV calling cannot be scored on them. Instead we use a controlled spike-in:

| Property | Value |
|---|---|
| Event types | DEL, INS, INV, DUP |
| Sizes per event | DEL 200 / 1k / 5k / 20k; INS 200 / 1k / 5k; INV 500 / 2k / 10k; DUP 500 / 1k / 3k |
| Events per (type, size) | 8 |
| Total events | **104** (DEL 32, INS 24, INV 24, DUP 24) |
| Designs (seeds) | 3 (`--seed 11 / 22 / 33`) |
| Placement | deterministic, evenly spread with jitter; injected high-to-low coordinate order |
| Truth coordinates | reported in the **original** reference frame, i.e. the frame reads are mapped to |

The donor reference is built from the project's `*_variant_ref.fa`, so each
simulated read set carries both the small-variant truth and the injected SVs.
Every tool - internal routes and external tools - is simulated from the same
donor for a given species and design.

Sanity checks run on every build:

* donor length minus variant-reference length equals the analytically
  expected net delta (for *E. coli* design seed 11: `-124000`, exact match);
* the injected tandem duplications are re-discoverable in the donor by a
  sequence-signature scan (`s[p:p+200] == s[p+size:p+size+200]`).

## 2. Calling

Identical command for every simulator and every design:

```
minimap2 -ax map-ont --secondary=no -t <threads> <original_ref> <reads>
  | samtools sort
sniffles -i aligned.bam -v sv_calls.vcf --threads <threads> --minsvlen 50 --mapq 20
```

* caller: Sniffles2
* minimum SV length: 50 bp
* minimum mapping quality: 20
* `BND`, `TRA` and `CNV` records are excluded from scoring.

## 3. Matching rule

One-to-one greedy matching, per contig and per SV type, best score first.
Type must be identical; a `DEL` call never matches a `DUP` truth event.

| SV type | Match condition |
|---|---|
| DEL, INV, DUP | reciprocal overlap >= 0.5 of the **shorter** interval |
| INS | \|ΔPOS\| <= 300 bp **and** \|ΔSVLEN\| <= max(100 bp, 50 % of truth SVLEN) |

Reciprocal overlap is `intersection / min(len_truth, len_call)`, so a short
call that sits entirely inside a long truth interval does not score unless it
covers at least half of the truth interval as well.

Reported metrics: precision, recall and F1 per SV type and overall, plus
Wilson 95 % intervals on precision and recall. Multiple spike-in designs are
reported as mean ± SD, and pooled TP/FP/FN with Wilson intervals are given for
the combined design set. Single-design point estimates are not used in the
manuscript.

## 4. Why a perfect deletion score is expected

Deletions in this panel are 200 bp - 20 kb, i.e. far above the 50 bp calling
floor, and their breakpoints are unambiguous in long reads. With
`precision = 1.0` (no spurious deletion calls) the score is driven purely by
recall, `F1 = 2R / (1 + R)`.

The same quantisation explains an otherwise surprising observation: two
independent designs can return *identical* per-type F1 values even though
their call sets differ line by line. Because precision is saturated at 1.0
and the recall denominator is the fixed number of truth events per type
(DEL 32, INS 24, INV 24, DUP 24), identical F1 values only require the same
integer TP count, not the same calls. We verified this directly: for PBSim3
designs 11 and 22 the simulated FASTQ, truth VCF, BAM and call VCF all have
different MD5 hashes and every call record differs, while both designs
recover DEL 32/32, DUP 7/24 and INV 16/24.

## 5. Known limitation

Tandem duplications are the weak class for **every** tool in this panel
(recall ~0.29 at 10x). Reads originating from the two identical copies map
ambiguously between them, which cancels the depth signature a caller would
otherwise use. Depth measured back against the donor reference over the
duplicated interval is 0.45-1.50x the genome mean rather than the ~2x a clean
duplication would give. This is a caller-side limitation of the benchmark
design, not a property of any individual simulator: all four tools are
reference-based samplers and draw reads from the whole donor sequence.
