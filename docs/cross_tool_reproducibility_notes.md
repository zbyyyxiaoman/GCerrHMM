# Cross-tool reproducibility notes

Findings that affect how the cross-tool comparison must be read. Each one was
confirmed against the raw files rather than inferred from summary tables.

## 1. PBSim3 errhmm mode writes no quality scores

PBSim3's native error-HMM mode (`--method errhmm`) emits `!` for every quality
character, i.e. **Phred 0 for every base**. Measured mean quality of the
10x simulated FASTQ (mean over all bases, Phred+33):

| Tool | mean Q | median Q | quality characters seen |
|---|---:|---:|---|
| GCerrHMM | 42.18 | 42 | `#` .. `~` |
| NanoSim | 42.16 | 40 | `"` .. `}` |
| PBSim3 sample | 42.92 | 45 | `"` .. `S` |
| **PBSim3 errhmm** | **0.00** | **0** | **only `!`** |
| badread | 28.25 | 30 | `"` .. `{` |

Consequence: any quality-aware variant caller discards the entire pileup.
`bcftools mpileup` defaults to `-Q 13`, so the PBSim3-errhmm run returns an
**empty VCF** and a nominal F1 of 0 against 5,082 truth variants. That zero is
an input-format artefact, not a measurement of simulation accuracy, so the
variant-calling cell for this tool is reported as `NA` with the reason
recorded, instead of `0.0`.

This is worth stating in the manuscript's reproducibility section: a user who
picks PBSim3's errhmm mode for a variant-calling benchmark will silently get no
calls unless they know to ignore the quality column.

## 2. The errhmm model dominates, so identity is genome-independent

PBSim3 errhmm applies one bundled model regardless of input genome. Base-level
identity therefore comes out nearly the same for two very different genomes:

* *E. coli*: 0.860524 (mismatch 0.0584, ins 0.0320, del 0.0491)
* *S. cerevisiae*: 0.860486 (mismatch 0.0589, ins 0.0316, del 0.0490)

The two values differ in the fifth decimal. This is expected for a fixed
external model and is not a copy or caching error; the underlying JSON files
and read sets differ.

## 3. Saturated precision quantises F1

When a tool makes no false positives for an SV class, precision saturates at
1.0 and `F1 = 2R / (1 + R)` depends only on the integer number of true
positives. With a fixed number of truth events per class, independent designs
can therefore return *identical* F1 values while their call sets differ
event-by-event. This was verified by MD5 comparison of the simulated reads,
the truth VCF, the BAM and the call VCF for two such designs. See
`sv_evaluation_protocol.md` for the full matching rules.
