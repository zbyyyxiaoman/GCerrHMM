# Human HiFi region strategy

Date: 2026-09-22

## Question

Can the human HiFi critical path be shortened by fetching only chromosome
21 reads instead of the full public run?

## Findings

1. The ENA BAM for `ERR13110527` is unaligned.
   `samtools view` reports `FLAG=4` and the header has no `@SQ` records.
   Region-based remote BAM slicing is therefore impossible.
2. The ENA FASTQ is a single gzip stream, not a BGZF file with indexed
   virtual offsets. A prefix of the stream is not a chromosome subset;
   it only makes a shallower genome-wide subsample.
3. `fastq-dump -N/-X` offers spot-range extraction, but a bounded test
   against the remote SRA access layer timed out after 300 seconds.
   It remains a possible fallback if NCBI remote access recovers.
4. A usable aligned route was found in the public GIAB S3 mirror:
   `HG002.Sequel.15kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam`.
   Remote `samtools view <S3-BAM> 21 | samtools fastq` is streaming
   chr21 reads at a usable rate. This avoids the 26-32 GB FASTQ
   download and provides a real diploid HG002 HiFi read set for the
   phasing and HiFi panels. The extracted reads are re-aligned to the
   project GRCh38 chr21 reference before use.

## Decision

Use the S3 aligned-BAM extraction as the primary human HiFi path. Keep the
old partial ENA file untouched as provenance. Do not claim a human HiFi
panel until the extracted 21 reads have completed re-alignment, training,
simulation, assembly, and phasing.
