# Panel-specific data preparation

The paper does not use whole-genome mouse or whole-genome human panels. Its
PANEL keys are:

| Panel | Paper source | Public quick path |
|---|---|---|
| `Ecoli` | RefSeq E. coli + SRR39619343 ONT | `bash reproduce.sh --quick-reproduce` |
| `Scerevisiae` | RefSeq S. cerevisiae + SRR39791672 ONT | same command with `SPECIES=Scerevisiae` |
| `Athaliana` | RefSeq A. thaliana + ERR5716408 ONT | SRA conversion is supported by `--prepare-demo` |
| `Dmelanogaster` | RefSeq D. melanogaster + SRR22071686 ONT | SRA conversion is supported by `--prepare-demo` |
| `Mmusculus_chr19` | GRCm39 chr19 + SRR14685232 ONT | requires explicit chr19 extraction |
| `Hsapiens_chr21` | HG002 chr21 30x ONT anchor | requires the HG002 chr21 span/anchor preparation |

## Mouse chr19

The full mouse reference/BAM must be reduced to chromosome 19 before the
quick demo. A suitable preparation is:

```bash
samtools faidx data/references/Mmusculus_ref.fa chr19 \
  > data/panels/mouse_chr19/references/Mmusculus_chr19_ref.fa

samtools view -b data/real_reads_verified/Mmusculus_ont_aligned.bam chr19 \
  | samtools sort -@ 8 -m 1G \
  > data/panels/mouse_chr19/real_reads_verified/Mmusculus_chr19_ont_aligned.bam
samtools index data/panels/mouse_chr19/real_reads_verified/Mmusculus_chr19_ont_aligned.bam

samtools fastq -@ 8 \
  data/panels/mouse_chr19/real_reads_verified/Mmusculus_chr19_ont_aligned.bam \
  | pigz -p 8 \
  > data/panels/mouse_chr19/real_reads_verified/Mmusculus_chr19_ont.fastq.gz
```

Then point the quick demo at the panel root:

```bash
PROJECT_DIR="$PWD/data/panels/mouse_chr19" \
SPECIES=Mmusculus_chr19 \
THREADS=16 JOBS=3 \
bash reproduce.sh --quick-demo
```

## Human chr21

The human paper panel is the HG002 chr21 30x real ONT anchor, not a
whole-genome ERR13491992 panel. Use the project's chr21 span/anchor script
to obtain the BAM and reference, then place the resulting files as:

```text
data/panels/human_chr21/
  data/references/Hsapiens_chr21_ref.fa
  data/real_reads_verified/Hsapiens_chr21_ont_aligned.bam
  data/real_reads_verified/Hsapiens_chr21_ont.fastq.gz
```

The FASTQ can be regenerated from the anchor BAM:

```bash
samtools fastq -@ 8 \
  data/panels/human_chr21/real_reads_verified/Hsapiens_chr21_ont_aligned.bam \
  | pigz -p 8 \
  > data/panels/human_chr21/real_reads_verified/Hsapiens_chr21_ont.fastq.gz
```

Then run:

```bash
PROJECT_DIR="$PWD/data/panels/human_chr21" \
SPECIES=Hsapiens_chr21 \
THREADS=16 JOBS=3 \
bash reproduce.sh --quick-demo
```

The frozen decision panel still requires the original coverage-matched
`*_Hsapiens_chr21_30x.json` files and is checked with:

```bash
bash reproduce.sh --audit-panel
```

## Important boundary

`prepare_demo_data.py` deliberately refuses whole-genome `Mmusculus` and
`Hsapiens` inputs. This prevents a reviewer from accidentally treating a
whole-genome run as the paper's chr19/chr21 panel.
