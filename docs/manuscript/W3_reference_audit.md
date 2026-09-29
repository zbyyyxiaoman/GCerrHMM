# Reference audit for W3.3

## Verification status

All bibliographic references were checked against Crossref on 2026-09-29.
Badread is a software repository rather than a DOI-bearing article; its URL
was checked directly and the access date is recorded in the reference list.

| No. | Reference | Verified identifier | Citation role |
|---:|---|---|---|
| 1 | NanoSim | 10.1093/gigascience/gix010 | Read-level simulator baseline |
| 2 | PBSim3 | 10.1093/nargab/lqac092 | HMM/error-profile simulator baseline |
| 3 | NPBSS | 10.1186/s12859-018-2208-0 | Empirical long-read simulator comparator |
| 4 | Meta-NanoSim | 10.1093/gigascience/giad013 | Metagenomic simulation extension |
| 5 | Squigulator | 10.1101/gr.278730.123 | Tunable signal-level simulation |
| 6 | Icarust | 10.1093/bioinformatics/btae141 | Adaptive-sampling simulation |
| 7 | Transformer signal simulator | 10.1093/bioinformatics/btae744 | Feed-forward transformer signal model |
| 8 | NanoSimFormer | 10.1093/bioinformatics/btag402 | Current transformer/basecaller-guided simulator |
| 9 | Context-aware simulation | 10.1093/gigascience/giag079 | Context-aware parameter optimisation |
| 10 | Platinum Pedigree | 10.1038/s41592-025-02750-y | Long-read benchmark reference |
| 11 | Taouk et al. | 10.64898/2026.05.06.723380 | 2026 ONT simulator benchmark preprint |
| 12 | Badread | GitHub software URL | Read-level simulator comparator |
| 13 | Minimap2 | 10.1093/bioinformatics/bty191 | Alignment engine |
| 14 | SAMtools | 10.1093/bioinformatics/btp352 | Alignment processing |
| 15 | Flye | 10.1038/s41587-019-0072-8 | ONT assembly |
| 16 | WhatsHap | 10.1089/cmb.2014.0157 | Diploid phasing |

## Citation-order correction

The previous W3.3 order had two problems:

1. references 4-8 were present in the list but not explicitly cited in the
   main text;
2. the list order did not follow first mention.

The revision renumbers all references by first mention:

* Background: 1-11;
* Results, reads-level comparison: 12 (Badread);
* Results, alignment: 13-14 (Minimap2, SAMtools);
* Results, assembly/phasing: 15-16 (Flye, WhatsHap).

All 16 references are now cited, and the reference list follows the same
order.

## Citation-mark formatting

BMC Bioinformatics uses numbered references and inline bracketed citations.
The manuscript therefore keeps inline forms such as `[1]`, `[1-3]` and
`[13, 14]` rather than converting them to superscripts. The numbers,
group ranges and first-mention order were checked after renumbering.
