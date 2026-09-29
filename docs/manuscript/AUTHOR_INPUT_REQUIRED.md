# Author-only fields required before submission

The W3 submission-format manuscript is complete except for facts that must
come from the authors. Do not invent these values. Fill
`docs/manuscript/author_metadata.template.json` and then run the metadata
application script.

Author names, affiliations and the corresponding-author designation are now
supplied for Boyang Zhang and Tao Jiang. Email addresses, ORCIDs, funding,
competing interests, CRediT contributions and the Zenodo DOI still require
confirmation before submission.

## Required fields

| Field | Where it is used | Status |
|---|---|---|
| Full author list, order and ORCID | title page and contributions | required |
| Affiliations and department/city/country | title page | required |
| Corresponding author name and email | title page and cover letter | required |
| Funding agency and grant number(s), or "none" | Declarations | required |
| Competing interests | Declarations and cover letter | required |
| Author contributions in CRediT format | Declarations | required |
| Acknowledgements, or "none" | Declarations | required |
| Repository URL | Availability of data and materials | required |
| Zenodo DOI | Availability of data and materials | required |
| CC BY or CC BY-NC-ND licence choice | submission system | required |
| Confirmation that all authors approved submission | cover letter | required |
| Confirmation that the manuscript is not published or under consideration elsewhere | cover letter | required |
| Optional suggested reviewers | cover letter | optional |

## Fill-in workflow

1. Copy `author_metadata.template.json` to `author_metadata.json`.
2. Fill every required value and set both confirmation flags to `true`.
3. Run:

```bash
python3 scripts/apply_author_metadata.py \
  --w3-docx /path/to/W3_submission.docx \
  --metadata /path/to/author_metadata.json \
  --output-docx /path/to/W3_submission_final.docx \
  --output-cover-letter /path/to/cover_letter.docx
```

The script fails loudly if a required value is missing, so an incomplete
metadata file cannot silently produce a submission-ready document.
