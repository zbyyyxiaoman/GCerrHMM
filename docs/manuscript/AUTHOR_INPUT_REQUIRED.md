# Author-only fields required before submission

The W3 submission-format manuscript is complete except for facts that must
come from the authors. Do not invent these values. Fill
`docs/manuscript/author_metadata.template.json` and then run the metadata
application script.

Author names, affiliations and the corresponding-author designation are now
supplied for Boyang Zhang and Tao Jiang. The manuscript now uses the default
statements "no competing interests", "funding not applicable", "not
applicable" for acknowledgements, and CRediT contributions based on the
roles supplied by the authors. Email addresses and ORCIDs remain
submission-system fields.

## Required fields

| Field | Where it is used | Status |
|---|---|---|
| Full author list and order | title page and contributions | supplied |
| Author ORCIDs | title page | required if available |
| Affiliations and department/city/country | title page | supplied |
| Corresponding author name and email | title page and cover letter | name/affiliation supplied; email required |
| Funding agency and grant number(s), or "none" | Declarations | supplied as "Not applicable"; confirm if funding exists |
| Competing interests | Declarations and cover letter | supplied as "The authors declare that they have no competing interests" |
| Author contributions in CRediT format | Declarations | supplied from the roles provided by the authors |
| Acknowledgements, or "none" | Declarations | supplied as "Not applicable" |
| Repository URL | Availability of data and materials | supplied |
| Zenodo DOI | Availability of data and materials | supplied |
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
