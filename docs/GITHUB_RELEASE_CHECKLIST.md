# GitHub release checklist

Before publishing the release repository:

1. Replace `OWNER` in the README badge URLs with the GitHub account or
   organisation.
2. Fill `repository-code`, `authors` and version metadata in `CITATION.cff`.
3. Confirm the MIT licence or replace it with the institution's approved
   licence before upload.
4. Confirm that no files under `data/`, `results/`, `logs/` or local caches
   are tracked.
5. Run:

```bash
PYTHONPATH=src python tests/run_tests.py
PYTHONPATH=src bash reproduce.sh --static-check
PYTHONPATH=src bash reproduce.sh --smoke
PYTHONPATH=src bash reproduce.sh --quick-demo  # requires the documented data layout
```

6. Rebuild the release directory and archive:

```bash
python scripts/package_github_release.py \
  --source . \
  --output dist/gcerrhmm_github_20260929 \
  --zip dist/GCerrHMM_github_20260929.zip
```

The packager starts from an empty staging directory, rejects private
usernames/hostnames/local paths in text and binary artifacts, and writes a
new `RELEASE_MANIFEST.sha256`.

7. Verify the clean package:

```bash
cd dist/gcerrhmm_github_20260929
PYTHONPATH=src python tests/run_tests.py
PYTHONPATH=src bash reproduce.sh --static-check
PYTHONPATH=src bash reproduce.sh --smoke
PATH_AUDIT_STRICT=1 bash reproduce.sh --audit-paths
```

8. Enable the `.github/workflows/smoke.yml` workflow.
9. Create a `v1.0.0` tag after the first successful workflow run.
10. Upload the generated archive to Zenodo and add the DOI to `CITATION.cff`
   and the manuscript availability statement.
11. Review the manuscript author-only fields in
   `docs/manuscript/AUTHOR_INPUT_REQUIRED.md`.
12. Keep frozen results append-only; add new timestamps instead of
    overwriting existing publication artifacts.
