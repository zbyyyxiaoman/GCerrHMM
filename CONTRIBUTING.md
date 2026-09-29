# Contributing

Issues and pull requests are welcome.

Before opening a pull request:

1. run `python tests/run_tests.py`;
2. run `bash reproduce.sh --static-check`;
3. keep generated files and private data out of the repository;
4. state the exact command and input files used for any benchmark change;
5. do not overwrite frozen artifacts; add a new timestamped output instead.

For a new metric or simulator, include:

- the input contract and expected output;
- a small deterministic smoke test;
- an explicit statement of which claim the metric can and cannot test.

