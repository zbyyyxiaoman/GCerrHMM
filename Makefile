.PHONY: install test static smoke quick-demo gc-demo audit-panel figures full-check package

install:
	python -m pip install -e .

test:
	python tests/run_tests.py

static:
	bash reproduce.sh --static-check

smoke:
	bash reproduce.sh --smoke

quick-demo:
	bash reproduce.sh --quick-demo

gc-demo:
	bash reproduce.sh --gc-demo

audit-panel:
	bash reproduce.sh --audit-panel

figures:
	bash reproduce.sh --framework-figures

full-check:
	bash reproduce.sh --full-check

package:
	python scripts/package_github_release.py

