## Summary of available make targets:
##
## make install      -- Install bubbles in the current virtual environment (using poetry)
## make pip_install  -- Install bubbles in the current virtual environment (using pip)
## make format       -- Run code formatter
## make lint         -- Run linter and type checker
## make tests        -- Run tests.
## make docs		 -- Build documentation.
##
## See also 'make -C requirements help' for information about the
## requirements files.
##
## This Makefile needs to be run inside an active virtual environment.
## If you don't have one, you can create one with:
##	python3 -m venv .venv
##	source .venv/bin/activate

ifndef VIRTUAL_ENV
$(error "This Makefile needs to be run inside a Python virtual environment.")
endif

install:
	poetry install

pip_install:
	pip install -r requirements.txt
	pip install -e .

format:
	black bubbles tests examples
	isort bubbles tests examples

lint:
	ruff bubbles tests/test_topology.py
	mypy bubbles tests/test_topology.py

tests:
	pytest tests

docs:
	cd docs && make html

.PHONY: format lint tests docs