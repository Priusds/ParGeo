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