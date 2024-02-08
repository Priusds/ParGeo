## This Makefile needs to be run inside an active virtual environment.
## If you don't have one, you can create one with:
##	python3 -m venv .venv
##	source .venv/bin/activate

ifndef VIRTUAL_ENV
$(error "This Makefile needs to be run inside a Python virtual environment.")
endif

# Install bubbles using poetry
install:
	poetry install

# Install bubbles using pip
pip_install:
	pip install -r requirements.txt
	pip install -e .

# Update requirements.txt
update_requirements:
	poetry export -f requirements.txt --output requirements.txt --with dev --without-hashes

# Run code formatter
format:
	black bubbles tests examples
	isort bubbles tests examples

#Run linter and type checker
lint:
	ruff bubbles tests/test_topology.py
	mypy bubbles tests/test_topology.py

# Run tests
tests:
	pytest tests

# Build documentation
docs:
	sphinx-apidoc -o docs bubbles && cd docs && make html

.PHONY: format lint tests docs
