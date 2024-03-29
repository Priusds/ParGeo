## This Makefile needs to be run inside an active virtual environment.
## If you don't have one, you can create one with:
##	python3 -m venv .venv
##	source .venv/bin/activate

ifndef VIRTUAL_ENV
$(error "This Makefile needs to be run inside a Python virtual environment.")
endif

.PHONY: format lint tests docs clean install pip_install update_requirements

# Install pargeo using poetry
install:
	poetry install --with dev

# Install pargeo using pip
pip_install:
	pip install -r requirements.txt
	pip install -e .

# Update requirements.txt
update_requirements:
	poetry export -f requirements.txt --output requirements.txt --with dev --with docs --without-hashes

# Run code formatter
format:
	ruff format pargeo tests examples

#Run linter and type checker
lint:
	ruff pargeo tests examples
	mypy pargeo tests examples
# Run tests
tests:
	pytest tests

# Build documentation
docs:
	sphinx-apidoc -o docs pargeo && cd docs && make html
html:
	rm -r docs/_build && cd docs && make html

# Delete all .GEO and .MSH files
clean:
	find . -name "*.geo_unrolled" -type f -delete
	find . -name "*.msh" -type f -delete


