## Summary of available make targets:
##
## make help         -- Display this message
## make -B venv      -- (Re)install all requisite Python packages
## make format       -- Run code formatter
## make lint         -- Run linter and type checker
## make tests        -- Run tests
## make licenses     -- Summarise licenses used by packages in your venv
##
## See also 'make -C requirements help' for information about the
## requirements files.
##
## This Makefile needs to be run inside a virtual environment.

ifndef VIRTUAL_ENV
ifdef CONDA_PREFIX
$(warning "For better compatibility, consider using a plain Python venv instead of Conda")
VIRTUAL_ENV := $(CONDA_PREFIX)
else
$(error "This Makefile needs to be run inside a virtual environment")
endif
endif

# The list of requirements/*.txt files to install in the venv. If you
# have a requirements/extra.txt, it is always installed.
REQUIREMENTS := base ci dev

help:
	@sed -nE 's/^## ?//p' $(MAKEFILE_LIST)

venv: $(VIRTUAL_ENV)/timestamp

$(VIRTUAL_ENV)/timestamp:
	pip install pip-tools
	pip-sync requirements/requirements.txt
	pip install -e .
	@touch $(VIRTUAL_ENV)/timestamp

format: venv
	isort bubbles tests
	black bubbles tests

lint: venv
	ruff bubbles tests
	mypy bubbles tests

tests: venv
	pytest -m '$(MARK)'

licenses: venv
	pip-licenses --from=mixed --order=license --summary

.PHONY: help venv format lint tests licenses
