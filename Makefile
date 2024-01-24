lint:
	ruff bubbles tests
	mypy bubbles tests

tests:
	pytest tests
