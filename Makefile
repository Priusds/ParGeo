format:
	black bubbles tests examples
	isort bubbles tests examples

lint:
	ruff bubbles examples tests
	mypy bubbles examples tests

tests:
	pytest tests
