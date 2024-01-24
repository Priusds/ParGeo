format:
	black bubbles tests
	isort bubbles tests

lint:
	ruff bubbles
	mypy bubbles

tests:
	pytest tests
