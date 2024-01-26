format:
	black bubbles tests examples
	isort bubbles tests examples

lint:
	ruff bubbles tests/test_topology.py
	mypy bubbles tests/test_topology.py

tests:
	pytest tests
