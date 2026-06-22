.PHONY: install install-dev install-editable format lint type-check test test-cov clean help

# Default target
help:
	@echo "PEvolution Development Commands"
	@echo "================================"
	@echo "make install           - Install production dependencies (uses pipenv or pip depending on environment)"
	@echo "make install-dev       - Install all dependencies including dev tools"
	@echo "make install-editable  - Install package in editable mode (pip install -e .)"
	@echo "make format            - Format code with ruff"
	@echo "make lint              - Check code quality with ruff"
	@echo "make type-check        - Run type checking with mypy"
	@echo "make test              - Run unit tests"
	@echo "make test-cov          - Run tests with coverage report"
	@echo "make clean             - Remove build and cache files"
	@echo "make help              - Show this help message"

install:
	pipenv install || pip install -r requirements.txt

install-dev:
	pipenv install --dev || (pip install -r requirements-dev.txt)

install-editable:
	pip install -e .

format:
	pipenv run ruff format --line-length 120 . || ruff format --line-length 120 .

lint:
	pipenv run ruff check . || ruff check .

type-check:
	pipenv run mypy Utils/ helpers/ || mypy Utils/ helpers/

test:
	pipenv run pytest test/ -v

test-cov:
	pipenv run pytest test/ -v --cov=Utils --cov=helpers --cov-report=html --cov-report=term || pytest test/ -v --cov=Utils --cov=helpers --cov-report=html --cov-report=term

clean:
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name .pytest_cache -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name .mypy_cache -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name .ruff_cache -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name htmlcov -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name Orthos -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name XML -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name .DS_Store -delete 2>/dev/null || true
	@echo "Clean complete"
