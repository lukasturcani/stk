.DEFAULT_GOAL := help

.PHONY: help
help: ## Show this help
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "\033[1m%-15s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

.PHONY: check
check: ## Check the source code
	-./docs/remake_modules.bash
	-make -C docs doctest
	-pytest
	-isort .
	-mypy src
	-flake8 src tests benchmarks
