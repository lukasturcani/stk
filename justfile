# List all recipes.
default:
  @just --list

# Build the docs.
docs:
  make -C docs html
  @echo Docs are in: $PWD/docs/build/html/index.html

# Install development environment.
dev:
  pip install -e '.[dev]'

# Run code checks.
check:
  #!/usr/bin/env bash

  error=0
  trap error=1 ERR

  echo
  (set -x; ruff . )

  echo
  ( set -x; black --check . )

  echo
  ( set -x; mypy src )

  echo
  ( set -x; pytest --cov=src --cov-report term-missing )

  echo
  ( set -x; make -C docs doctest )

  test $error = 0

# Auto-fix code issues.
fix:
  black .
  ruff --fix .

# Start a MongoDB instance in docker.
mongo:
  docker run -d --rm -p 27017:27017 --name mongo mongo:latest
