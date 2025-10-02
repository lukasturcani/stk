# List all commands.
default:
  @just --list

# Build docs.
docs:
  rm -rf docs/source/_autosummary
  uv run make -C docs html
  echo Docs are in $PWD/docs/build/html/index.html

# Do a dev install.
setup:
  uv sync --all-extras --dev

# Run code checks.
check:
  #!/usr/bin/env bash

  error=0
  trap error=1 ERR

  echo
  (set -x; uv run ruff check . )

  echo
  ( set -x; uv run ruff format --check . )

  echo
  ( set -x; uv run mypy src )

  echo
  ( set -x; uv run pytest --cov=src --cov-report term-missing )

  echo
  ( set -x; uv run make -C docs doctest )

  test $error = 0

# Auto-fix code issues.
fix:
  uv run ruff format .
  uv run ruff check --fix .

# Start a MongoDB instance in docker.
mongo:
  docker run -d --rm -p 27017:27017 --name mongo mongo:latest

# Build a release.
build:
  uv build
