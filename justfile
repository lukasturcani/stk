# Check the source code.
check:
    #!/usr/bin/env bash
    set -uxo pipefail
    err=0
    trap 'err=1' ERR

    ./docs/remake_modules.bash
    make -C docs doctest
    pytest
    isort .
    mypy src
    flake8 .
    black .

    test $err = 0
