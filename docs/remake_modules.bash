#!/usr/bin/env bash

# STK_DIR=$(readlink -e "$0")
STK_DIR=$(pwd)

rm -- "$STK_DIR"/docs/source/{stk.*.rst,stk.rst}
sphinx-apidoc -feEM -o "$STK_DIR/docs/source" "$STK_DIR/src/stk"
