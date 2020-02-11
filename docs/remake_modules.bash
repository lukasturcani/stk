STK_DIR=$(readlink -e "$0")
STK_DIR=${STK_DIR%/*/*}

rm -- "$STK_DIR/docs/source/{stk.*.rst,stk.rst}"
sphinx-apidoc -fe -o "$STK_DIR/docs/source" "$STK_DIR/src/stk"
