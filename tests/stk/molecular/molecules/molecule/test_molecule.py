import os
import pytest
import numpy as np
import itertools as it
import stk













def test_clone(molecule):
    clone = molecule.clone()
    atoms = it.zip_longest(clone.get_atoms(), molecule.get_atoms())
    for a1, a2 in atoms:
        assert is_equivalent_atom(a1, a2)

    bonds = it.zip_longest(clone.get_bonds(), molecule.get_bonds())
    for b1, b2 in bonds:
        assert is_equivalent_bond(b1, b2)
