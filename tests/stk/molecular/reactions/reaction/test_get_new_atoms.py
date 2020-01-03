import numpy as np
import itertools as it

from .utilities import is_equivalent_atom


def test_get_new_atoms(test_case):
    _test_get_new_atoms(
        reaction_result=test_case.reaction.get_result(),
        new_atoms=test_case.new_atoms,
    )


def _test_get_new_atoms(reaction_result, new_atoms):
    atoms = it.zip_longest(reaction_result.get_new_atoms(), new_atoms)
    for (atom1, position1), (atom2, position2) in atoms:
        is_equivalent_atom(atom1, atom2)
        assert np.allclose(position1, position1, atol=1e-32)
