import itertools as it

import numpy as np

from .utilities import is_equivalent_atom


def test_get_new_atoms(case_data):
    """
    Test that correct atoms are added by a :class:`.Reaction`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Includes the reaction to be tested and the atoms,
        which should be added.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_new_atoms(
        reaction_result=case_data.reaction.get_result(),
        new_atoms=case_data.new_atoms,
    )


def _test_get_new_atoms(reaction_result, new_atoms):
    """
    Test that the correct atoms are added by `reaction_result`.

    Parameters
    ----------
    reaction_result : :class:`.ReactionResult`
        The result of a :class:`.Reaction`.

    new_atoms : :class:`tuple` of :class:`.NewAtom`
        The atoms, which should be added.

    Returns
    -------
    None : :class:`NoneType`

    """

    for (atom1, position1), (atom2, position2) in it.zip_longest(
        reaction_result.get_new_atoms(),
        new_atoms,
    ):
        is_equivalent_atom(atom1, atom2)
        assert np.allclose(position1, position1, atol=1e-32)
