import itertools as it

from tests.utilities import is_equivalent_bond


def test_get_bonds(case_data):
    """
    Test :meth:`.Molecule.get_bonds`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        The test case. Holds the molecule to test and the correct
        bonds.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_bonds(case_data.molecule, case_data.bonds)


def _test_get_bonds(molecule, bonds):
    """
    Test :meth:`.Molecule.get_bonds`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    bonds : :class:`tuple` of :class:`.Bond`
        The correct bonds of `molecule`.

    Returns
    -------
    None : :class:`NoneType`

    """

    for bond1, bond2 in it.zip_longest(molecule.get_bonds(), bonds):
        is_equivalent_bond(bond1, bond2)
