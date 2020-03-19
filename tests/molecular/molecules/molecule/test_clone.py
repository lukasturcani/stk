from ..utilities import has_same_structure, is_equivalent_molecule


def test_clone(molecule):
    """
    Test :meth:`.Molecule.clone`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    Returns
    -------
    None : :class:`NoneType`

    """

    clone = molecule.clone()
    has_same_structure(molecule, clone)
    assert molecule is not clone
    is_equivalent_molecule(molecule, clone)
