from ..utilities import has_same_structure, is_clone


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
    is_clone(molecule, clone)
