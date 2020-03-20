import itertools as it


def test_clone(constructed_molecule):
    """
    Test :meth:`.ConstructedMolecule.clone`.

    Parameters
    ----------
    constructed_molecule : :class:`.ConstructedMolecule`
        The constructed molecule to test.

    Returns
    -------
    None : :class:`NoneType`

    """

    clone = constructed_molecule.clone()
    is_clone(constructed_molecule, clone)


def is_clone(molecule1, molecule2):
    for atom_info1, atom_info2 in it.zip_longest(
        molecule1.get_atom_infos(),
        molecule2.get_atom_infos(),
    ):
        assert atom_info1 is atom_info2
    for bond_info1, bond_info2 in it.zip_longest(
        molecule1.get_bond_infos(),
        molecule2.get_bond_infos(),
    ):
        assert bond_info1 is bond_info2
