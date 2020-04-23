import itertools as it

from ....utilities import is_equivalent_atom, normalize_ids


def test_get_atoms(case_data, get_atom_ids):
    """
    Test :meth:`.Molecule.get_atoms`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        The test case. Holds the molecule to test and the correct
        atoms.

    get_atom_ids : :class:`callable`
        Takes a single parameter, `molecule`, and returns a valid
        `atom_ids` parameter for :meth:`.Molecule.get_atoms`. This
        allows the testing of different values of this parameter.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_atoms(
        molecule=case_data.molecule,
        atoms=case_data.atoms,
        get_atom_ids=get_atom_ids,
    )


def _test_get_atoms(molecule, atoms, get_atom_ids):
    """
    Test :meth:`.Molecule.get_atoms`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    atoms : :class:`tuple` of :class:`.Atom`
        The correct atoms of `molecule`.

    get_atom_ids : :class:`callable`
        Takes a single parameter, `molecule`, and returns a valid
        `atom_ids` parameter for :meth:`.Molecule.get_atoms`. This
        allows the testing of different values of this parameter.

    Returns
    -------
    None : :class:`NoneType`

    """

    atom_ids = get_atom_ids(molecule)
    if atom_ids is None:
        atom_ids = range(molecule.get_num_atoms())

    for atom_id, atom in it.zip_longest(
        normalize_ids(molecule, atom_ids),
        molecule.get_atoms(get_atom_ids(molecule)),
    ):
        is_equivalent_atom(atoms[atom_id], atom)
