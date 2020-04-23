import itertools as it

from ...utilities import is_clone_generic_functional_group


def test_with_atoms(generic_functional_group, get_atom_map):
    """
    Test :meth:`.GenericFunctionalGroup.with_atoms`.

    Parameters
    ----------
    generic_functional_group : :class:`.GenericFunctionalGroup`
        The functional group to test.

    get_atom_map : :class:`callable`
        Takes a single parameter, `functional_group`, and returns a
        valid `atom_map` parameter for its
        :meth:`.GenericFunctionalGroup.with_atoms` method.
        This allows the testing of different values of this parameter.

    Returns
    -------
    None : :class:`NoneType`

    """

    # Save a clone to ensure that "generic_functional_group" is not
    # changed by the test, because it meant to be immutable.
    before = generic_functional_group.clone()
    _test_with_atoms(generic_functional_group, get_atom_map)
    is_clone_generic_functional_group(before, generic_functional_group)


def _test_with_atoms(functional_group, get_atom_map):
    """
    Test :meth:`.GenericFunctionalGroup.with_atoms`.

    Parameters
    ----------
    generic_functional_group : :class:`.GenericFunctionalGroup`
        The functional group to test.

    get_atom_map : :class:`callable`
        Takes a single parameter, `functional_group`, and returns a
        valid `atom_map` parameter for its
        :meth:`.GenericFunctionalGroup.with_atoms` method.
        This allows the testing of different values of this parameter.

    Returns
    -------
    None : :class:`NoneType`

    """

    atom_map = get_atom_map(functional_group)
    clone = functional_group.with_atoms(atom_map)

    is_modified_sequence(
        atoms1=functional_group.get_atoms(),
        atoms2=clone.get_atoms(),
        atom_map=atom_map,
    )
    is_modified_id_sequence(
        ids1=functional_group.get_atom_ids(),
        ids2=clone.get_atom_ids(),
        atom_map=atom_map,
    )

    is_modified_sequence(
        atoms1=functional_group.get_bonders(),
        atoms2=clone.get_bonders(),
        atom_map=atom_map,
    )
    is_modified_id_sequence(
        ids1=functional_group.get_bonder_ids(),
        ids2=clone.get_bonder_ids(),
        atom_map=atom_map,
    )

    is_modified_sequence(
        atoms1=functional_group.get_deleters(),
        atoms2=clone.get_deleters(),
        atom_map=atom_map,
    )
    is_modified_id_sequence(
        ids1=functional_group.get_deleter_ids(),
        ids2=clone.get_deleter_ids(),
        atom_map=atom_map,
    )

    is_modified_id_sequence(
        ids1=functional_group.get_placer_ids(),
        ids2=clone.get_placer_ids(),
        atom_map=atom_map,
    )
    is_modified_id_sequence(
        ids1=functional_group.get_core_atom_ids(),
        ids2=clone.get_core_atom_ids(),
        atom_map=atom_map,
    )


def is_modified_sequence(atoms1, atoms2, atom_map):
    for atom1, atom2 in it.zip_longest(atoms1, atoms2):
        assert atom2 is atom_map.get(atom1.get_id(), atom1)


def is_modified_id_sequence(ids1, ids2, atom_map):
    for id1, id2 in it.zip_longest(ids1, ids2):
        if id1 in atom_map:
            assert id2 == atom_map[id1].get_id()
        else:
            assert id2 == id1
