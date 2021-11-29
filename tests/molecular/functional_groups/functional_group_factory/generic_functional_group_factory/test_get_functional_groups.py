import itertools as it

from ..utilities import (
    are_clone_sequences,
    are_same_id_sequences,
    atom_id,
)


def test_get_functional_groups(case_data):
    """
    Test :meth:`.GenericFunctionalGroupFactory.get_functional_groups`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        The test case. Holds the factory, molecule and correct
        functional groups.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_functional_groups(
        factory=case_data.factory,
        molecule=case_data.molecule,
        functional_groups=case_data.functional_groups,
    )


def _test_get_functional_groups(factory, molecule, functional_groups):
    """
    Test :meth:`.GenericFunctionalGroupFactory.get_functional_groups`.

    Parameters
    ----------
    factory : :class:`.GenericFunctionalGroupFactory`
        The factory to test.

    molecule : :class:`.Molecule`
        The molecule to test.

    functional_groups : :class:`tuple`
        The correct :class:`.GenericFunctionalGroup` instances.

    Returns
    -------
    None : :class:`NoneType`

    """

    for expected_fg, fg in it.zip_longest(
        functional_groups,
        factory.get_functional_groups(molecule),
    ):
        are_clone_functional_groups(expected_fg, fg)


def are_clone_functional_groups(functional_group1, functional_group2):
    """
    Test if the functional groups are clones of each other.

    """

    assert functional_group1.__class__ is functional_group2.__class__

    are_clone_sequences(
        atoms1=sorted(functional_group1.get_atoms(), key=atom_id),
        atoms2=sorted(functional_group2.get_atoms(), key=atom_id),
    )
    are_same_id_sequences(
        ids1=sorted(functional_group1.get_atom_ids()),
        ids2=sorted(functional_group2.get_atom_ids()),
    )

    are_same_id_sequences(
        ids1=sorted(functional_group1.get_placer_ids()),
        ids2=sorted(functional_group2.get_placer_ids()),
    )

    are_same_id_sequences(
        ids1=sorted(functional_group1.get_core_atom_ids()),
        ids2=sorted(functional_group2.get_core_atom_ids()),
    )

    are_clone_sequences(
        atoms1=sorted(functional_group1.get_bonders(), key=atom_id),
        atoms2=sorted(functional_group2.get_bonders(), key=atom_id),
    )
    are_same_id_sequences(
        ids1=sorted(functional_group1.get_bonder_ids()),
        ids2=sorted(functional_group2.get_bonder_ids()),
    )

    are_clone_sequences(
        atoms1=sorted(functional_group1.get_deleters(), key=atom_id),
        atoms2=sorted(functional_group2.get_deleters(), key=atom_id),
    )
    are_same_id_sequences(
        ids1=sorted(functional_group1.get_deleter_ids()),
        ids2=sorted(functional_group2.get_deleter_ids()),
    )
