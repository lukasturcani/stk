import itertools as it
import pytest
import stk

from .utilities import (
    are_clone_sequences,
    atom_id,
    are_same_id_sequences,
)
from .case_data import CaseData


@pytest.mark.parametrize(
    argnames='case_data',
    argvalues=(
        CaseData(
            factory=stk.RingAmineFactory(),
            molecule=stk.BuildingBlock('NCC(Br)c1c(Br)cccc1'),
            functional_groups=(
                stk.RingAmine(
                    nitrogen=stk.N(0),
                    carbon1=stk.C(1),
                    carbon2=stk.C(2),
                    carbon3=stk.C(4),
                    hydrogen1=stk.H(11),
                    hydrogen2=stk.H(12),
                    hydrogen3=stk.H(15),
                ),
            ),
        ),
    )
)
def test_get_functional_groups(case_data):
    _test_get_functional_groups(
        factory=case_data.factory,
        molecule=case_data.molecule,
        functional_groups=case_data.functional_groups,
    )


def _test_get_functional_groups(factory, molecule, functional_groups):
    fgs = it.zip_longest(
        functional_groups,
        factory.get_functional_groups(molecule),
    )
    for expected_fg, fg in fgs:
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
