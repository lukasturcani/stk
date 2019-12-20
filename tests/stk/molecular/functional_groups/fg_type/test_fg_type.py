import itertools as it
import pytest
from collections import Counter
import stk


def is_equivalent_fg(fg1, fg2):
    atoms1 = Counter(
        (a.__class__, a.id) for a in fg1.atoms
    )
    atoms2 = Counter(
        (a.__class__, a.id) for a in fg2.atoms
    )

    bonders1 = Counter(
        (a.__class__, a.id) for a in fg1.bonders
    )
    bonders2 = Counter(
        (a.__class__, a.id) for a in fg2.bonders
    )

    deleters1 = Counter(
        (a.__class__, a.id) for a in fg1.deleters
    )
    deleters2 = Counter(
        (a.__class__, a.id) for a in fg2.deleters
    )

    return (
        atoms1 == atoms2
        and bonders1 == bonders2
        and deleters1 == deleters2
    )


class TestGetFunctionalGroups:
    def case1():
        amine = stk.functional_groups.fg_types['amine']
        molecule = stk.BuildingBlock('NCCN')
        functional_groups = {
            (0, 4, 5): stk.FunctionalGroup(
                atoms=(stk.N(0), stk.H(4), stk.H(5)),
                bonders=(stk.N(0), ),
                deleters=(stk.H(4), stk.H(5)),
                fg_type=amine,
            ),
            (3, 10, 11): stk.FunctionalGroup(
                atoms=(stk.N(3), stk.H(10), stk.H(11)),
                bonders=(stk.N(3), ),
                deleters=(stk.H(10), stk.H(11)),
                fg_type=amine,
            ),
        }
        return amine, molecule, functional_groups

    @pytest.mark.parametrize(
        argnames=(
            'fg_type',
            'molecule',
            'expected_functional_groups',
        ),
        argvalues=(
            case1(),
        )
    )
    def test(
        self,
        fg_type,
        molecule,
        expected_functional_groups,
    ):
        for fg in fg_type.get_functional_groups(molecule):
            key = tuple(sorted(fg.get_atom_ids()))
            other = expected_functional_groups[key]
            assert is_equivalent_fg(fg, other)
