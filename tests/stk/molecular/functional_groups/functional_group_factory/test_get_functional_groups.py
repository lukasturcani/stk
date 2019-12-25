import itertools as it
import pytest
import stk


class _TestCase:
    """
    A test case.

    Attributes
    ----------
    factory : :class:`.FunctionalGroupFactory`
        The factory being tested.

    molecule : :class:`.Molecule`
        The molecule passed to
        :meth:`.FunctionalGroupFactory.get_functional_groups`.

    functional_groups : :class:`tuple` of :class:`.FunctionalGroup`
        The expected functional groups.

    """

    def __init__(self, factory, molecule, functional_groups):
        self.factory = factory
        self.molecule = molecule
        self.functional_groups = functional_groups


@pytest.mark.parametrize(
    argnames='test_case',
    argvalues=(
        _TestCase(
            factory=stk.AmineFactory(),
            molecule=stk.BuildingBlock('NCCN'),
            functional_groups=(
                stk.Amine(
                    atoms=(stk.N(0), stk.H(4), stk.H(5)),
                    bonders=(stk.N(0), ),
                    deleters=(stk.H(4), stk.H(5)),
                ),
                stk.Amine(
                    atoms=(stk.N(3), stk.H(10), stk.H(11)),
                    bonders=(stk.N(3), ),
                    deleters=(stk.H(10), stk.H(11)),
                ),
            ),
        ),

        _TestCase(
            factory=stk.AmineFactory(num_deleters=1),
            molecule=stk.BuildingBlock('NCCN'),
            functional_groups=(
                stk.Amine(
                    atoms=(stk.N(0), stk.H(4), stk.H(5)),
                    bonders=(stk.N(0), ),
                    deleters=(stk.H(4), ),
                ),
                stk.Amine(
                    atoms=(stk.N(3), stk.H(10), stk.H(11)),
                    bonders=(stk.N(3), ),
                    deleters=(stk.H(10), ),
                ),
            ),
        ),

        _TestCase(
            factory=stk.AmineFactory(num_deleters=0),
            molecule=stk.BuildingBlock('NCCN'),
            functional_groups=(
                stk.Amine(
                    atoms=(stk.N(0), stk.H(4), stk.H(5)),
                    bonders=(stk.N(0), ),
                    deleters=(),
                ),
                stk.Amine(
                    atoms=(stk.N(3), stk.H(10), stk.H(11)),
                    bonders=(stk.N(3), ),
                    deleters=(),
                ),
            ),
        ),

        _TestCase(
            factory=stk.AmineFactory(),
            molecule=stk.BuildingBlock('CCCC'),
            functional_groups=(),
        ),

        _TestCase(
            factory=stk.SecondaryAmineFactory(),
            molecule=stk.BuildingBlock('CNCCNCC'),
            functional_groups=(
                stk.Amine(
                    atoms=(stk.N(1), stk.C(0), stk.C(2), stk.H(10)),
                    bonders=(stk.N(1), ),
                    deleters=(stk.H(10), ),
                ),
                stk.Amine(
                    atoms=(stk.N(4), stk.C(3), stk.C(5), stk.H(15)),
                    bonders=(stk.N(4), ),
                    deleters=(stk.H(15), ),
                ),
            ),
        ),

        _TestCase(
            factory=stk.AldehydeFactory(),
            molecule=stk.BuildingBlock('O=CCC=O'),
            functional_groups=(
                stk.Aldehyde(
                    atoms=(stk.O(0), stk.C(1), stk.H(5)),
                    bonders=(stk.C(1), ),
                    deleters=(stk.O(0), ),
                ),
                stk.Aldehyde(
                    atoms=(stk.O(4), stk.C(3), stk.H(8)),
                    bonders=(stk.C(3), ),
                    deleters=(stk.O(4), ),
                ),
            ),
        ),

        _TestCase(
            factory=stk.CarboxylicAcidFactory(),
            molecule=stk.BuildingBlock('O=C(O)CCCC(O)=O'),
            functional_groups=(
                stk.CarboxylicAcid(
                    atoms=(stk.O(0), stk.C(1), stk.O(2), stk.H(9)),
                    bonders=(stk.C(1), ),
                    deleters=(stk.O(2), stk.H(9)),
                ),
                stk.CarboxylicAcid(
                    atoms=(stk.O(7), stk.O(8), stk.C(6), stk.H(16)),
                    bonders=(stk.C(6), ),
                    deleters=(stk.O(7), stk.H(16)),
                ),
            ),
        ),

        _TestCase(
            factory=stk.AmideFactory(),
            molecule=stk.BuildingBlock('O=C(N)CCC(=O)N'),
            functional_groups=(
                stk.Amide(
                    atoms=(
                        stk.O(0),
                        stk.C(1),
                        stk.N(2),
                        stk.H(8),
                        stk.H(9),
                    ),
                    bonders=(stk.C(1), ),
                    deleters=(stk.N(2), stk.H(8), stk.H(9)),
                ),
                stk.Amide(
                    atoms=(
                        stk.C(5),
                        stk.O(6),
                        stk.N(7),
                        stk.H(14),
                        stk.H(15),
                    ),
                    bonders=(stk.C(5), ),
                    deleters=(stk.N(7), stk.H(14), stk.H(15)),
                ),
            ),
        ),

        _TestCase(
            factory=stk.ThioacidFactory(),
            molecule=stk.BuildingBlock('O=C(S)CC(S)=O'),
            functional_groups=(
                stk.Thioacid(
                    atoms=(stk.O(0), stk.C(1), stk.S(2), stk.H(7)),
                    bonders=(stk.C(1), ),
                    deleters=(stk.S(2), stk.H(7)),
                ),
                stk.Thioacid(
                    atoms=(stk.C(4), stk.S(5), stk.O(6), stk.H(10)),
                    bonders=(stk.C(4), ),
                    deleters=(stk.S(5), stk.H(10)),
                ),
            ),
        ),

        _TestCase(
            factory=stk.AlcoholFactory(),
            molecule=stk.BuildingBlock('OCCCO'),
            functional_groups=(
                stk.Alcohol(
                    atoms=(stk.O(0), stk.H(5)),
                    bonders=(stk.O(0), ),
                    deleters=(stk.H(5), ),
                ),
                stk.Alcohol(
                    atoms=(stk.O(4), stk.H(12)),
                    bonders=(stk.O(4), ),
                    deleters=(stk.H(12), ),
                ),
            ),
        ),

        _TestCase(
            factory=stk.ThiolFactory(),
            molecule=stk.BuildingBlock('SCCCS'),
            functional_groups=(
                stk.Thiol(
                    atoms=(stk.S(0), stk.H(5)),
                    bonders=(stk.S(0), ),
                    deleters=(stk.H(5), ),
                ),
                stk.Thiol(
                    atoms=(stk.S(4), stk.H(12)),
                    bonders=(stk.S(4), ),
                    deleters=(stk.H(12), ),
                ),
            ),
        ),

        _TestCase(
            factory=stk.FluoroFactory(),
            molecule=stk.BuildingBlock('FCC(F)CCF'),
            functional_groups=(
                stk.Fluoro(
                    atoms=(stk.F(0), stk.C(1)),
                    bonders=(stk.C(1), ),
                    deleters=(stk.F(0), ),
                ),
                stk.Fluoro(
                    atoms=(stk.F(3), stk.C(2)),
                    bonders=(stk.C(2), ),
                    deleters=(stk.F(3), ),
                ),
                stk.Fluoro(
                    atoms=(stk.F(6), stk.C(5)),
                    bonders=(stk.C(5), ),
                    deleters=(stk.F(6), ),
                ),
            ),
        ),

        _TestCase(
            factory=stk.BromoFactory(),
            molecule=stk.BuildingBlock('BrCC(Br)CCBr'),
            functional_groups=(
                stk.Bromo(
                    atoms=(stk.Br(0), stk.C(1)),
                    bonders=(stk.C(1), ),
                    deleters=(stk.Br(0), ),
                ),
                stk.Bromo(
                    atoms=(stk.Br(3), stk.C(2)),
                    bonders=(stk.C(2), ),
                    deleters=(stk.Br(3), ),
                ),
                stk.Bromo(
                    atoms=(stk.Br(6), stk.C(5)),
                    bonders=(stk.C(5), ),
                    deleters=(stk.Br(6), ),
                ),
            ),
        ),

        _TestCase(
            factory=stk.IodoFactory(),
            molecule=stk.BuildingBlock('ICC(I)CCI'),
            functional_groups=(
                stk.Iodo(
                    atoms=(stk.I(0), stk.C(1)),
                    bonders=(stk.C(1), ),
                    deleters=(stk.I(0), ),
                ),
                stk.Iodo(
                    atoms=(stk.I(3), stk.C(2)),
                    bonders=(stk.C(2), ),
                    deleters=(stk.I(3), ),
                ),
                stk.Iodo(
                    atoms=(stk.I(6), stk.C(5)),
                    bonders=(stk.C(5), ),
                    deleters=(stk.I(6), ),
                ),
            ),
        ),

        _TestCase(
            factory=stk.TerminalAlkyneFactory(),
            molecule=stk.BuildingBlock('C#CC#CC'),
            functional_groups=(
                stk.TerminalAlkyne(
                    atoms=(stk.C(0), stk.C(1), stk.H(5)),
                    bonders=(stk.C(0), ),
                    deleters=(stk.H(5), ),
                ),
            ),
        ),

        _TestCase(
            factory=stk.TerminalAlkyneFactory(delete_carbon=True),
            molecule=stk.BuildingBlock('C#CC#CC'),
            functional_groups=(
                stk.TerminalAlkyne(
                    atoms=(stk.C(0), stk.C(1), stk.H(5)),
                    bonders=(stk.C(1), ),
                    deleters=(stk.H(5), stk.C(0)),
                ),
            ),
        ),

        _TestCase(
            factory=stk.TerminalAlkeneFactory(),
            molecule=stk.BuildingBlock('C=CC=CC'),
            functional_groups=(
                stk.TerminalAlkene(
                    atoms=(stk.C(0), stk.C(1), stk.H(5), stk.H(6)),
                    bonders=(stk.C(1), ),
                    deleters=(stk.C(0), stk.H(5), stk.H(6)),
                ),
            ),
        ),

        _TestCase(
            factory=stk.BoronicAcidFactory(),
            molecule=stk.BuildingBlock('B(O)(O)CCB(O)O'),
            functional_groups=(
                stk.BoronicAcid(
                    atoms=(
                        stk.B(0),
                        stk.O(1),
                        stk.O(2),
                        stk.H(8),
                        stk.H(9),
                    ),
                    bonders=(stk.B(0), ),
                    deleters=(stk.O(1), stk.O(2), stk.H(8), stk.H(9)),
                ),
                stk.BoronicAcid(
                    atoms=(
                        stk.B(5),
                        stk.O(6),
                        stk.O(7),
                        stk.H(14),
                        stk.H(15)
                    ),
                    bonders=(stk.B(5), ),
                    deleters=(
                        stk.O(6),
                        stk.O(7),
                        stk.H(14),
                        stk.H(15),
                    ),
                ),
            ),
        ),

        _TestCase(
            factory=stk.DiolFactory(),
            molecule=stk.BuildingBlock('CC(O)C(O)CC'),
            functional_groups=(
                stk.Diol(
                    atoms=(
                        stk.C(1),
                        stk.O(2),
                        stk.C(3),
                        stk.O(4),
                        stk.H(11),
                        stk.H(13),
                    ),
                    bonders=(stk.O(2), stk.O(4)),
                    deleters=(stk.H(11), stk.H(13)),
                ),
            ),
        ),

        _TestCase(
            factory=stk.DifluoroFactory(),
            molecule=stk.BuildingBlock('CC(F)C(F)CC'),
            functional_groups=(
                stk.Difluoro(
                    atoms=(stk.C(1), stk.F(2), stk.C(3), stk.F(4)),
                    bonders=(stk.C(1), stk.C(3)),
                    deleters=(stk.F(2), stk.F(4)),
                ),
            ),
        ),

        _TestCase(
            factory=stk.DibromoFactory(),
            molecule=stk.BuildingBlock('CC(Br)C(Br)CC'),
            functional_groups=(
                stk.Dibromo(
                    atoms=(stk.C(1), stk.Br(2), stk.C(3), stk.Br(4)),
                    bonders=(stk.C(1), stk.C(3)),
                    deleters=(stk.Br(2), stk.Br(4)),
                ),
            ),
        ),

        _TestCase(
            factory=stk.RingAmineFactory(),
            molecule=stk.BuildingBlock('NCC(Br)c1c(Br)cccc1'),
            functional_groups=(
                stk.RingAmine(
                    atoms=(
                        stk.N(0),
                        stk.C(1),
                        stk.C(2),
                        stk.C(4),
                        stk.H(11),
                        stk.H(12),
                        stk.H(15),
                    ),
                    bonders=(stk.N(0), stk.C(2)),
                    deleters=(stk.H(11), stk.H(12), stk.H(15)),
                ),
            ),
        ),
    )
)
def test_get_functional_groups(test_case):
    _test_get_functional_groups(
        factory=test_case.factory,
        molecule=test_case.molecule,
        functional_groups=test_case.functional_groups,
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
    are_clone_sequences(
        atoms1=sorted(functional_group1.get_bonders(), key=atom_id),
        atoms2=sorted(functional_group2.get_bonders(), key=atom_id),
    )
    are_clone_sequences(
        atoms1=sorted(functional_group1.get_deleters(), key=atom_id),
        atoms2=sorted(functional_group2.get_deleters(), key=atom_id),
    )


def atom_id(atom):
    return atom.id


def are_clone_sequences(atoms1, atoms2):
    """
    Test if `atoms1` and `atoms2` are clones of each other.

    """

    for a1, a2 in it.zip_longest(atoms1, atoms2):
        assert a1 is not a2
        assert a1.id == a2.id
        assert a1.charge == a2.charge
        assert a1.__class__ is a2.__class__
