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
            factory=stk.PrimaryAminoFactory(),
            molecule=stk.BuildingBlock('NCCN'),
            functional_groups=(
                stk.PrimaryAmino(
                    nitrogen=stk.N(0),
                    hydrogen1=stk.H(4),
                    hydrogen2=stk.H(5),
                    atom=stk.C(1),
                    bonders=(stk.N(0), ),
                    deleters=(stk.H(4), stk.H(5)),
                ),
                stk.PrimaryAmino(
                    nitrogen=stk.N(3),
                    hydrogen1=stk.H(10),
                    hydrogen2=stk.H(11),
                    atom=stk.C(2),
                    bonders=(stk.N(3), ),
                    deleters=(stk.H(10), stk.H(11)),
                ),
            ),
        ),

        _TestCase(
            factory=stk.PrimaryAminoFactory(),
            molecule=stk.BuildingBlock('CCCC'),
            functional_groups=(),
        ),

        _TestCase(
            factory=stk.SecondaryAminoFactory(),
            molecule=stk.BuildingBlock('CNCCNCC'),
            functional_groups=(
                stk.SecondaryAmino(
                    nitrogen=stk.N(1),
                    hydrogen=stk.H(10),
                    atom1=stk.C(0),
                    atom2=stk.C(2),
                    bonders=(stk.N(1), ),
                    deleters=(stk.H(10), ),
                ),
                stk.SecondaryAmino(
                    nitrogen=stk.N(4),
                    hydrogen=stk.H(15),
                    atom1=stk.C(3),
                    atom2=stk.C(5),
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
                    carbon=stk.C(1),
                    oxygen=stk.O(0),
                    hydrogen=stk.H(5),
                    atom=stk.C(2),
                    bonders=(stk.C(1), ),
                    deleters=(stk.O(0), ),
                ),
                stk.Aldehyde(
                    carbon=stk.C(3),
                    oxygen=stk.O(4),
                    hydrogen=stk.H(8),
                    atom=stk.C(2),
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
                    carbon=stk.C(1),
                    oxygen1=stk.O(0),
                    oxygen2=stk.O(2),
                    hydrogen=stk.H(9),
                    atom=stk.C(3),
                    bonders=(stk.C(1), ),
                    deleters=(stk.O(2), stk.H(9)),
                ),
                stk.CarboxylicAcid(
                    carbon=stk.C(6),
                    oxygen1=stk.O(8),
                    oxygen2=stk.O(7),
                    hydrogen=stk.H(16),
                    atom=stk.C(5),
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
                    oxygen=stk.O(0),
                    carbon=stk.C(1),
                    nitrogen=stk.N(2),
                    hydrogen1=stk.H(8),
                    hydrogen2=stk.H(9),
                    atom=stk.C(3),
                    bonders=(stk.C(1), ),
                    deleters=(stk.N(2), stk.H(8), stk.H(9)),
                ),
                stk.Amide(
                    carbon=stk.C(5),
                    oxygen=stk.O(6),
                    nitrogen=stk.N(7),
                    hydrogen1=stk.H(14),
                    hydrogen2=stk.H(15),
                    atom=stk.C(4),
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
                    carbon=stk.C(1),
                    oxygen=stk.O(0),
                    sulfur=stk.S(2),
                    hydrogen=stk.H(7),
                    atom=stk.C(3),
                    bonders=(stk.C(1), ),
                    deleters=(stk.S(2), stk.H(7)),
                ),
                stk.Thioacid(
                    carbon=stk.C(4),
                    oxygen=stk.O(6),
                    sulfur=stk.S(5),
                    hydrogen=stk.H(10),
                    atom=stk.C(3),
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
                    oxygen=stk.O(0),
                    hydrogen=stk.H(5),
                    atom=stk.C(1),
                    bonders=(stk.O(0), ),
                    deleters=(stk.H(5), ),
                ),
                stk.Alcohol(
                    oxygen=stk.O(4),
                    hydrogen=stk.H(12),
                    atom=stk.C(3),
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
                    sulfur=stk.S(0),
                    hydrogen=stk.H(5),
                    atom=stk.C(1),
                    bonders=(stk.S(0), ),
                    deleters=(stk.H(5), ),
                ),
                stk.Thiol(
                    sulfur=stk.S(4),
                    hydrogen=stk.H(12),
                    atom=stk.C(3),
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
                    fluorine=stk.F(0),
                    atom=stk.C(1),
                    bonders=(stk.C(1), ),
                    deleters=(stk.F(0), ),
                ),
                stk.Fluoro(
                    fluorine=stk.F(3),
                    atom=stk.C(2),
                    bonders=(stk.C(2), ),
                    deleters=(stk.F(3), ),
                ),
                stk.Fluoro(
                    fluorine=stk.F(6),
                    atom=stk.C(5),
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
                    bromine=stk.Br(0),
                    atom=stk.C(1),
                    bonders=(stk.C(1), ),
                    deleters=(stk.Br(0), ),
                ),
                stk.Bromo(
                    bromine=stk.Br(3),
                    atom=stk.C(2),
                    bonders=(stk.C(2), ),
                    deleters=(stk.Br(3), ),
                ),
                stk.Bromo(
                    bromine=stk.Br(6),
                    atom=stk.C(5),
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
                    iodine=stk.I(0),
                    atom=stk.C(1),
                    bonders=(stk.C(1), ),
                    deleters=(stk.I(0), ),
                ),
                stk.Iodo(
                    iodine=stk.I(3),
                    atom=stk.C(2),
                    bonders=(stk.C(2), ),
                    deleters=(stk.I(3), ),
                ),
                stk.Iodo(
                    iodine=stk.I(6),
                    atom=stk.C(5),
                    bonders=(stk.C(5), ),
                    deleters=(stk.I(6), ),
                ),
            ),
        ),

        _TestCase(
            factory=stk.TerminalAlkyneFactory(),
            molecule=stk.BuildingBlock('C#CC#CC'),
            functional_groups=(
                stk.Alkyne(
                    carbon1=stk.C(0),
                    carbon2=stk.C(1),
                    atom1=stk.H(5),
                    atom2=stk.C(2),
                    bonders=(stk.C(1), ),
                    deleters=(stk.C(0), stk.H(5), ),
                ),
            ),
        ),

        _TestCase(
            factory=stk.TerminalAlkeneFactory(),
            molecule=stk.BuildingBlock('C=CC=CC'),
            functional_groups=(
                stk.Alkene(
                    carbon1=stk.C(0),
                    carbon2=stk.C(1),
                    atom1=stk.H(5),
                    atom2=stk.H(6),
                    atom3=stk.C(2),
                    atom4=stk.H(7),
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
                    boron=stk.B(0),
                    oxygen1=stk.O(1),
                    oxygen2=stk.O(2),
                    hydrogen1=stk.H(8),
                    hydrogen2=stk.H(9),
                    atom=stk.C(3),
                    bonders=(stk.B(0), ),
                    deleters=(stk.O(1), stk.O(2), stk.H(8), stk.H(9)),
                ),
                stk.BoronicAcid(
                    boron=stk.B(5),
                    oxygen1=stk.O(6),
                    oxygen2=stk.O(7),
                    hydrogen1=stk.H(14),
                    hydrogen2=stk.H(15),
                    atom=stk.C(4),
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
                    atom1=stk.C(1),
                    oxygen1=stk.O(2),
                    atom2=stk.C(3),
                    oxygen2=stk.O(4),
                    hydrogen1=stk.H(11),
                    hydrogen2=stk.H(13),
                    bonders=(stk.C(1), stk.C(3)),
                    deleters=(
                        stk.O(2), stk.O(4), stk.H(11), stk.H(13)
                    ),
                ),
            ),
        ),

        _TestCase(
            factory=stk.DifluoroFactory(),
            molecule=stk.BuildingBlock('CC(F)C(F)CC'),
            functional_groups=(
                stk.Difluoro(
                    atom1=stk.C(1),
                    fluorine1=stk.F(2),
                    atom2=stk.C(3),
                    fluorine2=stk.F(4),
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
                    atom1=stk.C(1),
                    bromine1=stk.Br(2),
                    atom2=stk.C(3),
                    bromine2=stk.Br(4),
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
    return atom.get_id()


def are_clone_sequences(atoms1, atoms2):
    """
    Test if `atoms1` and `atoms2` are clones of each other.

    """

    for a1, a2 in it.zip_longest(atoms1, atoms2):
        assert a1 is not a2
        assert a1.get_id() == a2.get_id()
        assert a1.get_charge() == a2.get_charge()
        assert a1.__class__ is a2.__class__
