import stk


def test_has_correct_molecular_graph() -> None:
    bb = stk.BuildingBlock.from_smiles("CN")
    expected_atoms = [
        stk.Atom(0, 6),
        stk.Atom(1, 7),
        stk.Atom(2, 1),
        stk.Atom(3, 1),
        stk.Atom(4, 1),
        stk.Atom(5, 1),
        stk.Atom(6, 1),
    ]
    assert bb.atoms == expected_atoms
    assert bb.integer_bonds == [
        stk.IntegerBond(
            atom1=expected_atoms[0],
            atom2=expected_atoms[1],
            order=1,
        ),
        stk.IntegerBond(
            atom1=expected_atoms[0],
            atom2=expected_atoms[2],
            order=1,
        ),
        stk.IntegerBond(
            atom1=expected_atoms[0],
            atom2=expected_atoms[3],
            order=1,
        ),
        stk.IntegerBond(
            atom1=expected_atoms[0],
            atom2=expected_atoms[4],
            order=1,
        ),
        stk.IntegerBond(
            atom1=expected_atoms[1],
            atom2=expected_atoms[5],
            order=1,
        ),
        stk.IntegerBond(
            atom1=expected_atoms[1],
            atom2=expected_atoms[6],
            order=1,
        ),
    ]
    assert bb.dative_bonds == []


def test_has_correct_functional_groups() -> None:
    bb = stk.BuildingBlock.from_smiles("BrCCBr", stk.bromo())
    assert bb.functional_groups == [
        stk.FunctionalGroup([1], [0]),
        stk.FunctionalGroup([2], [3]),
    ]
