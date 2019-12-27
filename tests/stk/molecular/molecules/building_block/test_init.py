import numpy as np
import itertools as it
import stk
import pytest
import rdkit.Chem.AllChem as rdkit







class _TestInitFromSmiles:
    def case1():
        expected_atoms = (
            stk.N(0),
            stk.C(1),
            stk.C(2),
            stk.N(3),
            stk.H(4),
            stk.H(5),
            stk.H(6),
            stk.H(7),
            stk.H(8),
            stk.H(9),
            stk.H(10),
            stk.H(11),
        )
        expected_bonds = (
            stk.Bond(stk.N(0), stk.C(1), 1),
            stk.Bond(stk.C(1), stk.C(2), 1),
            stk.Bond(stk.C(2), stk.N(3), 1),
            stk.Bond(stk.N(0), stk.H(4), 1),
            stk.Bond(stk.N(0), stk.H(5), 1),
            stk.Bond(stk.C(1), stk.H(6), 1),
            stk.Bond(stk.C(1), stk.H(7), 1),
            stk.Bond(stk.C(2), stk.H(8), 1),
            stk.Bond(stk.C(2), stk.H(9), 1),
            stk.Bond(stk.N(3), stk.H(10), 1),
            stk.Bond(stk.N(3), stk.H(11), 1),
        )
        expected_functional_groups = (
            stk.Amine(
                atoms=(stk.N(0), stk.H(4), stk.H(5)),
                bonders=(stk.N(0),),
                deleters=(stk.H(4), stk.H(5)),
            ),
            stk.Amine(
                atoms=(stk.N(3), stk.H(10), stk.H(11)),
                bonders=(stk.N(3), ),
                deleters=(stk.H(10), stk.H(11)),
            ),
        )
        return (
            'NCCN',
            [stk.AmineFactory()],
            expected_atoms,
            expected_bonds,
            expected_functional_groups,
        )

    @pytest.mark.parametrize(
        argnames=(
            'smiles',
            'functional_groups',
            'expected_atoms',
            'expected_bonds',
            'expected_functional_groups',
        ),
        argvalues=(
            case1(),
        ),
    )
    def test(
        self,
        smiles,
        functional_groups,
        expected_atoms,
        expected_bonds,
        expected_functional_groups,
    ):
        building_block = stk.BuildingBlock(smiles, functional_groups)
        atoms = it.zip_longest(
            building_block.get_atoms(),
            expected_atoms,
        )
        for a1, a2 in atoms:
            assert is_equivalent_atom(a1, a2)

        bonds = it.zip_longest(
            building_block.get_bonds(),
            expected_bonds,
        )
        for b1, b2 in bonds:
            assert is_equivalent_bond(b1, b2)

        fgs = it.zip_longest(
            building_block.get_functional_groups(),
            expected_functional_groups,
        )
        for fg1, fg2 in fgs:
            assert is_equivalent_fg(fg1, fg2)
