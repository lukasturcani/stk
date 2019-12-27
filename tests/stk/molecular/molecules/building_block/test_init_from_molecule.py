import stk
import pytest


class _TestInitFromMolecule:
    def case1():
        functional_groups = [stk.AmineFactory()]
        building_block = stk.BuildingBlock(
            smiles='NCCN',
            functional_groups=functional_groups,
        )
        return (
            building_block,
            functional_groups,
            building_block.atoms,
            building_block.bonds,
            building_block.func_groups,
        )

    def case2():
        functional_groups = [stk.BromoFactory()]
        molecule = stk.ConstructedMolecule(
            building_blocks=[
                stk.BuildingBlock('BrCCBr', functional_groups),
            ],
            topology_graph=stk.polymer.Linear('AA', 1),
        )
        expected_functional_groups = (
            stk.Bromo(
                atoms=(stk.C(1), stk.Br(2)),
                bonders=(stk.C(1), ),
                deleters=(stk.Br(2), ),
            ),
            stk.Bromo(
                atoms=(stk.C(8), stk.Br(7)),
                bonders=(stk.C(8), ),
                deleters=(stk.Br(7), ),
            ),
        )
        return (
            molecule,
            functional_groups,
            molecule.atoms,
            molecule.bonds,
            expected_functional_groups,
        )

    @pytest.mark.parametrize(
        argnames=(
            'molecule',
            'functional_groups',
            'expected_atoms',
            'expected_bonds',
            'expected_functional_groups',
        ),
        argvalues=(
            # case1(),
            # case2(),
        ),
    )
    def test(
        self,
        molecule,
        functional_groups,
        expected_atoms,
        expected_bonds,
        expected_functional_groups,
    ):
        building_block = stk.BuildingBlock.init_from_molecule(
            mol=molecule,
            functional_groups=functional_groups,
        )
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
