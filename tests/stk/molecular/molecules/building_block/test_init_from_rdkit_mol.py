import stk
import rdkit.Chem.AllChem as rdkit
import pytest


class _FunctionalGroup:
    def __init__(self, atom_ids, bonder_ids, deleter_ids):
        self.atom_ids = atom_ids
        self.bonder_ids = bonder_ids
        self.deleter_ids = deleter_ids


class TestInitFromRdkitMol:
    def case1():
        rdkit_molecule = rdkit.AddHs(rdkit.MolFromSmiles('NCCN'))
        rdkit.EmbedMolecule(rdkit_molecule, rdkit.ETKDGv2())

        functional_groups = [stk.AmineFactory()]

        expected_functional_groups = [
            _FunctionalGroup(
                atom_ids=(0, 4, 5),
                bonder_ids=(0, ),
                deleter_ids=(4, 5),
            ),
            _FunctionalGroup(
                atom_ids=(3, 10, 11),
                bonder_ids=(3, ),
                deleter_ids=(10, 11),
            ),
        ]
        return (
            rdkit_molecule,
            functional_groups,
            expected_functional_groups,
        )

    def case2():
        rdkit_molecule = rdkit.AddHs(rdkit.MolFromSmiles('NCCN'))
        rdkit.EmbedMolecule(rdkit_molecule, rdkit.ETKDGv2())
        return rdkit_molecule, None, []

    def case3():
        rdkit_molecule = rdkit.AddHs(rdkit.MolFromSmiles('NCCN'))
        rdkit.EmbedMolecule(rdkit_molecule, rdkit.ETKDGv2())
        return rdkit_molecule, [stk.AldehydeFactory()], []

    @pytest.mark.parametrize(
        argnames=(
            'rdkit_molecule',
            'functional_groups',
            'expected_functional_groups',
        ),
        argvalues=[
            case1(),
            case2(),
            case3(),
        ],
    )
    def test(
        self,
        rdkit_molecule,
        functional_groups,
        expected_functional_groups,
    ):
        molecule = stk.BuildingBlock.init_from_rdkit_mol(
            mol=rdkit_molecule,
            functional_groups=functional_groups,
        )
        assert rdkit_molecule.GetNumAtoms() == len(molecule.atoms)
        atoms = zip(molecule.atoms, rdkit_molecule.GetAtoms())
        for atom, rdkit_atom in atoms:
            assert atom.charge == rdkit_atom.GetFormalCharge()
            assert atom.atomic_number == rdkit_atom.GetAtomicNum()
            assert atom.mass == rdkit_atom.GetMass()

        assert len(molecule.bonds) == rdkit_molecule.GetNumBonds()
        bonds = zip(molecule.bonds, rdkit_molecule.GetBonds())
        for bond, rdkit_bond in bonds:
            assert bond.order == rdkit_bond.GetBondTypeAsDouble()
            assert bond.atom1.id == rdkit_bond.GetBeginAtomIdx()
            assert bond.atom2.id == rdkit_bond.GetEndAtomIdx()

        num_func_groups = len(molecule.func_groups)
        assert num_func_groups == len(expected_functional_groups)
        func_groups = zip(
            molecule.func_groups,
            expected_functional_groups,
        )
        for fg, expected_fg in func_groups:
            assert tuple(fg.get_atom_ids()) == expected_fg.atom_ids
            assert tuple(fg.get_bonder_ids()) == expected_fg.bonder_ids
            assert (
                tuple(fg.get_deleter_ids()) == expected_fg.deleter_ids
            )
