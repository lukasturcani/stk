import itertools as it
import stk
import pytest
import rdkit.Chem.AllChem as rdkit

from ..utilities import is_equivalent_atom, normalize_ids


@pytest.mark.parametrize(
    argnames=('molecule', 'atoms'),
    argvalues=(
        (
            stk.BuildingBlock('NCCN'),
            (
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
            ),
        ),
    ),
)
def test_get_atoms_1(molecule, atoms, get_atom_ids):
    atom_ids = get_atom_ids(molecule)
    if atom_ids is None:
        atom_ids = range(molecule.get_num_atoms())

    results = it.zip_longest(
        normalize_ids(molecule, atom_ids),
        molecule.get_atoms(get_atom_ids(molecule)),
    )
    for atom_id, atom in results:
        is_equivalent_atom(atoms[atom_id], atom)


@pytest.mark.slow
def test_get_atoms_2(tmp_path, case_data):
    # Visual inspection of the molecule can be useful.
    case_data.molecule.write(tmp_path / 'molecule.mol')
    _test_get_atoms_2(case_data.molecule, case_data.smiles)


def _test_get_atoms_2(molecule, smiles):
    # The basic idea is that you get a molecule and the canonical
    # smiles for that molecule. You then create a clone "molecule"
    # making sure that atoms have canonical ordering. This conversion
    # will only work correctly when get_atoms() and get_bonds() work
    # correctly - see get_rdkit_mol(). Finally, you compare that atoms
    # of the canonically ordered molecule to an rdkit molecule which
    # also has canonical ordering, and all the atoms should match.
    molecule = with_canonical_atom_ordering(molecule)
    expected = rdkit.AddHs(rdkit.MolFromSmiles(smiles))
    for atom1, atom2 in it.zip_longest(
        molecule.get_atoms(),
        expected.GetAtoms(),
    ):
        is_equivalent_rdkit_atom(atom1, atom2)


def with_canonical_atom_ordering(molecule):
    smiles = rdkit.MolToSmiles(get_rdkit_mol(molecule))
    # This print statement is useful for checking what went wrong
    # if the actual and expected smiles don't match.
    print(smiles.replace('([H])', '').replace('[H]', ''))
    return stk.BuildingBlock(smiles)


def get_rdkit_mol(molecule):
    rdkit_mol = rdkit.EditableMol(rdkit.Mol())
    for atom in molecule.get_atoms():
        rdkit_atom = rdkit.Atom(atom.get_atomic_number())
        rdkit_atom.SetFormalCharge(atom.get_charge())
        rdkit_mol.AddAtom(rdkit_atom)

    for bond in molecule.get_bonds():
        rdkit_mol.AddBond(
            beginAtomIdx=bond.get_atom1().get_id(),
            endAtomIdx=bond.get_atom2().get_id(),
            order=rdkit.BondType(bond.get_order())
        )
    return rdkit_mol.GetMol()


def is_equivalent_rdkit_atom(atom, rdkit_atom):
    assert atom.get_atomic_number() == rdkit_atom.GetAtomicNum()
    assert atom.get_id() == rdkit_atom.GetIdx()
    assert atom.get_charge() == rdkit_atom.GetFormalCharge()
