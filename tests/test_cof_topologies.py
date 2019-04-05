import stk
import os
from os.path import join
import rdkit.Chem.AllChem as rdkit

test_dir = 'cof_topology_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def test_honeycomb(amine2, aldehyde3):
    cof = stk.Periodic([amine2, aldehyde3], stk.Honeycomb())
    path = join(test_dir, 'honeycomb.sdf')
    cof.write(path)
    island = cof.island([3, 3, 1])
    rdkit.MolToMolFile(island, path.replace('.sdf', '_island.sdf'))

    assert cof.bonds_made == 4
    assert (cof.mol.GetNumAtoms() ==
            amine2.mol.GetNumAtoms()*3 +
            aldehyde3.mol.GetNumAtoms()*2 -
            cof.bonds_made*3)
    assert (cof.mol.GetNumBonds() ==
            amine2.mol.GetNumBonds()*3 +
            aldehyde3.mol.GetNumBonds()*2 -
            cof.bonds_made*2)
    assert cof.bb_counter[amine2] == 3
    assert cof.bb_counter[aldehyde3] == 2
    assert cof.topology == stk.Honeycomb()


def test_hexagonal(amine2, aldehyde6):
    cof = stk.Periodic([amine2, aldehyde6], stk.Hexagonal())
    path = os.path.join(test_dir, 'hexagonal.sdf')
    cof.write(path)
    island = cof.island([3, 3, 1])
    rdkit.MolToMolFile(island, path.replace('.sdf', '_island.sdf'))

    assert cof.bonds_made == 17
    assert (cof.mol.GetNumAtoms() ==
            amine2.mol.GetNumAtoms()*12 +
            aldehyde6.mol.GetNumAtoms()*4 -
            cof.bonds_made*3)
    assert (cof.mol.GetNumBonds() ==
            amine2.mol.GetNumBonds()*12 +
            aldehyde6.mol.GetNumBonds()*4 -
            cof.bonds_made*2)
    assert cof.bb_counter[amine2] == 12
    assert cof.bb_counter[aldehyde6] == 4
    assert cof.topology == stk.Hexagonal()


def test_square(amine2, aldehyde4):
    cof = stk.Periodic([amine2, aldehyde4], stk.Square())
    path = os.path.join(test_dir, 'square.sdf')
    cof.write(path)
    island = cof.island([3, 3, 1])
    rdkit.MolToMolFile(island, path.replace('.sdf', '_island.sdf'))

    assert cof.bonds_made == 2
    assert (cof.mol.GetNumAtoms() ==
            amine2.mol.GetNumAtoms()*2 +
            aldehyde4.mol.GetNumAtoms()*1 -
            cof.bonds_made*3)
    assert (cof.mol.GetNumBonds() ==
            amine2.mol.GetNumBonds()*2 +
            aldehyde4.mol.GetNumBonds()*1 -
            cof.bonds_made*2)
    assert cof.bb_counter[amine2] == 2
    assert cof.bb_counter[aldehyde4] == 1
    assert cof.topology == stk.Square()


def test_kagome(amine2, aldehyde4):
    cof = stk.Periodic([amine2, aldehyde4], stk.Kagome())
    path = os.path.join(test_dir, 'kagome.sdf')
    cof.write(path)
    island = cof.island([3, 3, 1])
    rdkit.MolToMolFile(island, path.replace('.sdf', '_island.sdf'))

    assert cof.bonds_made == 9
    assert (cof.mol.GetNumAtoms() ==
            amine2.mol.GetNumAtoms()*6 +
            aldehyde4.mol.GetNumAtoms()*3 -
            cof.bonds_made*3)
    assert (cof.mol.GetNumBonds() ==
            amine2.mol.GetNumBonds()*6 +
            aldehyde4.mol.GetNumBonds()*3 -
            cof.bonds_made*2)
    assert cof.bb_counter[amine2] == 6
    assert cof.bb_counter[aldehyde4] == 3
    assert cof.topology == stk.Kagome()


def test_boron_cof(diol2, boronic_acid4):
    cof = stk.Periodic([diol2, boronic_acid4], stk.Square())
    path = join(test_dir, 'boron.sdf')
    cof.write(path)
    island = cof.island([3, 3, 1])
    rdkit.MolToMolFile(island, path.replace('.sdf', '_island.sdf'))

    assert cof.bonds_made == 4
    assert (cof.mol.GetNumAtoms() ==
            diol2.mol.GetNumAtoms()*2 +
            boronic_acid4.mol.GetNumAtoms()*1 -
            cof.bonds_made*3)
    assert (cof.mol.GetNumBonds() ==
            diol2.mol.GetNumBonds()*2 +
            boronic_acid4.mol.GetNumBonds()*1 -
            cof.bonds_made*2)
    assert cof.bb_counter[diol2] == 2
    assert cof.bb_counter[boronic_acid4] == 1
    assert cof.topology == stk.Square()


def test_nolinkerhoneycomb(amine3, aldehyde3):
    cof = stk.Periodic([amine3, aldehyde3], stk.NoLinkerHoneycomb())
    path = join(test_dir, 'nolinkerhoneycomb.sdf')
    cof.write(path)
    island = cof.island([3, 3, 1])
    rdkit.MolToMolFile(island, path.replace('.sdf', '_island.sdf'))

    assert cof.bonds_made == 1
    assert (cof.mol.GetNumAtoms() ==
            amine3.mol.GetNumAtoms() +
            aldehyde3.mol.GetNumAtoms() -
            cof.bonds_made*3)
    assert (cof.mol.GetNumBonds() ==
            amine3.mol.GetNumBonds()*1 +
            aldehyde3.mol.GetNumBonds()*1 -
            cof.bonds_made*2)
    assert cof.bb_counter[amine3] == 1
    assert cof.bb_counter[aldehyde3] == 1
    assert cof.topology == stk.NoLinkerHoneycomb()
