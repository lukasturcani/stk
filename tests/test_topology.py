import networkx as nx
import rdkit.Chem as chem
import rdkit

from MMEA.classes import BuildingBlock, Linker, MacroMolecule, FourPlusSix, FGInfo
from .test_struct_unit import get_mol_file

bb_file = next(x for x in get_mol_file() 
                                    if 'amine3f_14.mol' in x)
lk_file = next(x for x in get_mol_file() 
                                    if 'aldehyde2f_3.mol' in x) 

double_bond_file = next(x for x in get_mol_file() 
                                    if 'double_bond_test.mol' in x) 

single_bond_file = next(x for x in get_mol_file() 
                                    if 'single_bond_test.mol' in x)
                                    
distance_test_file = next(x for x in get_mol_file() 
                                    if 'min_distance_test.mol' in x)
                                        
partner_pool1_file = next(x for x in get_mol_file() 
                                    if 'unpaired_diff_atoms_test1.mol' in x)


                                        
bb = BuildingBlock(bb_file)
lk = Linker(lk_file)   

def test_extract_heavy_atoms():

    mol = MacroMolecule.__new__(MacroMolecule)
    mol.building_blocks = (bb, lk)
    mol.topology = FourPlusSix(mol)
    mol.topology.place_mols()    
    
    heavy_graph = mol.get_heavy_as_graph()
    molecules = [x.nodes() for x in 
                        nx.connected_component_subgraphs(heavy_graph)]
    assert not all(len(mol) == 2 or len(mol) == 3 for mol in molecules)
        
    heavy_mols = mol.topology.extract_heavy_atoms(molecules)
    
    assert all(len(mol) == 2 or len(mol) == 3 for mol in heavy_mols)
    
def test_final_sub():
 
    
    mol = MacroMolecule.__new__(MacroMolecule)
    mol.building_blocks = (bb, lk)
    mol.topology = FourPlusSix(mol)
    mol.topology.place_mols()    
    mol.topology.join_mols()
    assert  not hasattr(mol, 'prist_mol')
    mol.topology.final_sub()
    assert hasattr(mol, 'prist_mol')
    
    for atom in mol.prist_mol.GetAtoms():
        assert atom.GetAtomicNum() not in FGInfo.heavy_atomic_nums

def test_determine_bond_type():
    mol = MacroMolecule.__new__(MacroMolecule)
    mol.topology = FourPlusSix(mol)
    mol.heavy_mol = chem.MolFromMolFile(double_bond_file)

    assert mol.topology.determine_bond_type(0,1) == rdkit.Chem.rdchem.BondType.DOUBLE    
    mol.heavy_mol = chem.MolFromMolFile(single_bond_file)
    assert mol.topology.determine_bond_type(0,1) == rdkit.Chem.rdchem.BondType.SINGLE

    
def test_min_distance_partner():
    mol = MacroMolecule.__new__(MacroMolecule)
    mol.topology = FourPlusSix(mol)
    mol.heavy_mol = chem.MolFromMolFile(distance_test_file)
    
    min_partner = mol.topology.min_distance_partner(0, [1,2,3,4])
    assert mol.heavy_mol.GetAtomWithIdx(min_partner).GetAtomicNum() == 7    
    
def test_unpaired_diff_element_atoms():
    mol = MacroMolecule.__new__(MacroMolecule)
    mol.topology = FourPlusSix(mol)
    mol.topology.paired = set()
    mol.topology.paired_mols = set()
    mol.heavy_mol = chem.MolFromMolFile(partner_pool1_file)    
    assert len(list(mol.topology.unpaired_diff_element_atoms(0,[[0],[1],[2],[3],[4],[5],[6],
                                                                [7],[8],[9],[10],[11],[12],[13],[14]]))) == 14
    assert len(list(mol.topology.unpaired_diff_element_atoms(6,[[0],[1],[2],[3],[4],[5],[6],
                                                                [7],[8],[9],[10],[11],[12],[13],[14]]))) == 6

    
    
    
    
    
    