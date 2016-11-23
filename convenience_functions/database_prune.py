import os
import shutil
import networkx as nx
import subprocess as sp
import rdkit.Chem as chem

from ..classes import StructUnit, FGInfo
from .convenience_functions import MolFileError

def fg_prune(input_folder, output_folder, fg, fg_num):
    for file_name in os.listdir(input_folder):
        path = os.path.join(input_folder, file_name)
        try:
            mol = StructUnit(path, minimal=True)
        
        except MolFileError as error:
            print('V3000 {}.'.format(path))
            continue
    
        mol.func_grp = next((x for x in 
                                FGInfo.functional_group_list if 
                                x.name == fg), None)
        mol.heavy_ids = []
        mol._generate_heavy_attrs()
        if len(mol.find_functional_group_atoms()) == fg_num:
            print('Moving {}.'.format(path))
            shutil.copy(path, output_folder)
            
def mol2_to_mol(input_folder, output_folder):
    input_list = ['babel', '-m', '-imol2', '' , '-omol', '', '-x3']
    for x in os.listdir(input_folder):
        input_list[3] = os.path.join(input_folder, x)
        input_list[5] = os.path.join(output_folder, x).replace('.mol2', '.mol')
        sp.call(input_list)
 
def fg_distance_prune(folder, fg):
    """
    Deletes molecules with fg seperated by 1 atom.
    
    """
    
    for file_name in os.listdir(folder):
        path = os.path.join(folder, file_name)
        try:
            mol = StructUnit(path, minimal=True)
            
        except MolFileError as error:
            print('V3000 {}.'.format(path))
            continue
    
        mol.func_grp = next((x for x in 
                                FGInfo.functional_group_list if 
                                x.name == fg), None)
        
        mol.heavy_ids = []
        mol._generate_heavy_attrs()
        
        g = mol.graph('heavy')
        if nx.shortest_path_length(g, *mol.heavy_ids) < 3:
            print('Removing {}.'.format(path))
            os.remove(path)

def substurct_prune(folder, substruct):
    substruct_mol = chem.MolFromSmiles(substruct)
    for file_name in os.listdir(folder):
        path = os.path.join(folder, file_name)
        mol = chem.MolFromMolFile(path)
        if mol.HasSubstructMatch(substruct_mol):
            print('Removing {}.'.format(path))
            os.remove(path)
            
    
            