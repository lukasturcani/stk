import os
import shutil
import networkx as nx
import subprocess as sp
import rdkit.Chem as chem

from ..classes import StructUnit, FGInfo
from .convenience_functions import MolFileError

def fg_prune(input_folder, output_folder, fg, fg_num):
    """
    Copies molecules with a given functional group between folders.
    
    Parameters
    ----------
    input_folder : str
        The full path of the folder holding V3000 .mol files. Any
        molecules with a functional group `fg` in this folder are copied
        to `output_folder`.
        
    output_folder : str
        The full path of the folder into which files with a given 
        functional group `fg` are copied to.
    
    fg : str
        The name of the functional group which the copied molecules must
        possess. The name must correspond to one of the name of a 
        functional group defined within 
        ``FGInfo.functional_groups_list``.
        
    fg_num : int
        The number of functional groups of type `fg` which the molecule
        must have in order to be copied.
        
    Returns
    -------
    None : NoneType
        
    """
    
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
    """
    Converts all .mol2 in a folder to V3000 .mol files.
    
    Parameters
    ----------
    input_folder : str
        The full path of the folder filled with .mol2 files.
    
    output_folder : str
        The full path of the folder where the V3000 .mol files should be
        placed.
        
    Returns
    -------
    None : NoneType
    
    """
    
    input_list = ['babel', '-m', '-imol2', '' , '-omol', '', '-x3']
    for x in os.listdir(input_folder):
        input_list[3] = os.path.join(input_folder, x)
        input_list[5] = os.path.join(output_folder, x).replace('.mol2', '.mol')
        sp.call(input_list)
 
def fg_distance_prune(folder, fg):
    """
    Deletes molecules with functional groups seperated by 1 atom.
    
    Parameters
    ----------
    folder : str
        The full path of the folder which holdes the molecules in a
        V3000 .mol format. The .mol files are removed from this folder.
        
    fg : str
        The name of the functional group.
    
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
    """
    Deletes molecules which contain the substructure `substruct`.
    
    Parameters
    ----------
    folder : str
        The full path of the folder from which the files are checked for
        substructure and deleted.
        
    substruct : str
        The smiles string of the substructure, which if present in a 
        molecule causes it to be deleted from `folder`.
    
    """
    
    substruct_mol = chem.MolFromSmiles(substruct)
    for file_name in os.listdir(folder):
        path = os.path.join(folder, file_name)
        mol = chem.MolFromMolFile(path)
        if mol.HasSubstructMatch(substruct_mol):
            print('Removing {}.'.format(path))
            os.remove(path)
            
    
            