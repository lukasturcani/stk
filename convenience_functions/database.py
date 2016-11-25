import os
import shutil
import networkx as nx
import rdkit.Chem as chem
import rdkit.Chem.AllChem as ac

from ..classes import StructUnit, FGInfo
from .convenience_functions import MolFileError
from ..optimization.macromodel.macromodel_opt_funcs import structconvert

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

def categorize(path, output_dir):
    try:
        fgs = ['amine', 'aldehyde']
        dirs = ['amines2f', 'amines3f', 'amines4f',
                'aldehydes2f', 'aldehydes3f', 'aldehydes4f']
      
        mol = StructUnit(path, minimal=True)
        
        for fg in fgs:
            mol.func_grp = next((x for x in 
                                    FGInfo.functional_group_list if 
                                    x.name == fg), None)
            mol.heavy_ids = []
            mol._generate_heavy_attrs()
            fg_n = str(len(mol.heavy_ids))
            folder = next((x for x in dirs if fg_n in x and fg in x), None)
            if folder is not None:
                shutil.copy(path, os.path.join(output_dir,folder))

    except:
        print('Failed with {}.'.format(path))


def neutralize(path, output_dir, macromodel_path):
    try:
        m = chem.MolFromMol2File(path)
        if m is None:
            new_path = path.replace('.mol2', '.mol')
            structconvert(path, new_path, macromodel_path)
            m = chem.MolFromMolFile(new_path, sanitize=False)
        
        m = chem.MolToSmiles(m)            
        m = m.replace("[H]", "").replace("H", "").replace("()", "")
        m = m.replace("+", "").replace("-", "")
        mol = chem.MolFromSmarts(m)
        mol.UpdatePropertyCache()
        mol = chem.AddHs(mol)
        ac.EmbedMolecule(mol)
        name = os.path.basename(path)
        chem.MolToMolFile(mol, os.path.join(output_dir, name),
                          forceV3000=True)
        return (0, None)
    except ValueError as ex:
        return (1,path)

    except Exception as ex:
        print('{} failed with exception {}.'.format(path, ex))
        return (2,path)
        












