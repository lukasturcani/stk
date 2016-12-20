import os
import shutil
import networkx as nx
import rdkit.Chem as chem
import rdkit.Chem.AllChem as ac
from functools import partial
from multiprocessing import Pool

from .classes import StructUnit, FGInfo
from .convenience_tools import MolFileError
from .optimization import *


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
            print('Copying {}.'.format(path))
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
import gc
def categorize(path, filename, output_dir):
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
                oname = os.path.join(output_dir,folder,filename)
                with open(oname, 'w') as f:
                    f.write(path)


    except Exception as ex:
        with open(r'C:\Users\lukas\Dropbox\fails{}.mol2'.format(filename[0]), 'a') as f:
            f.write(path)

def mol_file_iter(mol_file):
    mol_block = ''
    for line in mol_file:
        if 'ROOT' not in line:
            mol_block += line
        else:
            mol_block += line
            yield mol_block
            mol_block = ''
        
def categorize_folder(ifolder, ofolder):
    for n1, filename in enumerate(os.listdir(ifolder)):
        print(n1)
        with open(os.path.join(ifolder,filename), 'r') as f:
            for n2, mol_block in enumerate(mol_file_iter(f)):
                fn = "{}{}.mol2".format(n1,n2)
                categorize(mol_block, fn, ofolder)



        
def neutralize(path, output_dir, macromodel_path):
    try:
        m = chem.MolFromMol2File(path)
        if m is None:
            new_path = path.replace('.mol2', '.mol')
            structconvert(path, new_path, macromodel_path)
            m = chem.MolFromMolFile(new_path, sanitize=False)
        
        m = chem.RemoveHs(m)         
        name = os.path.basename(path).replace(".mol2", ".mol")
        b = chem.MolToMolBlock(m, forceV3000=True, kekulize=False)
        b = b.replace("CHG=-1", "").replace("CHG=1", "")
        with open(os.path.join(output_dir, name), 'w') as mf:
            mf.write(b)
        
        return (0, None)
    except ValueError as ex:
        return (1,path)

    except Exception as ex:
        print('{} failed with exception {}.'.format(path, ex))
        return (2,path)

def optimize_folder(path, macromodel_path):
    
    # First make a list holding all macromolecule objects to be 
    # optimzied. Because the objects need to be initialized from a .mol
    # file a StructUnit instance not a MacroMolecule instance is used.
    # minimal = True, because it is all that is needed to run 
    # optimziations, plus you want to avoid doing functional group
    # substitutions.
    names = [os.path.join(path, file_name) for file_name in 
             os.listdir(path) if file_name.endswith(".mol")]    
    macro_mols = [StructUnit(file_path, minimal=True) for file_path in 
                                                                  names]

    # .mol files often get moved around. Make sure the `prist_mol_file`
    # attribute is updated to account for this.
    for name, macro_mol in zip(names, macro_mols):
        macro_mol.prist_mol_file = name

    md_opt = partial(macromodel_md_opt, macromodel_path=macromodel_path, 
        timeout=False, temp=1000, confs=1000, eq_time=100, sim_time=10000)    
    
    with Pool() as p:
        p.map(md_opt, macro_mols)
