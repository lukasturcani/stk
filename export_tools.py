import os
import glob
import networkx as nx
import rdkit.Chem as chem
import rdkit.Chem.AllChem as ac
from functools import partial
from multiprocessing import Pool
import itertools as it

from .classes import StructUnit, FGInfo, Population, Cage
from .convenience_tools import MolFileError
from .optimization import *

def assign_keys(population):
    """
    Some old Cage .dmp files do not have `key` attribute so give it.

    Parameters
    ----------
    population : Population
        Assigns the `key` attribute to every MacroMolecule in this 
        population.

    Modifies
    --------
    population
        All the MacroMolecule instances in this Population instance
        have the `key` attribute assigned.
        
    Returns
    -------
    None : NoneType
    
    """
    
    for mem in population:
        mem.key = (str([mem.building_blocks, type(mem.topology)]) + 
                    str(mem.topology_args))
    

def redump_pop(*folders, ofolder):
    """
    Collects all the MacroMolecule .dmp files and creates a dump file.
    
    The function creates a population of all the .dmp files in `folders`
    and creates a dump file of that population. The dump file gets 
    placed in `ofolder`.
    
    Parameters
    ----------
    *folders : str
        The full paths of the folders holding the MacroMolecule .dmp
        files.
        
    ofolder : str
        The folder in which the dump file should be placed. 
    
    Modifies
    --------
    ofolder/pop_dump
        Creates a file holding the population dump at this location.    
    
    Returns
    -------
    None : NoneType
    
    """
    
    cage_paths = []
    for folder in folders:
        s = os.path.join(folder, '*.dmp')
        cage_paths.extend(glob.glob(s))
    
    pop = Population(*(Cage.load(x) for x in cage_paths))

    assign_keys(pop)      
    
    pop.dump(os.path.join(ofolder, 'pop_dump'))

def fg_prune(ifolder, fg, fg_num):
    """
    Deletes molecules without a given functional group from a folder.
    
    Parameters
    ----------
    ifolder : str
        The full path of the folder holding .mol files. Any
        molecules without functional group `fg` `fg_num` amount of times
        are deleted.

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
    
    for file_name in os.listdir(ifolder):
        if not file_name.endswith('.mol'):
            continue
        path = os.path.join(ifolder, file_name)
        mol = StructUnit(path, minimal=True)
    
        mol.func_grp = next((x for x in 
                                FGInfo.functional_group_list if 
                                x.name == fg), None)
        mol.heavy_ids = []
        mol._generate_heavy_attrs()
        if len(mol.find_functional_group_atoms()) != fg_num:
            print('Deleting {}.'.format(path))
            os.remove(path)

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
        if not file_name.endswith('.mol'):
            continue

        mol = StructUnit(path, minimal=True)

        mol.func_grp = next((x for x in 
                                FGInfo.functional_group_list if 
                                x.name == fg), None)
        
        mol.heavy_ids = []
        mol._generate_heavy_attrs()
        g = mol.graph('heavy')
        for start, end in it.combinations(mol.heavy_ids, 2):
            if nx.shortest_path_length(g, *mol.heavy_ids) < 3:
                print('Removing {}.'.format(path))
                os.remove(path)
                break

def substruct_prune(folder, substruct):
    """
    Deletes molecules which contain the substructure `substruct`.
    
    Parameters
    ----------
    folder : str
        The full path of the folder from which the files are checked for
        substructure and deleted.
        
    substruct : str
        The smarts string of the substructure, which if present in a 
        molecule causes it to be deleted from `folder`.
    
    """
    
    substruct_mol = chem.MolFromSmarts(substruct)
    for file_name in os.listdir(folder):
        if not file_name.endswith('.mol'):
            continue
        path = os.path.join(folder, file_name)
        mol = chem.MolFromMolFile(path, sanitize=False, removeHs=False)
        if mol.HasSubstructMatch(substruct_mol):
            print('Removing {}.'.format(path))
            os.remove(path)

def nonamine_groups(group_tuples, mol):
    invalid_groups = []
    for group in group_tuples:
        n = next(mol.GetAtomWithIdx(x) for x in group if 
                    mol.GetAtomWithIdx(x).GetAtomicNum() == 7)
        non_Hs_neighbors = 0
        for neighbor in n.GetNeighbors():
            if neighbor.GetAtomicNum() != 1:
                non_Hs_neighbors += 1
        
        if non_Hs_neighbors > 1:
            invalid_groups.append(group)
        
    return invalid_groups
            
def categorize(mol2block, filename, output_dir):

    try:
        fgs = [chem.MolFromSmarts('C(=O)N([H])[H]'),
               chem.MolFromSmarts('C(=O)[H]'),
               chem.MolFromSmarts('[N]([H])[H]')]
               

        dirs = ['amines2f', 'amines3f', 'amines4f',
                'aldehydes2f', 'aldehydes3f', 'aldehydes4f']

        mol = chem.MolFromMol2Block(mol2block, sanitize=False,
                                    removeHs=False)
        
        amide_matches = mol.GetSubstructMatches(fgs[0])
        if len(amide_matches) > 0:
            return 
            
        aldehyde_matches = mol.GetSubstructMatches(fgs[1])
        amine_matches = mol.GetSubstructMatches(fgs[2])
        remove = nonamine_groups(amine_matches, mol)
        amine_matches = [x for x in amine_matches if 
                         x not in remove]
        
        non0 = sum(1 for x in [amide_matches, aldehyde_matches,
                               amine_matches]
                        if len(x) > 0)
        
        if non0 > 1 or non0 == 0:
            return
        
            
        fg, fg_n = next((name, str(len(x))) for name, x in  [('amines', amine_matches), 
                                        ('aldehydes', aldehyde_matches)]
                     if len(x) > 0)
        

        
        folder = next((x for x in dirs if fg_n in x and fg in x), None)
        if folder is not None:
            oname = os.path.join(output_dir,folder,filename)
            with open(oname, 'w') as f:
                f.write(mol2block)

    except:

        with open('/home/lukas/database/fails{}.mol2'.format(filename[0]), 'a') as f:
            f.write(mol2block)

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
                fn = "{}_{}.mol2".format(n1,n2)
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
        timeout=False, temp=300, confs=500, eq_time=50, sim_time=500)    
    
    with Pool() as p:
        p.map(md_opt, macro_mols)
