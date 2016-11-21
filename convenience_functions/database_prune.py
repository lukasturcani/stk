import os

from ..classes import StructUnit, FGInfo
from .convenience_functions import MolFileError

def fg_prune(folder, fg, fg_num):
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
        if len(mol.find_functional_group_atoms()) != fg_num:
            print('Removing {}.'.format(path))
            os.remove(path)

