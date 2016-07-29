import numpy as np
from functools import wraps
from operator import attrgetter
import itertools
import weakref
import rdkit
from rdkit import Chem as chem
from rdkit.Chem import AllChem as ac
from collections import namedtuple
from operator import attrgetter
from copy import deepcopy
import os
import math

from ..convenience_functions import dedupe, flatten

class Cached(type):   
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)        
        self.__cache = weakref.WeakValueDictionary()
    
    def __call__(self, *args):
        if args in self.__cache.keys():
            return self.__cache[args]
        else:
            obj = super().__call__(*args)
            self.__cache[args] = obj
            return obj
      
FGInfo = namedtuple('FGInfo', ['name', 'smarts', 
                               'target_atomic_num', 'heavy_atomic_num',
                               'heavy_symbol'])
class StructUnit:
    functional_group_list = [
                        
                    FGInfo("aldehyde", "C(=O)[H]", 6, 39, "Y"), 
                    FGInfo("carboxylic acid", "C(=O)O[H]", 6, 40, "Zr"),
                    FGInfo("amide", "C(=O)N([H])[H]", 6, 41, "Nb"),
                    FGInfo("thioacid", "C(=O)S[H]", 6, 42, "Mo"),
                    FGInfo("alcohol", "O[H]", 8, 43, "Tc"),
                    FGInfo("thiol", "[S][H]", 16, 44, "Ru"),
                    FGInfo("amine", "[N]([H])[H]", 7, 45, "Rh"),    
                    FGInfo("nitroso", "N=O", 7, 46, "Pd"),
                    FGInfo("boronic acid", "[B](O[H])O[H]", 5, 47, "Ag")
                             
                             ]

    def __init__(self, prist_mol_file):
        self.prist_mol_file = prist_mol_file
        self.prist_mol = chem.MolFromMolFile(prist_mol_file, 
                                             sanitize=False, 
                                             removeHs=False)
                                             
        self.prist_smiles = chem.MolToSmiles(self.prist_mol, 
                                             isomericSmiles=True,
                                             allHsExplicit=True)

        self.func_grp = next((x for x in 
                                StructUnit.functional_group_list if 
                                x.name in prist_mol_file), None)
        
        self.generate_heavy_attrs()

    def generate_heavy_attrs(self):        
        func_grp_atom_ids = flatten(self.find_functional_group_atoms())       

        self.heavy_mol = deepcopy(self.prist_mol)      
        
        for atom_id in func_grp_atom_ids:
            atom = self.heavy_mol.GetAtomWithIdx(atom_id)
            if atom.GetAtomicNum() == self.func_grp.target_atomic_num:
                atom.SetAtomicNum(self.func_grp.heavy_atomic_num)
        
        heavy_file_name = list(os.path.splitext(self.prist_mol_file))
        heavy_file_name.insert(1,'HEAVY')
        heavy_file_name.insert(2, self.func_grp.name)
        self.heavy_mol_file = '_'.join(heavy_file_name)     
        
        chem.MolToMolFile(self.heavy_mol, self.heavy_mol_file,
                          includeStereo=True, kekulize=False,
                          forceV3000=True) 

        self.heavy_smiles = chem.MolToSmiles(self.heavy_mol, 
                                             isomericSmiles=True,
                                             allHsExplicit=True)        

    def find_functional_group_atoms(self):
        """


        """
        
        func_grp_mol = chem.MolFromSmarts(self.func_grp.smarts)
        return self.prist_mol.GetSubstructMatches(func_grp_mol)        

    def shift_heavy_mol(self, x, y, z):
        conformer = chem.Conformer(self.heavy_mol.GetConformer())
        
        for atom in self.heavy_mol.GetAtoms():            
            atom_id = atom.GetIdx()
            atom_position = conformer.GetAtomPosition(atom_id)
            
            new_x = atom_position.x + x
            new_y = atom_position.y + y
            new_z = atom_position.z + z
            
            new_coords = rdkit.Geometry.rdGeometry.Point3D(new_x, 
                                                           new_y, new_z)            
            
            conformer.SetAtomPosition(atom_id, new_coords)
        
        new_heavy = deepcopy(self.heavy_mol)
        new_heavy.RemoveAllConformers()
        new_heavy.AddConformer(conformer)
        return new_heavy        
        
    def get_heavy_coords(self):
        conformer = self.heavy_mol.GetConformer()
        for atom in self.heavy_mol.GetAtoms():        
            atom_position = conformer.GetAtomPosition(atom.GetIdx())
            yield atom_position.x, atom_position.y, atom_position.z
        
class BuildingBlock(StructUnit):
    pass
        
class Linker(StructUnit):
    pass

class Cage(metaclass=Cached):
    def __init__(self, *args):
        if len(args) == 3:
            self.testing_init(*args)
        if len(args) == 4:
            self.std_init(*args)

    def std_init(self, bb_file, lk_file, topology, full_path):
        self.bb = BuildingBlock(bb_file)
        self.lk = Linker(lk_file)        
        self.topology = topology(self)
        self.full_path = full_path
        
        self.topology.build_cage()
        
    def bb_only_init(self, ):
        pass
    def lk_only_init(self, ):
        pass
    
    def same_cage(self, other):
        return (self.bb == other.bb and self.lk == other.lk and 
                                    self.topology == other.topology)
        
    def __str__(self):
        return str(self.__dict__) + "\n"
    
    def __repr__(self):
        return str(self.__dict__) + "\n"

    """
    The following methods are inteded for convenience while 
    debugging or testing and should not be used during typical 
    execution of the program.
    
    """

    def testing_init(self, bb_str, lk_str, topology_str):
        self.bb = bb_str
        self.lk = lk_str
        self.topology = topology_str

    @classmethod
    def init_empty(cls):
        obj = cls()
        string = ['a','b','c','d','e','f','g','h','i','j','k','l','m',
                  'n','o', 'p','q','r','s','t','u','v','w','x','y','z']
        obj.bb = np.random.choice(string)
        obj.lk = np.random.choice(string)
        obj.fitness = abs(np.random.sample())
        return obj



        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        