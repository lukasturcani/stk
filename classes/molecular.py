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


_FGInfo = namedtuple('FGInfo', ['name', 'smarts', 
                                'target_atomic_num', 'heavy_atomic_num',
                                'target_symbol', 'heavy_symbol'])                                
class FGInfo(_FGInfo):
    """
    Contains key information for substitutions.
    
    
    
    """
    pass

class StructUnit:
    """
    Represents the building blocks of molecules examined by MMEA.
    
    ``Building blocks`` in this case refers to the smallest molecular 
    unit of the assembled molecules (such as cages) examined by MMEA. 
    This is not the be confused with building-blocks* of cages. 
    Building-blocks* of cages are examples of the ``building blocks`` 
    referred to here. To be clear, the ``StructUnit`` class represents 
    all building blocks of the molecules, such as linkers and 
    building-blocks* of cages.
    
    To avoid confusion, in the documentation general building blocks 
    represented by ``StructUnit`` are referred to as `building blocks`` 
    while building-blocks* of cages are always referred to as 
    ``building-blocks*``. 
    
    The goal of this class is the conveniently store and perform 
    operations on the building blocks of assembled molecules. The 
    class stores information regarding the rdkit instance of the 
    building block, the ``SMILES`` string and the location of its 
    ``.mol`` file. See the attributes section of this docstring for 
    more details.
    
    This class also takes care of perfoming substitutions of the 
    functional groups in the building blocks via the 
    `_generate_functional_group_atoms` method. This method is 
    automatically invoked by the initialzer, so each initialized
    instance of ``StructUnit`` should atomatically have all of the 
    attributes associated with the substituted version of the molecular 
    building block.
    
    More information regarding what operations the class supports can be
    found by examining the methods documented below. A noteworthy 
    example is the `shift_heavy_mol` method. This method is invoked
    by other processes in MMEA (such as in the creation of assembled 
    molecules - see the `place_mols` documentation of the ``Topolgy`` 
    class) and is generally very useful. Similar class such as 
    `set_heavy_mol_position` may be added in the future. Refer to the 
    documentation of `shift_heavy_mol` below for more details. Note that 
    this paragraph is not an exhaustive list of useful operations.
    
    The class is intended to be inherited from. As mentioned 
    before, ``StructUnit`` is a general building block. If one wants to
    represent a specific building block, such as a linker or 
    building-block* (of a cage) a new class should be created. This new
    class will will inherit ``StructUnit``. In this way, any operations 
    which apply generally to building blocks can be stored here and any
    which apply specifically to one kind of building block such as a 
    linker or building-block* can be placed within its own class.
    
    Consider a useful example of this approach. When setting the 
    coordinates of linkers or building-blocks* during assembly of a 
    cage, it is necessary to know if the molecule you are placing is a 
    building-block* or linker. This is because a building-block* will go 
    on vertex (in this example, this may or may not be generally true)
    and a linker will go on an edge. 
    
    Assume that there is a ``Linker`` and a ``BuildingBlock`` class 
    which inherit from ``StructUnit``. As luck would have it, these 
    classes are in fact implemented in MMEA. Even if nothing is present
    in the class definition itself, both classes will have all the 
    attributes and methods associated with ``StructUnit``. This means
    the positions of the rdkit molecules held in instances of those 
    classes can be shifted with the `shift_heavy_mol` method.
    
    By running:
    
        >>> isinstance(your_struct_unit_instance, Linker)
        
    you can determine if the molecule you are dealing with is an example
    of a building-block* or linker of a cage. As a result you can easily
    choose to run the correct function which shifts the coordinates 
    either to a vertex or an edge.
    
    A final note on the intended use. Each instance of an assembled 
    molecule class (such as an instance of the ``Cage`` class) will have
    one instance of each class derived from ``StructUnit`` at most. It 
    holds information which applies to every building-block* or linker
    present in a class. As a result it does not hold information 
    regarding how a individual building-blocks* and linkers are joined
    up in a cage. That is the cage's problem. Specifically cage's 
    `topology` attributes problem.
    
    In summary, the intended use of this class is to answer questions
    such as (not exhaustive):
        
        > What basic structural units were used in the assembly of this 
          cage?
        > Which functional group was substituted in building-blocks*
          of this cage? 
        > Which atom was substituted for which in the linker? (Note that
          this question is delegated to the ``FGInfo`` instance held in 
          the `func_grp` attribute of a ``StructUnit`` instance)
        > Where is the ``.mol`` file represnting a single 
          building-block* of the cage located?
        > Where is the ``.mol`` file represnting the a single 
          building-block* of the cage, after it has been substituted 
          with a heavy atom, located?
        > Give me an rdkit instance of the molecule which represents the
          building-block* of a cage. Before and after 
          it has been substituted.
        > Give me an rdkit instance of the molecule which represents a
          a single linker of a cage, at postion ``(x,y,z)``.
          
    Questions which this class should not answer include:
    
        > How many building-blocks* does this cage have?
        > What is the position of a linker within this cage?
        > Create a bond between this ``Linker`` and ``BuildingBlock``.

    Class attributes
    ----------------
    functional_groups_list : list of FGInfo instances
        This list holds all ``FGInfo`` instances used by MMEA. If a new
        functional group is to be used by MMEA, a new ``FGInfo`` 
        instance must be added to this list. For details on extending 
        the MMEA to use more functional groups consult the ``FGInfo`` 
        class string or the developer's guide.

    Attributes
    ----------
    prist_mol_file : str
        The full path of the ``.mol`` file (V3000) holding the 
        unsubstituted molecule. This is the only attribute which needs 
        to be provided to the initializer. The remaining attributes have 
        values derived from this ``.mol`` file.
        
    prist_mol : rdkit.Chem.rdchem.Mol
        This is an ``rdkit molecule type``. It is the rdkit instance
        of the molecule held in `prist_mol_file`.
        
    prist_smiles : str
        This string holds the ``SMILES`` code of the unsubstituted form
        of the molecule.
        
    heavy_mol_file : str
        The full path of the ``.mol`` file (V3000) holding the 
        substituted molecule. This attribute is initialized by the 
        initializer indirectly when it calls the `generate_heavy_attrs` 
        method. 
    
    heavy_mol : rdkit.Chem.rdchem.Mol
        The rdkit instance of the substituted molecule. Generated by 
        the initializer when it calls the `generate_heavy_attrs` method.
        
    heavy_smiles : str
        A string holding the ``SMILES`` code of the substituted version
        of the molecule.
    
    func_grp : FGInfo
        This attribute holds an instance of ``FGInfo``. The ``FGInfo``
        instance holds the information regarding which functional group
        was substituted in the pristine molecule and which atom was 
        substituted for which. Furthermore, it also holds the atomic 
        numbers of the atom which was substitued and the one used in its 
        palce. For details on how this information is stored see the 
        ``FGInfo`` class string.
    
    """
    
    functional_groups_list = [
                        
            FGInfo("aldehyde", "C(=O)[H]", 6, 39, "C", "Y"), 
            FGInfo("carboxylic acid", "C(=O)O[H]", 6, 40, "C", "Zr"),
            FGInfo("amide", "C(=O)N([H])[H]", 6, 41, "C", "Nb"),
            FGInfo("thioacid", "C(=O)S[H]", 6, 42, "C", "Mo"),
            FGInfo("alcohol", "O[H]", 8, 43, "O", "Tc"),
            FGInfo("thiol", "[S][H]", 16, 44, "S", "Ru"),
            FGInfo("amine", "[N]([H])[H]", 7, 45, "N", "Rh"),    
            FGInfo("nitroso", "N=O", 7, 46, "N", "Pd"),
            FGInfo("boronic acid", "[B](O[H])O[H]", 5, 47, "B", "Ag")
                             
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
        
        self._generate_heavy_attrs()

    def _generate_heavy_attrs(self):        
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
    """
    Holds information about the building-blocks* of a cage.
    
    """
    
    pass
        
class Linker(StructUnit):
    """
    Holds information about the likners of a cage.
    
    """
    
    pass

class Cage(metaclass=Cached):
    def __init__(self, *args):
        if len(args) == 3:
            self.testing_init(*args)
        if len(args) == 4:
            self.std_init(*args)

    def std_init(self, bb_file, lk_file, topology, prist_mol_file):
        self.bb = BuildingBlock(bb_file)
        self.lk = Linker(lk_file)        
        self.topology = topology(self)
        self.prist_mol_file = prist_mol_file
        
        heavy_mol_file = list(os.path.splitext(prist_mol_file))
        heavy_mol_file.insert(1,'HEAVY')        
        self.heavy_mol_file = '_'.join(heavy_mol_file) 
        
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



        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        