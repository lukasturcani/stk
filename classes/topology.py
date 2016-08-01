from rdkit import Chem as chem
import math
import numpy as np
import itertools as itertools
import re

from .molecular import FGInfo

class MolFileData:
    """
    Used to store data collected from a ``.mol`` file during assembly.
    
    This class is used during assembly to keep track of the content
    which should be written to the final ``.mol`` file of the assembled
    molecule.
    
    Attributes
    ----------
    mol_file_content : str
        This attribute holds a string with all of the ``.mol`` file 
        content. This content is edited to include the additional bonds
        need to assemble a larger molecule.
    
    count_line : str
        Represents the ``count line`` in the ``.mol`` file. This line
        holds the atom and bond number. The bond number needs to be
        edited during course of assembly.
        
    at_num : str
        This attribute holds the number of atoms defined in the ``.mol``
        file. It is held a string for convenience as it is read from a
        string and later written to one.

    bond_number : int
        Holds the number of bonds in the ``.mol`` file.

    atom1_list : list of ``Atom`` instances
        All the heavy atoms representing one of the functional groups
        are placed into ``Atom`` instances which are stored here. Each
        functional group has its own such attribute/list.
        
    atom2_list : list of ``Atom`` instances      
        All the heavy atoms representing one of the functional groups
        are placed into ``Atom`` instances which are stored here. Each
        functional group has its own such attribute/list.
    
    """
    
    __slots__ = ['mol_file_content', 'count_line', 'at_num', 
                 'bond_number', 'atom1_list', 'atom2_list']
    
    def __init__(self, mol_file_content, count_line, at_num, 
                 bond_number, atom1_list, atom2_list):
         self.mol_file_content = mol_file_content
         self.count_line = count_line
         self.at_num = at_num
         self.bond_number = bond_number
         self.atom1_list = atom1_list
         self.atom2_list = atom2_list
    
class Atom(object):
    """
    Used to store and manipulate data collected from a ``.mol`` file.    

    This class is used during assembly to keep track of which atoms need
    to be linked together by writing bonds to the ``.mol`` file. It is 
    also used to find the distances between heavy atoms in the ``.mol``
    file. This allows the creation of bonds between the correct atoms. A
    full list of operations performed by this class can be found by 
    examining its entire definition.
    
    Class attributes
    ----------------
    bb :
    
    bb_heavy_atoms_per_molecule :     
    
    lk :
    
    lk_heavy_atoms_per_molecule :
    
    linked_mols :    
    
    
    Attributes
    ----------
    element : 
    
    number :
    
    heavy_atom_num :
    
    x :
    
    y : 
    
    z :
    
    min_parter :
    
    min_distance :
    
    distances : set
    
    paired : bool

    """
    bb = None
    bb_heavy_atoms_per_molecule = None
    lk = None    
    lk_heavy_atoms_per_molecule = None
    linked_mols = []    
    
    def __init__(self, element, number, heavy_atom_num, x, y, z):
        self.element = element
        self.number = number
        self.heavy_atom_num = heavy_atom_num
        self.x = x
        self.y = y
        self.z = z
        self.min_partner = 10**9
        self.min_distance = 10**9
        self.distances = {}        
        self.paired = False

    @staticmethod
    def extract_mol_file_data(mol_file):
    
        mol_file_content = ""
        with open(mol_file, "r") as mol_file:
            
            atom1_element = None    
            atom2_element = None
            atom1_list = []
            atom2_list = []
        
            take_atom = False
            take_bond = False
            write_line = True
        
            type1_heavy_atom_number = 1
            type2_heavy_atom_number = 1
        
            for raw_line in mol_file:           
                line = raw_line.split()
                
                if "M  V30 END BOND" in raw_line:
                    write_line = False            
                    
                if "M  V30 COUNTS" in raw_line:
                    count_line = raw_line
                    at_num = line[3]            
                    
                if write_line == True:            
                    mol_file_content += raw_line
                        
                if "M  V30 BEGIN ATOM" in raw_line:
                    take_atom = True
                    continue
                
                if "M  V30 END ATOM" in raw_line:
                    take_atom = False
                    continue
                
                if "M  V30 BEGIN BOND" in raw_line:
                    take_bond = True
                    continue
                
                if take_atom == True:                 

                    atom_id = line[2]
                    atomic_symbol = line[3]
                    atom_x = float(line[4])
                    atom_y = float(line[5])
                    atom_z = float(line[6])

                    atom_is_heavy = atomic_symbol in Topology.heavy_symbols
                    atom1_found = atom1_element != None
                    atom2_found = atom2_element != None
                    
                    if atom_is_heavy and not atom1_found and not atom2_found:
                            
                            atom1_element = atomic_symbol
                            atom1_list.append(Atom(atomic_symbol, atom_id, 
                                              type1_heavy_atom_number, 
                                              atom_x, atom_y, atom_z))
                            
                            type1_heavy_atom_number += 1                    
                            continue
                        
                    if (atom_is_heavy and atom1_found and not atom2_found and 
                                    atomic_symbol != atom1_element):
                            
                            atom2_element = atomic_symbol
                            atom2_list.append(Atom(atomic_symbol, atom_id, 
                                              type2_heavy_atom_number, 
                                              atom_x, atom_y, atom_z))
                            type2_heavy_atom_number += 1                        
                            continue
                        
                    if (atom_is_heavy and atom1_found and 
                                    atom1_element == atomic_symbol):
                            
                            atom1_list.append(Atom(atomic_symbol, atom_id, 
                                              type1_heavy_atom_number, 
                                              atom_x, atom_y, atom_z))
                            type1_heavy_atom_number += 1
                            continue
                            
                    if (atom_is_heavy and atom2_found and 
                                    atom2_element == atomic_symbol):
                             
                            atom2_list.append(Atom(atomic_symbol, atom_id, 
                                                   type2_heavy_atom_number, 
                                              atom_x, atom_y, atom_z))
                            type2_heavy_atom_number += 1
                            continue  

                if take_bond == True and len(line) == 6:                   
                    bond_number = int(line[2])
        
        return MolFileData(mol_file_content, count_line, at_num, 
                           bond_number, atom1_list, atom2_list)
        
    def assign_molecule_number(self):       
        if self.element == Atom.bb:
            self.mol_number = math.ceil(self.heavy_atom_num / 
                                    Atom.bb_heavy_atoms_per_molecule)  
        if self.element == Atom.lk:
            self.mol_number = math.ceil(self.heavy_atom_num / 
                                    Atom.lk_heavy_atoms_per_molecule)
        
        return 1
    
    def distance(self, atom2):
        x_diff_sq = (self.x - atom2.x) ** 2
        y_diff_sq = (self.y - atom2.y) ** 2
        z_diff_sq = (self.z - atom2.z) ** 2
        r = np.sqrt(x_diff_sq + y_diff_sq + z_diff_sq)
        self.distances[atom2] = r
        atom2.distances[self] = r
        
        return r
        
    def pair_up_nitro(self):
        """
        This case is specific for the Nitroso cages, for which a double 
        counting for the atoms needs to be avoided.
        """
        if self.paired == True:
            return None
        
        if self.paired == False: 
            while True:
                min_partner = min(self.distances, 
                                  key=self.distances.get)
                
                if (min_partner.paired == False and 
                (self.mol_number, min_partner.mol_number) not in 
                Atom.linked_mols and 
                self.mol_number != min_partner.mol_number and 
                ((self.number, self.element, min_partner.number, 
                                              min_partner.element) 
                or (min_partner.number, min_partner.element, 
                    self.number, self.element)) not in 
                    Atom.linked_atoms):

                    self.min_partner = min_partner
                    self.paired = True
                    min_partner.min_partner = self
                    min_partner.paired = True
                    Atom.linked_mols.append((self.mol_number, 
                                            min_partner.mol_number))
                    Atom.linked_atoms.append((self.number, 
                                              self.element, 
                                              min_partner.number, 
                                              min_partner.element))
                    
                    return None
                else:
                    del self.distances[min_partner]

    def pair_up_v3_v2(self):
        if self.paired == True:
            return None
        
        if self.paired == False:
            
#                    print("  DISTANCES: \n  {} {}".format(self.number, self.distances.values()))
            min_partner = min(self.distances, 
                              key=self.distances.get)
#                    print("MIN PARTNER PAIRED? {}".format(min_partner.paired))
            
            if (min_partner.paired == False and 
            (self.mol_number, min_partner.mol_number) not in 
            Atom.linked_mols):
#                    if min_partner.paired == False:
                self.min_partner = min_partner
                self.paired = True
                min_partner.min_partner = self
                min_partner.paired = True
                Atom.linked_mols.append((self.mol_number, 
                                        min_partner.mol_number))
                Atom.linked_atoms.append((self.number, 
                                          self.element, 
                                          min_partner.number, 
                                          min_partner.element))
                
                return None
            else:
                del self.distances[min_partner]
                    
    def pair_up_v4_v2(self):
        if self.paired == True:
            return None
        
        if self.paired == False:
            while True:

                min_partner = min(self.distances, 
                                  key=self.distances.get)
                
                if (min_partner.paired == False and 
                (self.mol_number, min_partner.mol_number) not in 
                Atom.linked_mols):

                    self.min_partner = min_partner
                    self.paired = True
                    min_partner.min_partner = self
                    min_partner.paired = True
                    Atom.linked_mols.append((self.mol_number, 
                                        min_partner.mol_number))
                    Atom.linked_atoms.append((self.number, 
                                              self.element, 
                                             min_partner.number, 
                                           min_partner.element))
                    
                    return None
                else:
                    del self.distances[min_partner]        

    def pair_up_1plus1(self):
        if self.paired == True:
            return None
            
        if self.paired == False:
            while True:

                min_partner = min(self.distances, 
                                  key=self.distances.get)
                
                if (min_partner.paired == False and 
                ((self.number, self.element, min_partner.number, 
                  min_partner.element) 
                or (min_partner.number, min_partner.element, 
                    self.number, self.element)) not in 
                    Atom.linked_atoms):

                    self.min_partner = min_partner
                    self.paired = True
                    min_partner.min_partner = self
                    min_partner.paired = True
                    Atom.linked_mols.append((self.mol_number, 
                                        min_partner.mol_number))
                    Atom.linked_atoms.append((self.number, 
                                              self.element,
                                             min_partner.number, 
                                           min_partner.element))
                    
                    return None
                else:
                    del self.distances[min_partner]

    def pair_up_2plus2(self):
        """
        For 2+2 we have been experiencing ValueError for specific bb and lks.
        """
        
        if self.paired == True:
            return None

        if self.paired == False:
            min_partner = min(self.distances, 
                              key=self.distances.get)                    
            
            bonds_between_mols = (Atom.linked_mols.count(
            (self.mol_number, min_partner.mol_number)) + 
            Atom.linked_mols.count((min_partner.mol_number, 
                                            self.mol_number)))

            #In this topology each molecule cannot bind more than twice to an other molecule
            if (min_partner.paired == False and
                    bonds_between_mols <= 2):
                self.min_partner = min_partner
                self.paired = True
                min_partner.min_partner = self
                min_partner.paired = True
                Atom.linked_mols.append((self.mol_number, 
                                         min_partner.mol_number))
                Atom.linked_atoms.append((self.number, 
                                          self.element, 
                                          min_partner.number, 
                                          min_partner.element))
                
                return None
            else:
                del self.distances[min_partner]

    def pair_up_4plus4(self):
        if self.paired == True:
            return None
        
        if self.paired == False: 
            while True:
                min_partner = min(self.distances, 
                                  key=self.distances.get)
                
                if (min_partner.paired == False and 
                (self.mol_number, min_partner.mol_number) not in 
                Atom.linked_mols):

                    self.min_partner = min_partner
                    self.paired = True
                    min_partner.min_partner = self
                    min_partner.paired = True
                    Atom.linked_mols.append((self.mol_number, 
                                        min_partner.mol_number))
                    Atom.linked_atoms.append((self.number, 
                                              self.element, 
                                             min_partner.number, 
                                        min_partner.element))
                    
                    return None
                else:
                    del self.distances[min_partner]

class Topology:
    heavy_symbols = {x.heavy_symbol for x 
                        in FGInfo.functional_group_list}

    
    def __init__(self, cage):
        self.cage = cage
        

    def build_cage(self):
        self.place_mols()
        self.join_mols()
        self.final_sub()
        
    def build_heavy_cage(self):
        self.join_mols()
        
    def bulid_prist_cage(self):
        pass

    def join_mols(self):
        """
        takens an input file with disconnected moleucles and connects them
        """           
        
        
        mol_file_data = Atom.extract_mol_file_data(
                                               self.cage.heavy_mol_file)
            
        Atom.bb_heavy_atoms_per_molecule = self.heavy_atoms_per_bb  
        Atom.lk_heavy_atoms_per_molecule = self.heavy_atoms_per_lk   
        
        Atom.bb = self.cage.bb.func_grp.heavy_symbol
        Atom.lk = self.cage.lk.func_grp.heavy_symbol
        
        Atom.linked_mols = []
        Atom.linked_atoms = []
        
        for atom in itertools.chain(mol_file_data.atom1_list,
                                    mol_file_data.atom2_list):
            atom.assign_molecule_number()
        
        for atom1 in mol_file_data.atom1_list:
            for atom2 in mol_file_data.atom2_list:
                atom1.distance(atom2)   
    
        for atom1 in mol_file_data.atom1_list:
            for atom2 in mol_file_data.atom2_list:
                self.pair_up_func(atom1)
        
        for atom in Atom.linked_atoms:
            mol_file_data.bond_number += 1
            double_bond_present = [atom[1] in tup and atom[3] in tup for 
                                    tup in FGInfo.double_bond_combs]
            
            if True in double_bond_present:
                bond_order = "2"
            else:
                bond_order = "1"
                
            mol_file_data.mol_file_content += "M  V30 {2} {3} {0} {1}\n".format(
                              atom[0], atom[2], mol_file_data.bond_number, bond_order)
        
        
        mol_file_data.mol_file_content += ("M  V30 END BOND\nM  V30 END CTAB\nM"  
                                                                 "  END")

        p = re.compile(r" VAL=.")
        mol_file_data.mol_file_content = re.sub(p, "", mol_file_data.mol_file_content)

        mol_file_data.mol_file_content = mol_file_data.mol_file_content.replace(mol_file_data.count_line, 
             "M  V30 COUNTS {0} {1} 0 0 0\n".format(mol_file_data.at_num,mol_file_data.bond_number))
        
        new_mol_file_name = self.cage.heavy_mol_file
        new_mol_file = open(new_mol_file_name, "w")
        new_mol_file.write(mol_file_data.mol_file_content)
        new_mol_file.close()
        
        
    def final_sub(self):
        with open(self.cage.heavy_mol_file, "r") as add:
            new_file= ""
            for line in add:
                line = line.replace(self.cage.bb.func_grp.heavy_symbol,
                                    self.cage.bb.func_grp.target_symbol)
                line = line.replace(self.cage.lk.func_grp.heavy_symbol,
                                    self.cage.lk.func_grp.target_symbol)
                new_file += line
    
        with open(self.cage.prist_mol_file, "w") as f:
            f.write(new_file)
        

      
class FourPlusSix(Topology):
    def __init__(self, cage):
        super().__init__(cage)
        self.heavy_atoms_per_bb = 3
        self.heavy_atoms_per_lk = 2
        
        self.bb_num = 4
        self.lk_num = 6
        
        self.pair_up_func = Atom.pair_up_v4_v2
        
    def place_mols(self):
        bb_placement = self.place_bbs()
        lk_placement = self.place_lks()
        self.cage.heavy_mol = chem.CombineMols(bb_placement, 
                                                   lk_placement) 
                                                   
        chem.MolToMolFile(self.cage.heavy_mol, self.cage.heavy_mol_file,
                          includeStereo=True, kekulize=False,
                          forceV3000=True)
        
    
    def place_bbs(self):
        position1 = self.cage.bb.shift_heavy_mol(0,50,0)
        position2 = self.cage.bb.shift_heavy_mol(0,0,25)
        position3 = self.cage.bb.shift_heavy_mol(25,0,-25)
        position4 = self.cage.bb.shift_heavy_mol(-25,0,-25)

        combined_mol = chem.CombineMols(position1, position2)
        combined_mol2 = chem.CombineMols(position3, position4)
        return chem.CombineMols(combined_mol, combined_mol2)       
        
    def place_lks(self):
        position1 = self.cage.lk.shift_heavy_mol(0, 0, -25)
        position2 = self.cage.lk.shift_heavy_mol(-12.5, 0, 0)
        position3 = self.cage.lk.shift_heavy_mol(12.5, 0, 0)
        
        position4 = self.cage.lk.shift_heavy_mol(-12.5, 25, -12.5)
        position5 = self.cage.lk.shift_heavy_mol(12.5, 25, -12.5)
        position6 = self.cage.lk.shift_heavy_mol(0, 25, 25)

        combined_mol = chem.CombineMols(position1, position2)
        combined_mol2 = chem.CombineMols(position3, position4)        
        combined_mol3 = chem.CombineMols(position5, position6)
        
        combined_mol4 = chem.CombineMols(combined_mol, combined_mol2)
        return chem.CombineMols(combined_mol3, combined_mol4)