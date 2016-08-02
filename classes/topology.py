from rdkit import Chem as chem
import math
import numpy as np
import itertools as itertools
import re

from MMEA.classes.molecular import FGInfo, BuildingBlock, Linker


class Topology:
    """
    Represents the topology of an assembled molecule.
    
    The ``Topology`` class is concerned with how individual building 
    blocks are placed and connected in space to form an assembled 
    molecule used by MMEA. It also takes care of assmebling these
    molecules, though some tasks are delegated to the ``Atom`` class for
    this.
    
    This class directly defines any operations and attributes that are
    needed by any topology, be it a tetrahedron, octahedron or even a
    polymer. However, this class is not used directly by MMEA. It is
    intended to be inherited from. Any individual within MMEA will have
    a `topology` attribute which refers to an instnace of a class 
    derived from this. Derived classes of ``Topology`` define things
    specific to that one topology. For example the number of functional
    groups in building block molecules used to make that specific
    topology. Additionally, each derived class of ``Topology`` must 
    define which ``pair_up`` function defined in the ``Atom`` class it
    uses for pairing up its heavy atoms.
    
    Finally, each class derived from ``Topology`` must define methods
    which place building blocks in the correct positions in the ``.mol``
    file.
    
    Class attributes
    ----------------
    heavy_symbols : set of str
        The atomic symbols of heavy elements used in substitutions by
        MMEA.
    
    Attributes
    ----------
    macro_mol : MacroMolecule
        The macromolecule instance which has the given topology. This 
        provides easy access to the cages attributes to the ``Topology`` 
        instance.
    
    """
    heavy_symbols = {x.heavy_symbol for x 
                        in FGInfo.functional_group_list}
    
    def __init__(self, macro_mol):
        self.macro_mol = macro_mol
        

    def build(self):
        """
        Creates ``.mol`` files of the heavy and pristine macromolecules.
        
        Modifies
        --------
        .mol files
            This function modifies 2 ``.mol`` files. It creates the 
            ``.mol`` files holding the assembled heavy and pristine cage
            molecules.
            
        Returns
        -------
        None : NoneType
        
        """
        
        # This function places the individual building block molecules
        # into a single ``.mol`` file. These molecules should be placed
        # on a given set of vertices or edges, depending on the topology
        # desired. Bonds are then created between the placed molecules.
        # This is done using the heavy atoms as identifiers. As a 
        # result, this creates the heavy atom substituted version of the 
        # cage. To produce the pristine verion of the cage, 
        # ``final_sub`` is called which replaces the heavy atoms with
        # their pristine counterparts / functional groups.
        self.place_mols()
        self.join_mols()
        self.final_sub()

    def join_cage_mols(self):
        """
        takens an input file with disconnected moleucles and connects them
        """           
        
        
        mol_file_data = Atom.extract_mol_file_data(
                                          self.macro_mol.heavy_mol_file)

        mol_file_content = mol_file_data.mol_file_content        
        
        Atom.bb_heavy_atoms_per_molecule = self.heavy_atoms_per_bb  
        Atom.lk_heavy_atoms_per_molecule = self.heavy_atoms_per_lk   
        
        cage_bb = next(x for x in self.macro_mol.building_blocks if
                                    isinstance(x, BuildingBlock))
        cage_lk = next(x for x in self.macro_mol.building_blocks if
                                    isinstance(x, Linker))
        
        Atom.bb = cage_bb.func_grp.heavy_symbol
        Atom.lk = cage_lk.func_grp.heavy_symbol
        
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
                
            mol_file_content += "M  V30 {2} {3} {0} {1}\n".format(
                              atom[0], atom[2], 
                              mol_file_data.bond_number, bond_order)
        
        
        mol_file_content += "M  V30 END BOND\nM  V30 END CTAB\nM  END"

        p = re.compile(r" VAL=.")
        mol_file_content = re.sub(p, "", mol_file_content)

        mol_file_content = mol_file_content.replace(
                                            mol_file_data.count_line, 
                                "M  V30 COUNTS {0} {1} 0 0 0\n".format(
                                             mol_file_data.at_num,
                                             mol_file_data.bond_number))
        
        new_mol_file_name = self.macro_mol.heavy_mol_file
        new_mol_file = open(new_mol_file_name, "w")
        new_mol_file.write(mol_file_content)
        new_mol_file.close()

    def join_polymer_mols(self):
        pass
        
        
    def final_sub(self):
        """
        Replaces heavy atoms with functional group atoms they represent.        
        
        Once a heavy cage has been assembled the pristine cage is formed
        by replacing the heavy atoms. This function does this by editing
        the heavy ``.mol`` file and changing the appropriate symbols.
        
        Modifies
        --------
        .mol file
            This function creates a ``.mol`` file holding the pristine
            assembled cage.
        
        Returns
        -------
        None : NoneType
        
        """

        cage_bb = next(x for x in self.macro_mol.building_blocks if
                                    isinstance(x, BuildingBlock))
        cage_lk = next(x for x in self.macro_mol.building_blocks if
                                    isinstance(x, Linker))

        
        # Read the assembled heavy cage file, line by line. Any time
        # a heavy atomic symbol is found, replace it with the 
        # corresponding light element. Save all lines to a string and
        # write this string to a new file. Do not change the original
        # file because it is meant to hold the heavy cage. 
        with open(self.macro_mol.heavy_mol_file, "r") as add:
            new_file= ""
            for line in add:
                line = line.replace(cage_bb.func_grp.heavy_symbol,
                                    cage_bb.func_grp.target_symbol)
                line = line.replace(cage_lk.func_grp.heavy_symbol,
                                    cage_lk.func_grp.target_symbol)
                new_file += line
    
        with open(self.macro_mol.prist_mol_file, "w") as f:
            f.write(new_file)

      
class FourPlusSix(Topology):
    """
    Defines the tetrahedral, 4+6, topology.

    This is a topology of cages where 4 building-blocks* are placed on
    vertices and 6 linkers are placed on the edges between them. This
    class defines functions which places these molecules in the correct
    positions in the ``.mol`` file.

    Attributes
    ----------
    This class also inhertis all the attributes of the ``Topology`` 
    class.
    
    heavy_atoms_per_bb : int (default=3, do not change)
        The number of heavy atoms (functional groups) in a 
        building-block* molecule used for this topolgy.
    
    heavy_atoms_per_lk : int (default=2, do not change)
        The number of heavy atoms (functional groups) in a linker 
        molecule used for this topology.
   
    bb_num : int (default=4, do not change)
        The number of building-block* molecules used in this topology.
        
    lk_num : int (default=6, do not change)
        The number of linker molecules used in this topology.
   
    pair_up_func : function (default=Atom.pair_up_v4_v2, do not change)
        The function used to find atoms in different building-block* and
        linker molecules which need to have a bond created between them.
   
    """
    
    def __init__(self, macro_mol):
        super().__init__(macro_mol)
        self.heavy_atoms_per_bb = 3
        self.heavy_atoms_per_lk = 2
        
        self.bb_num = 4
        self.lk_num = 6
        
        self.pair_up_func = Atom.pair_up_v4_v2
        self.join_mols = self.join_cage_mols
        
    def place_mols(self):
        """
        Places all building block molecules on correct coordinates.

        The building block molecules are placed in their appropraite 
        positions on the topology. This means that the building-blocks* 
        are placed on vertices and linkers on edges. This function
        only places the molecules, it does not join them. It saves the 
        structure to a ``.mol`` file and an rdkit molecule instance for
        later use. Both are placed in the appropriate attribtes of the 
        ``Cage`` instance the topology is describing.
        
        Modifies
        --------
        self.macro_mol.heavy_mol : 
        
        .mol file
            
        
        """
        
        # Building-blocks* and linkers have different rules for 
        # placement so each is treated by its own function.
        bb_placement = self.place_bbs()
        lk_placement = self.place_lks()
        self.macro_mol.heavy_mol = chem.CombineMols(bb_placement, 
                                                   lk_placement) 
                                                   
        chem.MolToMolFile(self.macro_mol.heavy_mol, 
                          self.macro_mol.heavy_mol_file,
                          includeStereo=True, kekulize=False,
                          forceV3000=True)
        
    
    def place_bbs(self):
        """
        Places all building-blocks* on their appropriate coordinates.
        
        This function creates a new rdkit molecule made up of 4
        building-block* molecules. Each of the building-block* molecules
        is placed on the vertex of a tetrahedron. The molecules are held
        within the same rdkit molecule instance but are otherwise 
        unconnected.
        
        Returns
        -------
        rdkit.Chem.rdchem.Mol
            An rdkit molecule which is made up of the building-block*
            molecules all placed on the vertices of a tetrahedron. The
            individual molecules are unconnected.
        
        """
        cage_bb = next(x for x in self.macro_mol.building_blocks if 
                                        isinstance(x, BuildingBlock))

        position1 = cage_bb.shift_heavy_mol(0,50,0)
        position2 = cage_bb.shift_heavy_mol(0,0,25)
        position3 = cage_bb.shift_heavy_mol(25,0,-25)
        position4 = cage_bb.shift_heavy_mol(-25,0,-25)

        combined_mol = chem.CombineMols(position1, position2)
        combined_mol2 = chem.CombineMols(position3, position4)
        return chem.CombineMols(combined_mol, combined_mol2)       
        
    def place_lks(self):
        """
        Places all linkers on their appropriate coordinates.

        This function creates a new rdkit molecule made up of 6 linker
        molecules. Each of the linker molecules is placed on the edge of 
        a tetrahedron. The molecules are held within the same rdkit 
        molecule instance but are otherwise unconnected.
        
        Returns
        -------
        rdkit.Chem.rdchem.Mol
            An rdkit molecule which is made up of the linker molecules 
            all placed on the edges of a tetrahedron. The individual 
            molecules are unconnected.
        
        """
        
        cage_lk = next(x for x in self.macro_mol.building_blocks if 
                                                isinstance(x, Linker)) 
                                                
        position1 = cage_lk.shift_heavy_mol(0, 0, -25)
        position2 = cage_lk.shift_heavy_mol(-12.5, 0, 0)
        position3 = cage_lk.shift_heavy_mol(12.5, 0, 0)
        
        position4 = cage_lk.shift_heavy_mol(-12.5, 25, -12.5)
        position5 = cage_lk.shift_heavy_mol(12.5, 25, -12.5)
        position6 = cage_lk.shift_heavy_mol(0, 25, 25)

        combined_mol = chem.CombineMols(position1, position2)
        combined_mol2 = chem.CombineMols(position3, position4)        
        combined_mol3 = chem.CombineMols(position5, position6)
        
        combined_mol4 = chem.CombineMols(combined_mol, combined_mol2)
        return chem.CombineMols(combined_mol3, combined_mol4)
        
class BlockCopolymer(Topology):
    def __init__(self, macro_mol, repeating_unit):
        super().__init__(macro_mol)
        self.repeating_unit = repeating_unit
        self.join_mols = join_polymer_mols
        
    def place_mols(self):
        pass
    
    
        
from MMEA.classes.mol_reader import Atom        