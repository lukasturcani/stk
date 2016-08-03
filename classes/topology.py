import rdkit
from rdkit import Chem as chem
from rdkit.Chem import AllChem as ac
import math
import numpy as np
import itertools as itertools
import re
import networkx as nx
from functools import partial
import heapq

from MMEA.classes.molecular import FGInfo, BuildingBlock, Linker
from MMEA.convenience_functions import flatten

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
    
    Attributes
    ----------
    macro_mol : MacroMolecule
        The macromolecule instance which has the given topology. This 
        provides easy access to the cages attributes to the ``Topology`` 
        instance.
    
    """

    
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



    def join_mols(self):
        heavy_graph = self.macro_mol.get_heavy_as_graph()
        molecules = [x.nodes() for x in 
                    nx.connected_component_subgraphs(heavy_graph)]
        
        heavy_mols = self.extract_heavy_atoms(molecules)
        self.pair_up(heavy_mols)            


    def extract_heavy_atoms(self, molecules):
        heavy_mols = []
        for molecule in molecules:
            heavy_mol = []
            for atom_id in molecule:
                if self.macro_mol.heavy_mol.GetAtomWithIdx(atom_id).GetAtomicNum() in FGInfo.heavy_atomic_nums:
                    heavy_mol.append(atom_id)
            heavy_mols.append(heavy_mol) 
        return heavy_mols


        
        
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
        
        
        # Add Hydrogens to the pristine version of the molecule and
        # ensure this updated molecule is added to the ``.mol`` file as 
        # well. The ``GetSSSR`` function and optimization ensure that
        # the added Hydrogen atoms are placed in reasonable positions.
        # The ``GetSSSR`` function itself is just prerequisite for
        # running the optimization.
        self.macro_mol.prist_mol = chem.Mol(self.macro_mol.heavy_mol)
        
        for atom in self.macro_mol.prist_mol.GetAtoms():
            atomic_num = atom.GetAtomicNum()
            if atomic_num in FGInfo.heavy_atomic_nums:
                target_atomic_num = next(x.target_atomic_num for x in 
                                    FGInfo.functional_group_list if 
                                    x.heavy_atomic_num == atomic_num)                
                atom.SetAtomicNum(target_atomic_num)
            
        
        
        self.macro_mol.prist_mol = chem.AddHs(self.macro_mol.prist_mol)
        chem.GetSSSR(self.macro_mol.prist_mol)
        ac.MMFFOptimizeMolecule(self.macro_mol.prist_mol)

    def pair_up_diff_element_atoms(self, heavy_mols):
        editable_mol = chem.EditableMol(self.macro_mol.heavy_mol)
        
        self.paired = set()
        self.paired_mols = set()
        
        for atom_id in flatten(heavy_mols):
            if atom_id not in self.paired:
                
                partner_pool = self.unpaired_diff_element_atoms(atom_id, 
                                                         heavy_mols)
                partner = self.min_distance_partner(atom_id, 
                                                    partner_pool)
                
                                    
                bond_type = self.determine_bond_type(atom_id, partner)
                editable_mol.AddBond(atom_id, partner, bond_type)

                atom_mol_num = next(heavy_mols.index(x) for x in heavy_mols if atom_id in x)
                partner_mol_num = next(heavy_mols.index(x) for x in heavy_mols if partner in x)                 
                
                self.paired.add(atom_id)
                self.paired.add(partner)
                self.paired_mols.add(str(sorted((atom_mol_num, partner_mol_num))))
                
        self.macro_mol.heavy_mol = editable_mol.GetMol()

    def pair_up_polymer(self, heavy_mols):
        editable_mol = chem.EditableMol(self.macro_mol.heavy_mol)
        
        self.paired = set()
        self.paired_mols = set()
        
        distances = sorted(self.macro_mol.get_heavy_atom_distances())
        num_bonds_made= 0
        for _, atom1_id, atom2_id in distances:
            atom1_mol =  next(heavy_mols.index(x) for x in heavy_mols if atom1_id in x)           
            atom2_mol =  next(heavy_mols.index(x) for x in heavy_mols if atom2_id in x)
            mol_pair = str(sorted((atom1_mol, atom2_mol)))
            
            if (atom1_id not in self.paired and atom2_id not in self.paired and 
                mol_pair not in self.paired_mols and atom1_mol != atom2_mol):
                    bond_type= self.determine_bond_type(atom1_id, atom2_id)
                    editable_mol.AddBond(atom1_id, atom2_id, bond_type)
                    self.paired.add(atom1_id)
                    self.paired.add(atom2_id)
                    self.paired_mols.add(mol_pair)
                    
                    num_bonds_made += 1                    
                    if num_bonds_made >= len(self.repeating_unit)-1:
                        break
        
        self.macro_mol.heavy_mol = editable_mol.GetMol()
        chem.MolToMolFile(self.macro_mol.heavy_mol, 'hi.mol')
        raise Exception()
        

    def unpaired_diff_element_atoms(self, atom1_id, heavy_mols):
        for atom2_id in flatten(heavy_mols):
            atom1 = self.macro_mol.heavy_mol.GetAtomWithIdx(atom1_id)
            atom2 = self.macro_mol.heavy_mol.GetAtomWithIdx(atom2_id)
            
            atom1_mol =  next(heavy_mols.index(x) for x in heavy_mols if atom1_id in x)           
            atom2_mol =  next(heavy_mols.index(x) for x in heavy_mols if atom2_id in x)
            mol_pair = str(sorted((atom1_mol, atom2_mol)))
            
            if (atom1.GetAtomicNum() != atom2.GetAtomicNum() and 
                atom2_id not in self.paired and mol_pair not in self.paired_mols):
                yield atom2_id


                        
    def min_distance_partner(self, atom_id, partner_pool):
        distance_func = partial(self.macro_mol.heavy_distance, atom_id)
        return min(partner_pool, key=distance_func)

    def determine_bond_type(self, atom1_id, atom2_id):
        atom1 = self.macro_mol.heavy_mol.GetAtomWithIdx(atom1_id)
        atom1_atomic_n = atom1.GetAtomicNum()
        atom2 = self.macro_mol.heavy_mol.GetAtomWithIdx(atom2_id)
        atom2_atomic_n = atom2.GetAtomicNum()
                
        
        atom1_symbol = next(x.heavy_symbol for x in 
                            FGInfo.functional_group_list if 
                            atom1_atomic_n == x.heavy_atomic_num)
        atom2_symbol = next(x.heavy_symbol for x in 
                            FGInfo.functional_group_list if 
                            atom2_atomic_n == x.heavy_atomic_num)        
        
        double_bond_present = (atom1_symbol in tup and atom2_symbol in tup
                                for tup in FGInfo.double_bond_combs)
        if True in double_bond_present:
            return rdkit.Chem.rdchem.BondType.DOUBLE
        else:
            return rdkit.Chem.rdchem.BondType.SINGLE


      
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
        self.pair_up = self.pair_up_diff_element_atoms

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

    keys = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"    
    
    def __init__(self, macro_mol, repeating_unit):
        super().__init__(macro_mol)
        self.repeating_unit = repeating_unit
        self.pair_up = self.pair_up_polymer
        self.monomer_keys = {}
        for key, monomer in zip(BlockCopolymer.keys, self.macro_mol.building_blocks):
            self.monomer_keys[key] = monomer
    
    
    def place_mols(self):
        distance = 30
        self.macro_mol.heavy_mol = chem.Mol()
        for index, key in enumerate(self.repeating_unit):
            shifted_monomer = self.monomer_keys[key].shift_heavy_mol(distance*index,0,0)
            self.macro_mol.heavy_mol = chem.CombineMols(shifted_monomer, self.macro_mol.heavy_mol)
        
      