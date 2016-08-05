import rdkit
from rdkit import Chem as chem
from rdkit.Chem import AllChem as ac
import networkx as nx
from functools import partial

from MMEA.classes.molecular import FGInfo, BuildingBlock, Linker
from MMEA.convenience_functions import flatten

class Topology:
    """
    Represents the topology of an assembled molecule.
    
    The ``Topology`` class is concerned with how individual building 
    blocks are placed and connected in space to form an assembled 
    molecule used by MMEA. It also takes care of assembling these
    molecules. The general process of building macromolecules is 
    discussed in detail in the `build` method documentation.
    
    This class directly defines any operations and attributes that are
    needed by any topology, be it a tetrahedron, octahedron or even a
    polymer. However, this class is not used directly by MMEA. It is
    intended to be inherited from. Any individual within MMEA will have
    a `topology` attribute which refers to an instance of a class 
    derived from this. Derived classes of ``Topology`` define things
    specific to that one topology. For example, each derived class must 
    define which ``pair_up`` function defined in the ``Topology`` class 
    it uses for pairing up its heavy atoms. This is done by placing
    the function in the `pair_up` attribute of the derived class. See
    the included derived classes as examples. In addition, each class 
    derived from ``Topology`` must define methods which place building 
    blocks in the correct positions, such as chosen edges or vertices.

    Instances of this class should not be created directly. Only via a
    derived class.
    
    Attributes
    ----------
    macro_mol : MacroMolecule
        The ``MacroMolecule`` instance which has this topology. This 
        gives easy access to the macromolecule's attributes to the 
        ``Topology``instance.
        
    paired : set of ints
        This attribute is created and used during assembly by some pair
        up functions. Not all topolgies will need to do this. It is a
        set of atom ids which have already had a bond added to them
        during assembly.
        
    paired_mols : set of tuples of ints
        This attribute is created and used during assembly by some pair
        up functions. Not all topologies will need to do this. It is a 
        set of tuples. The tuples are sorted so that (1,2) and (2,1) are
        added as the same pairing. Note that using ``sorted`` on tuples
        returns a list, so it must be reconverted to a tuple before
        being added to the set. This is because lists are not hashable
        and cannot be added to sets a result.
        
    """

    
    def __init__(self, macro_mol):
        self.macro_mol = macro_mol
        

    def build(self):
        """
        Creates rdkit instances the heavy and pristine macromolecules.
        
        This function also places the created rdkit instances in the
        `prist_mol` and `heavy_mol` attributes of `self.macro_mol`.
        `self.macro_mol` is the ``MacroMolecule`` instance holding the 
        ``Topology`` instance carrying out `build`.
        
        To carry out `build` an instance of a class derived from 
        ``Topology`` must be used. This is because instances of such
        classes define a `pair_up` attribute during initialization.
        (This should be done by default, not passed as an argument to
        the initializer.) The `pair_up` attribute holds the pair up 
        function defined within ``Topology``, which should be used by
        `build`. (`pair_up` is used within the `join_mols` subroutine of
        `build`.)
        
        Modifies
        --------
        self.macro_molecule.heavy_mol
            Adds an rdkit instance of the heavy assembled molecule to 
            this attribute.
            
        self.macro_molecule.prist_mol
            Adds an rdkit instance of the pristine assembled molecule to
            this attribute.
            
        Returns
        -------
        None : NoneType
        
        """
        
        # This function places the individual building block molecules
        # into a single rdkit molecule instance. These molecules should 
        # be placed on a given set of vertices or edges, depending on 
        # the topology desired. Bonds are then created between the 
        # placed molecules. This is done using the heavy atoms as 
        # identifiers. As a result, this creates the heavy atom s
        # substituted version of the macromolecule. To produce the 
        # pristine verion of the molecule, `final_sub` is called which 
        # replaces the heavy atoms with their pristine counterparts / 
        # functional groups.
        self.place_mols()
        self.join_mols()
        self.final_sub()

    def join_mols(self):
        """
        Joins the disconnected building blocks of a macromolecule.
        
        Before this function is called the ``MacroMolecule`` instance 
        which holds a given topology, `self`, should have an rdkit 
        molecule instance in its `heavy_mol` attribute. At this point,
        the macromocule in this rdkit instance will consist of various
        building blocks on the edges and vertices of the defined 
        topology. This function joins them up into a single molecule.
        
        This function should not be altered when extending MMEA. It is
        designed so that it automatically calls the functions 
        appropriate for a given topology. See `build` method 
        documentation for more details on how this works.
        
        Modifies
        --------
        self.macro_molecule.heavy_mol
            The ``MacroMolecule`` instance which holds `self` has the 
            joined up rdkit version of the heavy macromolecule added to 
            its `heavy_mol` attribute.
        
        """
        
        # Get a mathematical graph representing the disconnected 
        # heavy molecule.
        heavy_graph = self.macro_mol.get_heavy_as_graph()
        
        # Use the graph to generate a list of lists. Each sublist is
        # a collection of atom ids all belonging to the same molecule.
        # x.nodes() generates the atom ids for each of the sublists. The
        # function ``connected_component_subgraphs`` generates the 
        # subgraphs representing the disconnected molecules.           
        molecules = [x.nodes() for x in 
                    nx.connected_component_subgraphs(heavy_graph)]
        
        # Get rid of all the atom ids which do not belong to heavy 
        # molecules, but keep data organised into sublists representing
        # molecules.
        heavy_mols = self.extract_heavy_atoms(molecules)
        
        # Run the pair up function chosen for the given topology. This
        # is set in the derived classes initializer. Only the atom ids
        # of heavy atoms are provided as bonds are only made between 
        # these. Keeping the data organised into sublists representing
        # molecules prevents bonds from being created between heavy
        # atoms on the same molecule by the pair up functions.
        self.pair_up(heavy_mols)            


    def extract_heavy_atoms(self, molecules):
        """
        Returns atom ids of heavy atoms, grouped by molecule.
        
        Parameters
        ----------
        molecules : iterable of iterables of ints
            This is list of the form [[1,2], [3,4,5], [9,8,10], [7,6]].
            Each sublist represents the atom ids belonging to the same
            building block molecule.
        
        Returns
        -------
        list of lists of ints
            The returned list has the form [[1,2], [3,5,7], [9,13]]. 
            Each sublist represents the atom ids belonging to the same
            building block molecule. Only heavy atom ids are present in
            the returned list.
        
        """
        
        # In essence, create a new sublist for each original sublist in 
        # in `molecules`. Copy atom ids from the original sublists into
        # the new subslists only if the atom ids belong to heavy
        # molecules.        
        
        
        # Create a new list which will holds the sublists holding heavy
        # atom ids. 
        heavy_mols = []
        
        # Iterate through each sublist in `molecules`. For each such 
        # sublist, ``molecule``, create a list, ``heavy_mol``, which 
        # will hold the atom ids found in ``molecule`` if the atom ids
        # belong to a heavy atom. 
        
        # Once all the atom ids in ``molecule`` have been checked, add
        # the created ``heavy_mol`` variable to ``heavy_mols``. Once
        # all sublists in ``molecules`` have been checked, return 
        # ``heavy_mols``.            
        for molecule in molecules:
            heavy_mol = []
            for atom_id in molecule:
                # Get the corresponding rdkit atom instance.
                atom = self.macro_mol.heavy_mol.GetAtomWithIdx(atom_id)              
                atom_n = atom.GetAtomicNum()  
                
                # Checks that the atom is heavy.
                if atom_n in FGInfo.heavy_atomic_nums:
                    # Heavy ids get added to the new sublist.
                    heavy_mol.append(atom_id)
            
            # After going through all the atom ids in the original 
            # sublist, add the new sublist to the collection of heavy
            # id sublists.
            heavy_mols.append(heavy_mol) 
        
        return heavy_mols
        
    def final_sub(self):
        """
        Replaces heavy atoms with functional group atoms they represent.        
        
        Once a heavy cage has been assembled the pristine macromolecule 
        is formed by replacing the heavy atoms and adding Hydrogens 
        where appropriate.
        
        Modifies
        --------
        self.macro_molecule.prist_mol
            Creates this attribute. It holds the rdkit instances of the
            assembled pristine macromolecule.
        
        Returns
        -------
        None : NoneType
        
        """

        # Add Hydrogens to the pristine version of the molecule. The 
        # ``GetSSSR`` function and optimization ensure that the added 
        # Hydrogen atoms are placed in reasonable positions. The 
        # ``GetSSSR`` function itself is just prerequisite for
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
        """
        Pairs up atoms of different elements in different molecules.
        
        For each atom this function finds the closest heavy atom which
        statisfies the following criteria:
            > it has not had a bond added to it yet
            > it has a different element
            > it is on a different molecule
            > the molecules of the two atoms are not bonded together
        
        Parameters
        ----------
        heavy_mols : iterable of iterables of ints
            This is list of the form [[1,2], [3,4,5], [9,8,10], [7,6]].
            Each sublist represents the atom ids of heavy atoms 
            belonging to the same building block molecule.
        
        Modifies
        --------
        self.macro_mol.heavy_mol
            Adds bonds between atoms in the rdkit molecule instance held
            by this attribute.
        
        Returns
        -------
        None : NoneType
        
        """
        
        # Create a ``EditableMol`` which allows for the adding of bonds.
        editable_mol = chem.EditableMol(self.macro_mol.heavy_mol)
        
        # Keep track of which atoms and molecules have been joined.
        self.paired = set()
        self.paired_mols = set()
        
        # Go through all the atom ids of heavy atoms and if that id 
        # belongs to an atom that has not yet been paired, generate a 
        # pool of atom ids of possible partners for pairing. The pool
        # consists of atom ids that fulfil the 4 criteria outlined in
        # the docstring. From these the atom closest to the atom being
        # iteratred through is chosen and a bond between the two is 
        # made. The sets tracking pairings are updated.
        for atom_id in flatten(heavy_mols):
            
            # Atom is not yet paired.
            if atom_id not in self.paired:
                
                # Generate the potential partner ids.
                partner_pool = self.unpaired_diff_element_atoms(atom_id, 
                                                         heavy_mols)
                
                # Find the closest atom from the partner pool.                                       
                partner = self.min_distance_partner(atom_id, 
                                                    partner_pool)
                
                # Check if the bond created should be double or single.
                bond_type = self.determine_bond_type(atom_id, partner)
                # Add the bond.                
                editable_mol.AddBond(atom_id, partner, bond_type)
                
                # Find that molecule number of the molecule that atom
                # and its partner belong to. The molecule number is just
                # the index of the sublist containing the id within 
                # `heavy_mols`.
                atom_mol_num = next(heavy_mols.index(x) for x in 
                                            heavy_mols if atom_id in x)
                
                partner_mol_num = next(heavy_mols.index(x) for x in 
                                            heavy_mols if partner in x)                 
                
                # Update the pair tracking sets.
                self.paired.add(atom_id)
                self.paired.add(partner)
                self.paired_mols.add(str(sorted((atom_mol_num, 
                                                 partner_mol_num))))
        
        # Turn the ``EditableMol`` into an rdkit molecule instance and
        # place that into the ``MacroMolecule``'s attribute.        
        self.macro_mol.heavy_mol = editable_mol.GetMol()

    def pair_up_polymer(self, heavy_mols):
        """
        Pairs up monomer units of a polymer.

        The monomer units should be placed in a straight line. This
        functions creates bonds between heavy atoms of any kind as long
        they are in separate molecules and the two molecules are not
        already joined. This function creates the shortest possible 
        bonds which satisfy this criteria first. It only creates N-1 
        bonds where N is number of monomers joined up in a repeating 
        unit of the polymer. For example if the repeating unit is
        ``AAA`` or ``ABC`` then N = 3 in both cases, and 2 bonds will be
        formed.
        
        Parameters
        ----------
        heavy_mols : iterable of iterables of ints
            This is list of the form [[1,2], [3,4,5], [9,8,10], [7,6]].
            Each sublist represents the atom ids of heavy atoms 
            belonging to the same building block molecule.

        Modifies
        --------    
        self.macro_mol.heavy_mol
            Adds bonds between atoms in the rdkit molecule instance held
            by this attribute.
        
        """
        
        # Create a ``EditableMol`` which allows for the adding of bonds.        
        editable_mol = chem.EditableMol(self.macro_mol.heavy_mol)
        
        # Keep track of which atoms and molecules have been joined.
        self.paired = set()
        self.paired_mols = set()
        
        # Get the distances between all heavy atoms. Starting with the 
        # shortest distance, check if the atoms are on the same molecule
        # is yes move on to the next shortest distance and try again. If
        # not, create a bond between the two atoms. If the number of
        # bonds created is greater than or equal to N - 1 stop.        
        distances = sorted(self.macro_mol.get_heavy_atom_distances())
        num_bonds_made= 0
        for _, atom1_id, atom2_id in distances:
            
            # The molecule number is just the index of the sublist 
            # holding the atom id in `heavy_mols`.           
            atom1_mol =  next(heavy_mols.index(x) for x in 
                                            heavy_mols if atom1_id in x)           
            
            atom2_mol =  next(heavy_mols.index(x) for x in 
                                            heavy_mols if atom2_id in x)
            mol_pair = str(sorted((atom1_mol, atom2_mol)))
            
            # Make sure the atoms or molecules are not already paired.            
            atom1_not_paired = atom1_id not in self.paired
            atom2_not_paired = atom2_id not in self.paired
            mols_not_paired = mol_pair not in self.paired_mols       
            # Check that you're not dealing with the same atom.            
            not_same_atom = atom1_mol != atom2_mol
            if (atom1_not_paired and atom2_not_paired and 
                 mols_not_paired and not_same_atom):
                    
                    # Check if double or single bond should be added.
                    bond_type= self.determine_bond_type(atom1_id, 
                                                        atom2_id)
                    # Add the bond.                    
                    editable_mol.AddBond(atom1_id, atom2_id, bond_type)
                    
                    # Update the tracking sets.                    
                    self.paired.add(atom1_id)
                    self.paired.add(atom2_id)
                    self.paired_mols.add(mol_pair)
                    
                    # Increment the bond count.
                    num_bonds_made += 1           
                    # Break if enough bonds have been made.
                    if num_bonds_made >= len(self.repeating_unit)-1:
                        break

        # Turn the ``EditableMol`` into an rdkit molecule instance and
        # place that into the ``MacroMolecule``'s attribute.        
        self.macro_mol.heavy_mol = editable_mol.GetMol()



    def unpaired_diff_element_atoms(self, atom1_id, heavy_mols):
        """
        Yield atoms forming the partner pool of another atom.        
        
        The atom ids yielded from `heavy_mols` belong to atoms which 
        satisfy the following criteria:
            > have a different element to the atom with `atom1_id`
            > are in a different molecule to the atom with `atom1_id`
            > have not yet been bonded to another heavy atom
            > belong to a molecule which has not been bonded to the 
              molecule containing the atom with `atom1_id`
                
        Parameters
        ----------
        atom1_id : int
            The id of the atom whose partner pool should be generated.
            
        heavy_mols : iterable of iterables of ints
            This is list of the form [[1,2], [3,4,5], [9,8,10], [7,6]].
            Each sublist represents the atom ids of heavy atoms 
            belonging to the same building block molecule.
        
        Yields
        ------
        int
            The id of the next atom within the partner pool.
        
        """
        
        for atom2_id in flatten(heavy_mols):
            # Get the rdkit atom instance with the given atom ids.
            atom1 = self.macro_mol.heavy_mol.GetAtomWithIdx(atom1_id)
            atom2 = self.macro_mol.heavy_mol.GetAtomWithIdx(atom2_id)

            # Get the molecule numbers of the molecules holding atoms
            # with those ids. The molecule number is just the index of 
            # the sublist holding the atom id in `heavy_mols`.             
            atom1_mol =  next(heavy_mols.index(x) for x in heavy_mols 
                                                       if atom1_id in x)           
            atom2_mol =  next(heavy_mols.index(x) for x in heavy_mols 
                                                       if atom2_id in x)
            mol_pair = str(sorted((atom1_mol, atom2_mol)))
            
            # Check if the conditions described in the docstring are
            # satisfied. If yes, yield the atom.
            if (atom1.GetAtomicNum() != atom2.GetAtomicNum() and 
                atom2_id not in self.paired and 
                mol_pair not in self.paired_mols and 
                atom1_id != atom2_id):
                
                yield atom2_id
              
    def min_distance_partner(self, atom_id, partner_pool):
        """
        Return the closest atom from `partner_pool`.
        
        A number of atom ids is supplied and this function finds the
        atom id of the atom closest to the atom who's id is supplied in 
        `atom_id`.

        Parameters
        ----------
        atom_id : int
            The id of an atom whose minimum distance partner needs to be
            found.
            
        partner_pool : iterable of iterables of ints
        
        Returns
        -------
        int
            An id from `partner_pool`. It belongs to the atom which is  
            shortest distance away from the atom, whose id was
            supplied in `atom_id`.
        
        """
        
        # ``partial`` creates a function which can be used as a key by
        # the ``min`` function. A function can only be used as a key if 
        # it takes a single argument. ``partial`` fill the first 
        # argument of `self.macro_mol.heavy_distance` with `atom_id`
        # which means that the output of the min function is the atom id
        # of the atom which is the closest to `atom_id`.
        distance_func = partial(self.macro_mol.heavy_distance, atom_id)
        return min(partner_pool, key=distance_func)

    def determine_bond_type(self, atom1_id, atom2_id):
        """
        Returns the bond order to be formed between the atoms.
        
        Some atoms will need to have a double bond created between them.
        This is defined in the `FGInfo.double_bond_combs` list. If the
        atom ids provided as paramters belong to atoms of elements
        found in this list, the rdkit double bond type will be returned.
        If not the rdkit single bond type will be returned. These types
        are needed when adding bonds using ``EditableMol`` instances.
        
        Parameters
        ----------
        atom1_id : int
            The id number of the first atom.
        
        atom2_id : int
            The id number of the second atom.
        
        Returns
        -------
        rdkit.Chem.rdchem.BondType.SINGLE
            If the combination of heavy atoms passed as arguments is not 
            in `FGInfo.double_bond_combs`.
        
        rdkit.Chem.rdchem.BondType.DOUBLE
            If the combination of heavy atoms passed as arguments is in
            `FGInfo.double_bond_combs`.
            
        """
        
        # Get the atomic numbers of the of the atoms whose atom ids were
        # supplied as arguments. Then use `FGInfo.functional_group_list`
        # attribute to find the atomic symbols. If the atomic symbols
        # for ma tuple in `FGInfo.double_bond_combs` return a rdkit 
        # double bond type. If they do not, return a rdkit single bond 
        # type.
        
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
        
        double_bond_present = ((atom1_symbol, atom2_symbol) == tup or 
                               (atom2_symbol, atom1_symbol) == tup for 
                               tup in FGInfo.double_bond_combs)
        
        if True in double_bond_present:
            return rdkit.Chem.rdchem.BondType.DOUBLE
        else:
            return rdkit.Chem.rdchem.BondType.SINGLE

class FourPlusSix(Topology):
    """
    Defines the tetrahedral, 4+6, topology.

    This is a topology of cages where 4 building-blocks* are placed on
    vertices and 6 linkers are placed on the edges between them. This
    class defines functions which place these molecules in the correct
    positions within an rdkit instance. The rdkit instance is stored in 
    the `heavy_mol` attribute of a ``Cage`` instance.

    Attributes
    ----------
    This class also inhertis all the attributes of the ``Topology`` 
    class. Only the attribute `pair_up` must be defined, which defines
    which ``pair up`` inhertied from ``Topology`` should be used for
    pairing up the linkers and building-blocks*. This attribute is
    default initialized and cannot be set during initialization or run
    time of MMEA. It is hard coded.
       
    pair_up : function (default = self.pair_up_diff_element_atoms)
        The function used to find atoms in different building-block* and
        linker molecules which need to have a bond created between them.
        It also creates the bonds between them.
        
    """
    
    def __init__(self, macro_mol):
        super().__init__(macro_mol)        
        self.pair_up = self.pair_up_diff_element_atoms

    def place_mols(self):
        """
        Places all building block molecules on correct coordinates.

        The building block molecules are placed in their appropriate 
        positions based on the topology. This means that the 
        building-blocks* are placed on vertices and linkers on edges. #
        This function only places the molecules, it does not join them. 
        It saves the structure a rdkit molecule instance for later use. 
        This rdkit instace is placed in the `heavy_mol` attribute of the 
        ``Cage`` instance the topology is describing.
        
        Modifies
        --------
        self.macro_mol.heavy_mol
            Places an rdkit instance with disconnected building blocks
            placed on edges and vertices in this attribute.    
        
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
 
class EightPlusTwelve(FourPlusSix):
    
    def place_bbs(self):
        cage_bb = next(x for x in self.macro_mol.building_blocks if 
                                        isinstance(x, BuildingBlock))

        position1 = cage_bb.shift_heavy_mol(0,0,0)
        position2 = cage_bb.shift_heavy_mol(50,0,0)
        position3 = cage_bb.shift_heavy_mol(0,50,0)
        position4 = cage_bb.shift_heavy_mol(50,50,0)

        combined_mol = chem.CombineMols(position1, position2)
        combined_mol2 = chem.CombineMols(position3, position4)
        one = chem.CombineMols(combined_mol, combined_mol2)   
        
        position1 = cage_bb.shift_heavy_mol(0,0,50)
        position2 = cage_bb.shift_heavy_mol(50,0,50)
        position3 = cage_bb.shift_heavy_mol(0,50,50)
        position4 = cage_bb.shift_heavy_mol(50,50,50)

        combined_mol = chem.CombineMols(position1, position2)
        combined_mol2 = chem.CombineMols(position3, position4)
        two = chem.CombineMols(combined_mol, combined_mol2)          
        
        return chem.CombineMols(one, two)
        
    def place_lks(self):
        cage_lk = next(x for x in self.macro_mol.building_blocks if 
                                        isinstance(x, Linker))

        position1 = cage_lk.shift_heavy_mol(25,0,0)
        position2 = cage_lk.shift_heavy_mol(0,25,0)
        position3 = cage_lk.shift_heavy_mol(25,50,0)
        position4 = cage_lk.shift_heavy_mol(50,25,0)

        combined_mol = chem.CombineMols(position1, position2)
        combined_mol2 = chem.CombineMols(position3, position4)
        one = chem.CombineMols(combined_mol, combined_mol2)   
        
        position1 = cage_lk.shift_heavy_mol(25,0,50)
        position2 = cage_lk.shift_heavy_mol(0,25,50)
        position3 = cage_lk.shift_heavy_mol(25,50,50)
        position4 = cage_lk.shift_heavy_mol(50,25,50)

        combined_mol = chem.CombineMols(position1, position2)
        combined_mol2 = chem.CombineMols(position3, position4)
        two = chem.CombineMols(combined_mol, combined_mol2)          

        three = chem.CombineMols(one, two) 

        position1 = cage_lk.shift_heavy_mol(0,0,25)
        position2 = cage_lk.shift_heavy_mol(50,0,25)
        position3 = cage_lk.shift_heavy_mol(0,50,25)
        position4 = cage_lk.shift_heavy_mol(50,50,25)

        combined_mol = chem.CombineMols(position1, position2)
        combined_mol2 = chem.CombineMols(position3, position4)
        four = chem.CombineMols(combined_mol, combined_mol2)         
        
        
        return chem.CombineMols(three, four)
       
class BlockCopolymer(Topology):
    """
    A class for describing the repeating units of polymers.    
    
    """
    
    keys = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"    
    
    def __init__(self, macro_mol, repeating_unit):
        super().__init__(macro_mol)
        self.repeating_unit = repeating_unit
        self.pair_up = self.pair_up_polymer
        self.monomer_keys = {}
        for key, monomer in zip(BlockCopolymer.keys, 
                                self.macro_mol.building_blocks):
            self.monomer_keys[key] = monomer
    
    
    def place_mols(self):
        """
        Places monomer in a line, seperated by an equal distance.        
        
        """
        
        distance = 30
        self.macro_mol.heavy_mol = chem.Mol()
        for index, key in enumerate(self.repeating_unit):
            shifted_monomer = self.monomer_keys[key].shift_heavy_mol(
                                                    distance*index,0,0)
            self.macro_mol.heavy_mol = chem.CombineMols(shifted_monomer, 
                                            self.macro_mol.heavy_mol)
        



      