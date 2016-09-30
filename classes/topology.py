import rdkit
from rdkit import Chem as chem
import networkx as nx
from functools import partial
import numpy as np
import itertools
from collections import deque
from scipy.spatial.distance import euclidean

from .molecular import FGInfo, StructUnit3, StructUnit2
from ..pyWindow import window_sizes
from ..convenience_functions import (flatten, normalize_vector, 
                                     vector_theta, atom_vdw_radii,
                                     rotation_matrix_arbitrary_axis,
                                     LazyAttr)

class Vertex:
    """
    Used to represent the vertices of Cage polyhedra.

    This class stores information about the vertices which make up a 
    Cage's structure.
    
    Attributes
    ----------
    x : float
        The x position of the vertex.
    
    y : float
        The y position of the vertex.    
        
    z : float
        The z position of the vertex.    
    
    coord : numpy.array of floats
        A numpy array which holds the x, y and z coordinates of the
        vertex, in that order.
    
    edges : list of Edge instances
        This list holds the Edge instances which represent the edges
        connected to a partical vertex within a given cage structure.
    
    heavy_ids : list of ints
        This list holds the ids of the heavy atoms which belong to the
        building block placed on a particular vertex. The ids correspond
        to the id of the heavy atoms in the cage molecule. This means
        they correspond to the ids of atoms in the `heavy_mol` attribute
        of a ``Cage`` instance.
    
    """

    __slots__ = ['x', 'y', 'z', 'coord', 'edges', 'heavy_ids']    
    
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.coord = np.array([x,y,z])
        self.edges = []
        self.heavy_ids = []
        
    def place_mol(self, building_block):
        """
        Places a building-block* on the coords of the vertex.
        
        The orientation of the building-block* is aligned with 2
        parameters. Firstly the normal of the plane of heavy atoms of
        the building-block* is aligned with the normal of the plane
        formed by the edges connected to the vertex. Because the normal
        of the plane of heavy atoms always points in the direction of
        the building_block*'s centroid, this alignment causes the bulk
        of the building-block* molecule to point away from the center
        of the cage.
        
        Secondly, the building-block* is rotated so that a heavy atom is 
        aligned perfectly with an edge. This reduces the rms distance
        between the edges and heavy atoms to some extent.                
        
        Parameters
        ----------
        building_block : StructUnit3
            The building-block* molecule to be placed on a vertex.
        
        Modifies
        --------
        building_block.heavy_mol : rdkit.Chem.rdchem.Mol
            The conformer of the rdkit instance in this attribute is 
            modified as per the description in the docstring.
            
        Returns
        -------
        rdkit.Chem.rdchem.Mol
            The rdkit instance holding the building-block* molecule with
            the coordinates placed on the vertex and orientation set as
            described in the docstring.
        
        """
        
        # The method first aligns the normal of the heavy atom plane to
        # the normal of the edge plane. This means the bulk of the 
        # building-block* is always pointed away from the center of the
        # molecule.
        building_block.set_heavy_mol_orientation(
                                               self.edge_plane_normal())    
        
        # Next, the building-block* must be rotated so that one of the 
        # heavy atoms is perfectly aligned with one of the edges. This
        # is a multi-step process:
        #   1) Place the centroid of the heavy atoms at the origin.
        #   2) Place the centroid of the edges connected to the vertex 
        #      at the origin.
        #   3) Rotate the building-block* by some amount `theta`, so
        #      so that one of the heavy atoms is perfectly aligned with
        #      one of the edges. The axis of rotation is the normal to 
        #      the plane of heavy atoms.
        # 
        # The rotation is carried out via matrices. This means a
        # coordinate matrix of atoms in the heavy molecule is generated
        # and modified.
        
        # Set the centroid of the heavy atoms at the origin.
        building_block.set_heavy_atom_centroid([0,0,0])
        # Get the coordinate of the atom which is to be aligned with an
        # edge.
        atom_coord = building_block.atom_coords('heavy',
                                            building_block.heavy_ids[0])

        # Get the coordinates of all the edges and translate the 
        # centroid to the origin.
        edge_coord_mat = self.edge_coord_matrix() - self.edge_centroid()
        edge_coord = np.array(edge_coord_mat[0,:])[0]
        
        # Get the angle between an edge and the atom.        
        theta = vector_theta(edge_coord, atom_coord)
        
        # Get the rotation matrix necessary to do the rotation of 
        # `theta` about the normal to the plane.
        rot_mat = rotation_matrix_arbitrary_axis(theta, 
                                               self.edge_plane_normal())
        # Apply the rotation to the positions of the atoms in the heavy
        # molecule and get a new position matrix which holds their
        # coordinates coordinates after the rotation.
        pos_mat = building_block.position_matrix('heavy')
        new_pos_mat = np.dot(rot_mat, pos_mat)
        # Update the positions in the rdkit instance in `heavy_mol`.
        building_block.set_position_from_matrix('heavy', new_pos_mat)
        
        # Finally the well orientated building-block* is placed on the
        # coords of the vertex.
        building_block.set_heavy_atom_centroid(self.coord)
     
        return building_block.heavy_mol

    def edge_plane_normal(self):
        """
        Return the normal of the plane formed by the connected edges.
        
        The normal is set such that it always points away from the 
        origin.        
        
        Returns
        -------
        numpy.array
            A normalized vector which defines the normal pointed away
            from the origin.        
        
        """
        # Get two of the direction vectors running between the edges.
        v1, v2 = itertools.islice(self.edge_direction_vectors(), 2)  
        # To get the normal to the plane get the cross product of these
        # vectors. Normalize it.        
        normal = normalize_vector(np.cross(v1, v2))
        
        # To check that the normal is pointing away from the center of
        # cage, find the angle, `theta`, between it and one of the
        # position vectors of the edges on the plane. Assuming that the
        # center of the cage is at the origin, which it should be as 
        # this is specified in the documentation, if the angle between
        # the normal the position vector is less than 90 degrees they
        # point in the same general direction. If the angle is greater
        # than 90 degrees it means that they are pointing in opposite 
        # directions. If this is the case make sure to multiply the 
        # nomral by -1 in all axes so that it points in the correct 
        # direction while still acting as the normal to the plane.
        theta = vector_theta(normal, self.edges[0].coord) 
        
        if theta > np.pi/2:
            normal = np.multiply(normal, -1)
        
        return normal
    
    def edge_plane(self):
        """
        Return coefficients of plane of edges connected to the vertex.
        
        A plane is defined by the scalar plane equation,
            
            ax + by + cz = d.
        
        This method returns the a, b, c and d coefficients of this 
        equation for the plane formed by the connected edges. The 
        coefficents a, b and c decribe the normal vector to the plane.
        The coefficent d is found by substituting these coefficients
        along with the x, y and z variables in the scalar equation and
        solving for d. The variables x, y and z are substituted by the
        coordinate of some point on the plane. For example, the position
        of one of the heavy atoms.
        
        Returns
        -------
        numpy.array
            This array has the form [a, b, c, d] and represents the 
            scalar equation of the plane formed by the heavy atoms.
        
        References
        ----------
        http://tutorial.math.lamar.edu/Classes/CalcIII/EqnsOfPlanes.aspx  
        
        """
        
        heavy_coord = self.edges[0].coord
        d = np.multiply(np.sum(np.multiply(self.edge_plane_normal(), 
                                           heavy_coord)), -1)
        return np.append(self.edge_plane_normal(), d)
        
    def edge_direction_vectors(self):
        """
        Yields direction vectors between edges connected to the vertex.
        
        Yields
        ------
        numpy.array
            A normalized direction vector running from one edge 
            connected to the vertex to another.        
        
        """

        for edge1, edge2 in itertools.combinations(self.edges, 2):
            yield normalize_vector(edge1.coord-edge2.coord)
    
    def edge_coord_matrix(self):
        """
        Return matrix holding coords of edges joined to the vertex.        

        Returns
        -------
        numpy.matrix
            The matrix is n x 3, where n is the number of edges
            connected to the vertex. The row holds the x, y and z
            coordinates, respectively.
        
        """
        
        coords = []
        for edge in self.edges:
            coords.append(edge.coord)
        return np.matrix(coords)
        
    def edge_centroid(self):
        """
        Returns the centroid of the edges connected to the vertex.
        
        Returns
        -------
        numpy.array
            An array which holds the x, y and z positions of the
            centroid of the edges connected to the vertex.
        
        """
        
        # The connected edges are held in the `edges`. To get the
        # centroid, add up all the x, y and z coordinates (separately) 
        # and divide each sum by the number of edges. 
        return sum(edge.coord for edge in self.edges) / len(self.edges)
        

class Edge:
    """
    Used to represent the edges of Cage polyhedra.

    This class stores information about the edges which make up a Cage's 
    structure.
    
    Attributes
    ----------
    v1 : Vertex
        The first vertex which an edge is connected to.
    
    v2 : Vertex
        The second vertex which an edge is connected to.   
    
    coord : numpy.array of floats
        A numpy array which holds the x, y and z coordinates of the
        edge, in that order. It corresponds to the midpoint of the two
        vertices which the edge connects.
    
    direction : numpy.array
        This vector represents the orientation of the edge. It is a 
        normalized direction vector which runs from `v2` to `v1`.
    
    heavy_ids : list of ints
        This list holds the ids of the heavy atoms which belong to the
        building block placed on a particular edge. The ids correspond
        to the id of the heavy atoms in the cage molecule. This means
        they correspond to the ids of atoms in the `heavy_mol` attribute
        of a ``Cage`` instance.
    
    """
    
    __slots__ = ['v1', 'v2', 'coord', 'direction', 'heavy_ids']    
    
    def __init__(self, v1, v2):
        self.v1 = v1
        self.v2 = v2
        self.coord = np.divide(np.add(v1.coord, v2.coord), 2)
        self.direction = normalize_vector(v1.coord - v2.coord)
        v1.edges.append(self)
        v2.edges.append(self)
        
        self.heavy_ids = []
        
    def place_mol(self, linker):
        """
        Places a linker molecule on the coordinates of an edge.
        
        It also orientates the linker so that a the heavy atoms sit
        exactly on the edge and bulk of the linker points away from the 
        center of the cage.

        Parameters
        ----------
        linker : StructUnit2
            The linker which is to be placed and orientated as described
            in the docstring.
        
        Modifies
        --------
        linker.heavy_mol : rdkit.Chem.rdchem.Mol
            The conformer of the rdkit instance in this attribute is 
            modified as per the description in the docstring. 
       
        Returns
        -------
        rdkit.Chem.rdchem.Mol
            The rdkit instance holding the linker molecule with the
            coordinates placed on the edge and orientation set as
            described in the docstring.

        """

        # First the centroid of the heavy atoms is placed on the
        # position of the edge, then the direction of the linker is 
        # aligned with the direction of the edge.
        linker.set_heavy_atom_centroid(self.coord)
        linker.set_heavy_mol_orientation(np.multiply(self.direction,
                                         np.random.choice([1,-1])))

        # Ensure the centroid of the linker is placed on the outside of 
        # the cage.
        linker.minimize_theta(self.coord, self.direction)
        
        return linker.heavy_mol

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
    derived class. Multiple inheritance can be useful when creating 
    derived classes. For example, all topologies describing cages will
    share some characteristics. This means a class ``CageTopology`` can
    be created which holds all information required by all cage 
    topologies. This class, ``CageTopology``, will inherit ``Topology``. 
    A specific cage topology such as ``FourPlusSix`` or 
    ``EightPlusTwelve`` will then inherit ``CageTopology`` and add any 
    information specific to that one topology.
    
    Extending MMEA: Adding new topologies
    -------------------------------------
    > Cages
    To add a new cage topology a new class should be created, named
    after the topology. This class should inhertic the ``CageTopology``
    class. This will give access to various methods which are necessary
    for dealing with any cage molecule. See the documenation of 
    ``CageTopology`` for more details.
    
    The new class will only need to have five class attributes added:
        1) a list called `vertices` 
        2) a list called `edges`
        3) `n_windows` which holds the number of windows the cage 
           topology has
        4) `n_window_types` which holds the number of different window
           types. For example, if `n_window_types` is 2 then the 
           topology will have two kinds of windows, each with a 
           different expected size even in a perfectly symmetrical case. 
           Windows of the same type are expected to be of the same size.
        
    The `vertices` list holds instances of the class ``Vertex``. Each
    instance represents a vertex of a cage and needs to be initialized
    with the coordinates of that vertex. Vertices of a cage are where
    building-blocks* of cages are placed.
    
    The `edges` list holds instances of the class ``Edge``. Each
    instance represents an edge of a cage and needs to be initialized
    with two instnaces of the ``Vertex`` class. The ``Vertex`` instances
    should be held in the `vertices` list mentioned above. These are 
    the two vertices which the edge connects. Linkers of cages are 
    placed on edges. The edge instances automatically derive their 
    positions from the vertices supplied during initialization.

    The vertices need to be positioned such that the center of the
    topology is at the origin.
    
    Attributes
    ----------
    macro_mol : MacroMolecule
        The ``MacroMolecule`` instance which has this topology. This 
        gives easy access to the macromolecule's attributes to the 
        ``Topology``instance.
        
    bonds_made : int
        The number of bonds created during assembly. This should be
        incremened for each new bond made during assembly. Used in some
        fitness functions.
    
    paired : set of ints
        This attribute is created and used during assembly by some pair
        up functions. Not all topolgies will need to do this and as a 
        result not all instances will have this attribute. It is a set
        of atom ids which have already had a bond added to them during
        assembly.
        
    paired_mols : set of tuples of ints
        This attribute is created and used during assembly by some pair
        up functions. Not all topologies will need to do this and as a 
        result not all instnaces will have this attribute. It is a set 
        of tuples. The tuples are sorted so that (1,2) and (2,1) are
        added as the same pairing. Note that using ``sorted`` on tuples
        returns a list, so it must be reconverted to a tuple before
        being added to the set. This is because lists are not hashable
        and cannot be added to sets a result.
        
    """
    
    def __init__(self, macro_mol):
        self.macro_mol = macro_mol
        self.bonds_made = 0         
        
    def build(self):
        """
        Creates rdkit instances of heavy and pristine macromolecules.
        
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
#        chem.MolToMolFile(self.macro_mol.heavy_mol, '1heavy.mol')
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
        heavy_graph = self.macro_mol.graph('heavy')
        
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
        
        This function is only usable during assembly. After assembly
        all heavy atoms will be on the same molecule. This makes
        grouping the heavy atoms by molecule somewhat nonsensical.
        
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
        
        self.macro_mol.prist_mol = chem.Mol(self.macro_mol.heavy_mol)
        
        for atom in self.macro_mol.prist_mol.GetAtoms():
            atomic_num = atom.GetAtomicNum()
            if atomic_num in FGInfo.heavy_atomic_nums:
                target_atomic_num = next(x.target_atomic_num for x in 
                                    FGInfo.functional_group_list if 
                                    x.heavy_atomic_num == atomic_num)                
                atom.SetAtomicNum(target_atomic_num)
                atom.UpdatePropertyCache()
        
        # Updating the property cache recalculates valencies. This means
        # Hydrogen atoms are added in places where they are missing.
        self.macro_mol.prist_mol.UpdatePropertyCache()
        self.macro_mol.prist_mol = chem.AddHs(self.macro_mol.prist_mol,
                                              addCoords=True)


    def pair_up_edges_with_vertices(self, *args):


        editable_mol = chem.EditableMol(self.macro_mol.heavy_mol)
        
        for edge in self.edges:
            for atom_id in edge.heavy_ids:
                atom_coord = self.macro_mol.atom_coords('heavy', 
                                                        atom_id)                
                
                
                distance = partial(euclidean, atom_coord)
                vertex = min([edge.v1.coord, edge.v2.coord], 
                             key=distance)
                             
                vertex = next(x for x in [edge.v1, edge.v2] if np.array_equal(vertex, x.coord))
                partner = self.min_distance_partner(atom_id, 
                                                    vertex.heavy_ids)
                vertex.heavy_ids.remove(partner)
                
                bond_type = self.determine_bond_type(atom_id, partner)
                # Add the bond.                
                editable_mol.AddBond(atom_id, partner, bond_type)
                self.bonds_made += 1
                
        self.macro_mol.heavy_mol = editable_mol.GetMol()

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
        distance_func = partial(self.macro_mol.atom_distance, 
                                'heavy', atom_id)
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


    def place_mols(self):
        """
        Places all building block molecules on correct coordinates.

        The building block molecules are placed in their appropriate 
        positions based on the topology. This means that the 
        building-blocks* are placed on vertices and linkers on edges.
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
        
        self.macro_mol.heavy_mol = chem.Mol()
        
        lk = next(x for x in self.macro_mol.building_blocks 
                        if isinstance(x, StructUnit2))
                            
        bb = next(x for x in self.macro_mol.building_blocks 
                        if isinstance(x, StructUnit3))
                            
        for edge in self.edges:
            self.macro_mol.heavy_mol = chem.CombineMols(
                                        self.macro_mol.heavy_mol, 
                                        edge.place_mol(lk))
                                        
            heavy_ids = deque(maxlen=2)
            for atom in self.macro_mol.heavy_mol.GetAtoms():
                if atom.GetAtomicNum() in FGInfo.heavy_atomic_nums:
                    heavy_ids.append(atom.GetIdx())
            
            edge.heavy_ids = list(heavy_ids)

        for vertex in self.vertices:
            self.macro_mol.heavy_mol = chem.CombineMols(
                                        self.macro_mol.heavy_mol, 
                                        vertex.place_mol(bb))
            heavy_ids = deque(maxlen=3)
            for atom in self.macro_mol.heavy_mol.GetAtoms():
                if atom.GetAtomicNum() in FGInfo.heavy_atomic_nums:
                    heavy_ids.append(atom.GetIdx())
            
            vertex.heavy_ids = list(heavy_ids)

class CageTopology(Topology):
    """
    A topology class which cage topologies should inherit.
        
    Attributes
    ----------
    In addition to all the attributes defined within ``Topology`` this
    class has the following attributes:

    pair_up : function object (default = pair_up_edges_with_vertices)
        This is the function which pairs up molecules placed using the
        ``Vertex`` and ``Edge`` classes. This should be how cage
        topologies should be defined.    
    
    """
    
    def __init__(self, macro_mol):
        Topology.__init__(self, macro_mol)        
        self.pair_up = self.pair_up_edges_with_vertices
        

    @LazyAttr
    def windows(self):
        """
        
        Returns
        -------
        None : NoneType
            If the function for finding windows and their sizes
            found fewer than the required number of windows or
            if it failed for some other reason.
            
        list of floats
            Each float in the list represents the size of a window in
            the cage. If the window finding function found more than
            the expected number of windows, only the largest n windows
            are returned. Where n is the number of expected windows.
        
        """
        
        all_windows = window_sizes(self.macro_mol.prist_mol_file)

        if all_windows is None:          
            return None
        
        return sorted(all_windows, reverse=True)[:self.n_windows]
            

   
    def cavity_size(self):
        """
        Returns the size of the cage cavity.

        Returns
        -------
        float
            The size of the cage cavity.        
        
        """
        
        center_of_mass = self.macro_mol.center_of_mass('prist')
        min_dist = min((euclidean(coord, center_of_mass) -
        atom_vdw_radii[self.macro_mol.atom_symbol('prist', atom_id)]) 
                                 for atom_id, coord in 
                                 self.macro_mol.all_atom_coords('prist'))
        return 2 * abs(min_dist)    

    def window_difference(self, default=500):
        """
        The total difference between all cage size.
        
        Every combination of windows is considered and all the size
        differences are summed and returned.

        Parameters
        ----------
        default : float or int (default = 500)
            The number returned if no windows were found in the cage.
        
        Returns
        -------
        default : float or int
            If the `windows` attribute is ``None``. This happens when
            the window finidning algorithm fails. In these cases the
            `default` value is returned.
            
        float
            The total difference of window size when considering every
            combination of windows.
                           
        """
        
        if self.windows is None:
            return default
        
    
        # Cluster the windows into groups so that only size differences
        # between windows of the same type are taken into account. To do
        # this, first sort the windows by size. If two windows types are
        # present split the windows at the two groups at the point where
        # the window sizes have the biggest difference. If there are
        # three types split it at the two biggest differences and so on.
                
        windows = np.array(self.windows)
        
        diffs = list(abs(np.ediff1d(windows)))
        sorted_diffs = sorted(diffs, reverse=True)
        
        # Get indices of where the list should be split.
        split = []
        for x in range(self.n_window_types-1):
            i = diffs.index(sorted_diffs[x]) + 1
            split.append(i)
    
        # Get the sub-lists.
        og = list(windows)
        clusters = []
        for i in sorted(split, reverse=True):
            clusters.append(og[i:])            
            og = og[:i]

        if self.n_window_types == 1:
            clusters.append(og)

        # After this sum the differences in each group and then sum the
        # group totals.
        diff_sums = []
        for cluster in clusters:
            diff_sum = sum(abs(w1 - w2) for w1, w2 in 
                                    itertools.combinations(cluster, 2))

            diff_num = sum(1 for _ in 
                itertools.combinations(cluster, 2))
            
            diff_sums.append(diff_sum / diff_num)
            
        return sum(diff_sums)

class TwoPlusThree(CageTopology):
    vertices = [Vertex(0,0,50), Vertex(0,0,-50)]
    edges = [Edge(vertices[0], vertices[1]),
             Edge(vertices[0], vertices[1]),
             Edge(vertices[0], vertices[1])]
             
    edges[0].coord = np.array([-50,
                              -10*np.sqrt(3),
                                0])

    edges[1].coord = np.array([50,
                              -10*np.sqrt(3),
                                0])

    edges[2].coord = np.array([0,
                              10*np.sqrt(3),
                                0])

class FourPlusSix(CageTopology):
    """
    Defines the tetrahedral, 4+6, topology.

    This is a topology of cages where 4 building-blocks* are placed on
    vertices and 6 linkers are placed on the edges between them. This
    class defines functions which place these molecules in the correct
    positions within an rdkit instance. The rdkit instance is stored in 
    the `heavy_mol` attribute of a ``Cage`` instance.
        
    """
    
    
    # Vertices of a tetrahdron so that origin is at the origin. Source:
    # http://tinyurl.com/lc262h8.
    vertices = [Vertex(100,0,-100/np.sqrt(2)), 
                Vertex(-100,0,-100/np.sqrt(2)), 
                Vertex(0,100,100/np.sqrt(2)), 
                Vertex(0,-100,100/np.sqrt(2))]
        
    edges = [Edge(v1,v2) for v1, v2 in 
                itertools.combinations(vertices, 2)]  
    
    n_windows = 4
    n_window_types = 1    
    
class EightPlusTwelve(CageTopology):
    """
    Defines a square-like topology.    
    
    """
    
    vertices = [Vertex(-50, 50, -50), 
                Vertex(-50, -50, -50), 
                Vertex(50, 50, -50), 
                Vertex(50, -50, -50),

                Vertex(-50, 50, 50), 
                Vertex(-50, -50, 50), 
                Vertex(50, 50, 50), 
                Vertex(50, -50, 50)]
        
    edges = [Edge(vertices[0], vertices[2]), 
             Edge(vertices[0], vertices[1]),
             Edge(vertices[1], vertices[3]),
             Edge(vertices[2], vertices[3]),
             
             Edge(vertices[4], vertices[6]), 
             Edge(vertices[4], vertices[5]),
             Edge(vertices[5], vertices[7]),
             Edge(vertices[6], vertices[7]),


             Edge(vertices[0], vertices[4]), 
             Edge(vertices[1], vertices[5]),
             Edge(vertices[2], vertices[6]),
             Edge(vertices[3], vertices[7])]  
    
    n_windows = 6  
    n_window_types = 1
    
class SixPlusNine(CageTopology):
    
    # source: http://eusebeia.dyndns.org/4d/prism3
    vertices = [Vertex(-50,-50/np.sqrt(3),-50), 
                Vertex(-50,-50/np.sqrt(3),50),
                Vertex(50,-50/np.sqrt(3),-50),
                Vertex(50,-50/np.sqrt(3),50),
                Vertex(0, 100/np.sqrt(3),-50),
                Vertex(0, 100/np.sqrt(3), 50)]
                
    edges = [Edge(vertices[0], vertices[1]),
             Edge(vertices[0], vertices[2]),
             Edge(vertices[2], vertices[3]),
             Edge(vertices[1], vertices[3]),
             Edge(vertices[0], vertices[4]),
             Edge(vertices[2], vertices[4]),
             Edge(vertices[1], vertices[5]),
             Edge(vertices[3], vertices[5]),
             Edge(vertices[4], vertices[5])]
    
    n_windows = 5
    n_window_types = 1

class Dodecahedron(CageTopology):
    
    # Source: http://tinyurl.com/h2dl949
    phi = (1 + np.sqrt(5))/2
    x = 50
    vertices = [Vertex(x* phi, 0.0, x/phi), 
                Vertex(x*-phi, 0.0, x/phi), 
                Vertex(x*-phi, 0.0, x/-phi),
                Vertex(x* phi, 0.0, x/-phi), 
                
                Vertex(x/ phi, x* phi, 0.0), 
                Vertex(x/ phi, x*-phi, 0.0),
                Vertex(x/-phi, x*-phi, 0.0), 
                Vertex(x/-phi, x* phi, 0.0), 
                Vertex(0.0, x/ phi, x* phi),
                Vertex(0.0, x/ phi, x*-phi), 
                Vertex(0.0, x/-phi, x*-phi), 
                Vertex(0.0, x/-phi, x* phi),
  
                Vertex( x, x, x), 
                Vertex( x,-x, x), 
                Vertex(-x,-x, x), 
                Vertex(-x, x, x), 
                Vertex(-x, x,-x), 
                Vertex( x, x,-x),
                Vertex( x,-x,-x), 
                Vertex(-x,-x,-x)]

    A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T = vertices
    edges = [Edge(A,N), Edge(A,M), Edge(A,D), Edge(B,O), Edge(B,P),
             Edge(B,C), Edge(C,T), Edge(C,Q), Edge(D,S), Edge(D,R),
             Edge(E,M), Edge(E,H), Edge(E,R), Edge(F,G), Edge(F,S),
             Edge(F,N), Edge(G,O), Edge(G,T), Edge(H,P), Edge(H,Q),
             Edge(I,L), Edge(I,M), Edge(I,P), Edge(J,K), Edge(J,R),
             Edge(J,Q), Edge(K,S), Edge(K,T), Edge(L,O), Edge(L,N)]    
    
    n_windows = 12
    n_window_types = 1


    
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
        



      