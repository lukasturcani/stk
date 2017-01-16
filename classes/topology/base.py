import numpy as np
import itertools
from collections import Counter
import rdkit
import rdkit.Chem as chem
from scipy.spatial.distance import euclidean

from ..fg_info import double_bond_combs
from ...convenience_tools import (centroid, vector_theta,
                                      rotation_matrix_arbitrary_axis,
                                      normalize_vector)       

class Vertex:
    """
    Used to represent the vertices of Cage polyhedra.

    This class stores information about the vertices which make up a 
    Cage's structure.
    
    Attributes
    ----------     
    coord : numpy.array of floats
        A numpy array which holds the x, y and z coordinates of the
        vertex, in that order.
    
    connected : list of Edge or Vertex instances
        This list holds the Edge or Vertex instances which represent the 
        edges or vertices connected to the `self` vertex.
    
    bonder_ids : list of ints
        This list holds the ids of atoms which belong to the building 
        block placed on the vertex and will form bonds.
    
    atom_position_pairs : list of tuples of form (int, Vertex)
        Each atom which  is paired to a neighboring vertex or 
        edge first. Only after this is the atom - atom pairing 
        performed. The atom - edge/vertex pairing is stored here. The
        int represents the id of the atom.
        
    distances : list of tuples of form (float, int, int)
        After bonder atoms have been associated with vertices to which
        they join, the idividual atoms are paired up. To do this the 
        distance between every bonder atom on the paired vertex and the
        bonder atom which is paired to the vertex is found. This
        information is stored here where float is the distance, the
        first int is the bonder atom and the second int is the bonder
        atom on the vertex paired to the first atom.
    
    """ 
    
    def __init__(self, x, y, z):
        self.coord = np.array([x,y,z])
        self.connected = []
        self.bonder_ids = []
        self.atom_position_pairs = []
        self.distances = []
        
    @classmethod
    def vertex_init(cls, *vertices):
        obj = cls(*centroid(*(v.coord for v in vertices)))
        obj.connected.extend(vertices)
        for v in vertices:
            v.connected.append(obj)
        return obj
        
    def place_mol(self, building_block):
        """
        Places a building-block* on the coords of the vertex.
        
        The orientation of the building-block* is aligned with 2
        parameters. Firstly the normal of the plane of bonder atoms of
        the building-block* is aligned with the normal of the plane
        formed by the edges connected to the vertex. Because the normal
        of the plane of bonder atoms always points in the direction of
        the building_block*'s centroid, this alignment causes the bulk
        of the building-block* molecule to point away from the center
        of the cage.
        
        Secondly, the building-block* is rotated so that a bonder atom 
        is aligned perfectly with an edge. This reduces the rms distance
        between the edges and bonder atoms to some extent.                
        
        Parameters
        ----------
        building_block : StructUnit3
            The building-block* molecule to be placed on a vertex.
        
        Modifies
        --------
        building_block.mol : rdkit.Chem.rdchem.Mol
            The conformer of the rdkit instance in this attribute is 
            modified as per the description in the docstring.
            
        Returns
        -------
        rdkit.Chem.rdchem.Mol
            The rdkit instance holding the building-block* molecule with
            the coordinates placed on the vertex and orientation set as
            described in the docstring.
        
        """
        # Flush the list of data from previous molecules.
        self.distances = []
      
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
        theta = vector_theta(normal, self.connected[0].coord) 
        
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

        for edge1, edge2 in itertools.combinations(self.connected, 2):
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
        for edge in self.connected:
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
        return sum(edge.coord for edge in self.connected) / len(self.connected)
        

class Edge(Vertex):
    """
    Used to represent the edges of Cage polyhedra.

    This class stores information about the edges which make up a Cage's 
    structure.
    
    Attributes
    ----------

    direction : numpy.array
        This vector represents the orientation of the edge. It is a 
        normalized direction vector which runs from `v2` to `v1`.

        
        
    """  
    
    def __init__(self, v1, v2):
        Vertex.__init__(self, *centroid(v1.coord, v2.coord))
        self.direction = normalize_vector(v1.coord - v2.coord)
        self.connected.extend([v1, v2])
        v1.connected.append(self)
        v2.connected.append(self)
        
        
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
        
        # Flush the lists from data of previous molecules.
        self.distances = []
        
        # First the centroid of the heavy atoms is placed on the
        # position of the edge, then the direction of the linker is 
        # aligned with the direction of the edge.
        linker.set_heavy_atom_centroid(self.coord)
        
        flip = np.random.choice([1,-1])                
        linker.set_heavy_mol_orientation(np.multiply(self.direction,
                                                     flip))

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
    
    bb_counter : Counter
        A counter keeping track of how much of each building block was
        used during assembly. The ``StructUnit`` instance acts as the
        counter key.
        
        
    """
    
    def __init__(self, macro_mol):
        self.macro_mol = macro_mol
        self.bonds_made = 0
        self.bb_counter = Counter()         
        
    def build(self):
        """
        Assembles rdkit instances of the macromolecules.
        
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
        self.macro_mol.mol
            Adds an rdkit instance of the assembled molecule into this 
            attribute.
            
        self.bonds_made : int
            This counter is updated with each bond made during assembly.            
            
        self.bb_counter : Counter
            The counter is updated with the number of building blocks of
            each StructUnit that make up the MacroMolecule.               
            
        Returns
        -------
        None : NoneType
        
        """

        self.place_mols()
        self.join_mols()

    def join_mols(self):
        """
        Joins up the separate building blocks which form the molecule.

        Modifies
        --------
        self.macro_mol : rdkit.Chem.rdchem.Mol
            Joins up the separate building blocks in this macromolecule.

        Returns
        -------
        None : NoneType    
        
        """
        
        editable_mol = chem.EditableMol(self.macro_mol.mol)
        
        # This loop finds all the distances between an atom paired with
        # a postion and all other atoms at the paired position.
        for position in self.positions_A:
            for atom_id, vertex in position.atom_position_pairs:
                # Get all the distances between the atom and the other
                # bonding atoms on the vertex. Store this information on 
                # the vertex.
                for atom2_id in vertex.bonder_ids:
                    distance = self.macro_mol.atom_distance(atom_id, 
                                                            atom2_id)
                    position.distances.append((distance, 
                                             atom_id, atom2_id))

        # This loop creates bonds between atoms at two different
        # positions so that each atom only bonds once and so that the
        # total length of all bonds made is minimzed.   
        paired = set()        
        for position in self.positions_A:
            for _, atom1_id, atom2_id in sorted(position.distances):
                if atom1_id in paired or atom2_id in paired:
                    continue            
                
                bond_type = self.determine_bond_type(atom1_id, atom2_id)
                # Add the bond.                
                editable_mol.AddBond(atom1_id, atom2_id, bond_type)
                self.bonds_made += 1
                paired.add(atom1_id)
                paired.add(atom2_id)
                
        self.macro_mol.mol = editable_mol.GetMol()    

    def pair_bonders_with_positions(self, vertex):
        """
        Matches atoms with the closest building block position.
        
        After a building block is placed on a position, each atom
        which forms a bond must be paired with the location of the 
        building block to which it bonds. This function matches atoms
        and positions so that each is only present in one pairing and so
        that the total distance of the pairings is minimized.

        Parameters
        ----------
        vertex : Vertex
            The position at which all the atoms being paired are 
            located.

        Modifies
        --------
        vertex.atom_position_pairs : list of tuples of (int, Vertex)
            Adds a tuples to this list represnting the id of the atom
            and position which were paired.
        
        Returns
        -------
        None : NoneType
        
        """

        # This loop looks at each atom which forms a new bond and all
        # the positions (not atoms) to which it may end up bonding. It
        # finds the distances of all the options.
        distances = []
        for bonder_id in vertex.bonder_ids:
            atom_coord = self.macro_mol.atom_coords(bonder_id)
            for position in vertex.connected:                
                distance = euclidean(atom_coord, position.coord)
                distances.append((distance, bonder_id, position))
 
        # Sort the pairings of atoms with potential bonding position,
        # smallest first. 
        distances.sort()
        
        # This loop looks at all the potential pairings of atoms to
        # positions. It pairs the shortest combinations of atoms and
        # positions, making sure that each atom and position is only
        # paired once. The pairings are saved to the `
        # atom_positions_pairs` attribute of the position on which all
        # the bonder atoms are placed.
        paired_pos = set()
        paired_ids = set()
        vertex.atom_position_pairs = []
        for _, bonder_id, pos in distances:
            if bonder_id in paired_ids or pos in paired_pos:
                continue
            vertex.atom_position_pairs.append((bonder_id, pos))
            paired_ids.add(bonder_id)
            paired_pos.add(pos)

    def determine_bond_type(self, atom1_id, atom2_id):
        """
        Returns the bond order to be formed between the atoms.
        
        Some atoms will need to have a double bond created between them.
        This is defined in the `double_bond_combs` list. If the
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
            in `double_bond_combs`.
        
        rdkit.Chem.rdchem.BondType.DOUBLE
            If the combination of heavy atoms passed as arguments is in
            `double_bond_combs`.
            
        """
        
        # Get the functional groups of the of the atoms whose atom ids 
        # were supplied as arguments. If the groups form a tuple in 
        # `double_bond_combs` return a rdkit double bond type. If they 
        # do not, return a rdkit single bond type.
        
        atom1 = self.macro_mol.mol.GetAtomWithIdx(atom1_id)
        atom1_grp = atom1.GetProp('fg')
        atom2 = self.macro_mol.mol.GetAtomWithIdx(atom2_id)
        atom2_grp= atom2.GetProp('fg')      
        
        double_bond_present = ((atom1_grp, atom2_grp) == tup or 
                               (atom2_grp, atom1_grp) == tup for 
                               tup in double_bond_combs)
        
        if any(double_bond_present):
            return rdkit.Chem.rdchem.BondType.DOUBLE
        else:
            return rdkit.Chem.rdchem.BondType.SINGLE