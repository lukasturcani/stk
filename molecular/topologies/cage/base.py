import itertools
from collections import deque
from scipy.spatial.distance import euclidean
import numpy as np
import rdkit.Chem.AllChem as rdkit

from ..base import Topology
from ....convenience_tools import (centroid, vector_theta,
                                      rotation_matrix_arbitrary_axis,
                                      normalize_vector)

class WindowError(Exception):
    def __init__(self, message):
        self.message = message

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
        This list holds the Edge or Vertex instances which represent
        the edges or vertices connected to the `self` vertex.

    bonder_ids : list of ints
        This list holds the ids of atoms which belong to the building
        block placed on the vertex and will form bonds. The ids reflect
        those in the macromolecule not in the building block itself.

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
        """
        Intializes the Vertex from  a list of other vertices.

        This initalizer automatically calculates the position of the
        vertex from the positions of the vertices provided to the
        initializer. Its position is set to the centroid of the
        provided vertices.

        The `connected` attributes of all the involved vertices are
        also updated.

        Parameters
        ----------
        vertices : tuple
            A tuple of Vertex objects.

        Returns
        -------
        Vertex
            The initialized vertex.

        """

        # Get the position of this Vertex, as the centroid of the
        # supplied vertices.
        obj = cls(*centroid(*(v.coord for v in vertices)))
        # Update the `connected` attributes.
        obj.connected.extend(vertices)
        for v in vertices:
            v.connected.append(obj)
        return obj

    def place_mol(self, building_block, aligner=0):
        """
        Place a StructUnit3 building block on the coords of the vertex.

        The orientation of the building block is aligned with 2
        parameters. Firstly, the normal of the plane of bonder atoms of
        the building block is aligned with the normal of the plane
        formed by the edges connected to the vertex. Because the normal
        of the plane of bonder atoms always points in the direction of
        the building block's centroid, this alignment causes the bulk
        of the building block  molecule to point away from the center
        of the cage.

        Secondly, the building block is rotated so that a bonder atom
        is aligned perfectly with an edge. This reduces the rms
        distance between the edges and bonder atoms to some extent,
        probably.

        Parameters
        ----------
        building_block : StructUnit3
            The building block molecule to be placed on a vertex.

        aligner : int (default = 0)
            The index of the atom within `bonder_ids` which is to be
            aligned with an edge.

        Modifies
        --------
        building_block.mol : rdkit.Chem.rdchem.Mol
            The conformer of the rdkit instance in this attribute is
            modified as per the description in the docstring.

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            The rdkit instance holding the building block molecule with
            the coordinates placed on the vertex and orientation set as
            described in the docstring.

        """
        # Flush the list of data from previous molecules.
        self.distances = []

        # The method first aligns the normal of the bonder atom plane
        # to the normal of the edge plane. This means the bulk of the
        # building block is always pointed away from the center of the
        # molecule.
        building_block.set_orientation2(self.edge_plane_normal())

        # Next, the building block must be rotated so that one of the
        # bonder atoms is perfectly aligned with one of the edges. This
        # is a multi-step process:
        #   1) Place the centroid of the bonder atoms at the origin.
        #   2) Place the centroid of the edges connected to the vertex
        #      at the origin.
        #   3) Rotate the building block by some amount `theta`, so
        #      so that one of the bonder atoms is perfectly aligned
        #      with
        #      one of the edges. The axis of rotation is the normal to
        #      the plane of bonder atoms.
        #
        # The rotation is carried out via matrices. This means a
        # coordinate matrix of atoms in the molecule is generated
        # and modified.

        # Set the centroid of the bonder atoms to the origin.
        building_block.set_bonder_centroid([0,0,0])
        # Get the coordinate of the atom which is to be aligned with an
        # edge.
        atom_coord = building_block.atom_coords(
                                 building_block.bonder_ids[aligner])

        # Get the coordinates of all the edges and translate the
        # centroid to the origin.
        edge_coord_mat = (self.edge_coord_matrix() -
                                                self.edge_centroid())
        edge_coord = np.array(edge_coord_mat[0,:])[0]

        # Get the angle between an edge and the atom.
        theta = vector_theta(edge_coord, atom_coord)

        # Get the rotation matrix necessary to do the rotation of
        # `theta` about the normal to the plane.
        rot_mat = rotation_matrix_arbitrary_axis(theta,
                                           self.edge_plane_normal())
        # Apply the rotation to the positions of the atoms and get a
        # position matrix of the new coordinates.
        pos_mat = building_block.position_matrix()
        new_pos_mat = np.dot(rot_mat, pos_mat)
        # Update the atomic positions in the building block.
        building_block.set_position_from_matrix(new_pos_mat)

        # Finally the well orientated building-block* is placed on the
        # coords of the vertex.
        building_block.set_bonder_centroid(self.coord)
        return building_block.mol

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
        coordinate of some point on the plane. For example, the
        position of one of the connected edges.

        Returns
        -------
        numpy.array
            This array has the form [a, b, c, d] and represents the
            scalar equation of the plane formed by the connected edges.

        References
        ----------
        https://tinyurl.com/okpqv6

        """

        edge_coord = self.edges[0].coord
        d = np.multiply(np.sum(np.multiply(self.edge_plane_normal(),
                                           edge_coord)), -1)
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
        return (sum(edge.coord for edge in self.connected) /
                len(self.connected))

class Edge(Vertex):
    """
    Used to represent the edges of Cage polyhedra.

    This class stores information about the edges which make up a
    Cage's structure.

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


    def place_mol(self, linker, flip=True):
        """
        Places a linker molecule on the coordinates of an edge.

        It also orientates the linker so that a the bonder atoms sit
        exactly on the edge and the bulk of the linker points away from
        the center of the cage.

        Parameters
        ----------
        linker : StructUnit2
            The linker which is to be placed and orientated as
            described in the docstring.

        flip : bool (default = True)
            If ``True`` the linker has a 50% chance of being rotated
            so that its direction vector runs in the opposite direction
            to the edge's direction vector.

        Modifies
        --------
        linker.mol : rdkit.Chem.rdchem.Mol
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

        # First the centroid of the bonder atoms is placed on the
        # position of the edge, then the direction of the linker is
        # aligned with the direction of the edge.
        linker.set_bonder_centroid(self.coord)

        # If `flip` is ``True`` the linker rotated at randomly chosen
        # to align either with the edge direction vector or its
        # opposite direction.
        f = 1
        if flip:
            f = np.random.choice([1,-1])
        linker.set_orientation2(np.multiply(self.direction, f))

        # Ensure the centroid of the linker is placed on the outside of
        # the cage.
        linker.minimize_theta(self.coord, self.direction)

        return linker.mol

class _CageTopology(Topology):
    """
    A topology class which cage topologies should inherit.

    Attributes
    ----------
    In addition to all the attributes defined within ``Topology`` this
    class has the following attributes:

    alignment : list of ints or None
        This length of this list must be equal to the number of
        vertices in the cage. When cages are built one of the bonder
        atoms is aligned with an edge during placement. The int
        indicates which bonder atom is aligned. The int corresponds to
        an index in `bonder_ids`.

        If ``None`` the first atom in `bonder_ids` is always aligned.

    """

    def __init__(self, alignment=None):
        super().__init__()
        self.alignment = alignment

    def join_mols(self, macro_mol):
        """
        Joins up the separate building blocks which form the molecule.

        Parameters
        ----------
        macro_mol : MacroMolecule
            The macromolecule being assembled.

        Modifies
        --------
        macro_mol.mol : rdkit.Chem.rdchem.Mol
            Joins up the separate building blocks in this
            macromolecule.

        macro_mol.bonds_made : int
            Places the number of bonds made during assembly into this
            attribute.

        Returns
        -------
        None : NoneType

        """

        editable_mol = rdkit.EditableMol(macro_mol.mol)
        macro_mol.bonds_made = 0

        # This loop finds all the distances between an atom paired with
        # a postion and all other atoms at the paired position.
        for position in self.positions_A:
            for atom_id, vertex in position.atom_position_pairs:
                # Get all the distances between the atom and the other
                # bonding atoms on the vertex. Store this information
                # on the vertex.
                for atom2_id in vertex.bonder_ids:
                    distance = macro_mol.atom_distance(atom_id,
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

                bond_type = self.determine_bond_type(macro_mol,
                                                     atom1_id,
                                                     atom2_id)
                # Add the bond.
                editable_mol.AddBond(atom1_id, atom2_id, bond_type)
                macro_mol.bonds_made += 1
                paired.add(atom1_id)
                paired.add(atom2_id)

        macro_mol.mol = editable_mol.GetMol()

    def pair_bonders_with_positions(self, macro_mol, vertex):
        """
        Matches atoms with the closest building block position.

        After a building block is placed on a position, each atom
        which forms a bond must be paired with the location of the
        building block to which it bonds. This function matches atoms
        and positions so that each is only present in one pairing and
        so that the total distance of the pairings is minimized.

        Parameters
        ----------
        macro_mol : MacroMolecule
            The macromolecule being buit.

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
            atom_coord = macro_mol.atom_coords(bonder_id)
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

    def place_mols(self, macro_mol):
        """
        Places all building block molecules on correct coordinates.

        The building block molecules are placed in their appropriate
        positions based on the topology. It does not join them.

        Parameters
        ----------
        macro_mol : MacroMolecule
            The macromolecule being built.

        Modifies
        --------
        macro_mol.mol
            An rdkit instance of the macromolecule with disconnected
            building blocks is placed in this attribute.

        macro_mol.bb_counter : Counter
            The counter is updated with the number of building blocks
            of each type used to form the macromolecule.

        Returns
        -------
        None : NoneType

        """

        macro_mol.mol = rdkit.Mol()

        # Get the StructUnit instances of the building blocks.
        bb1, bb2 = macro_mol.building_blocks
        # Get the number of functional groups in each building block.
        n_fg1 = len(bb1.functional_group_atoms())
        n_fg2 = len(bb2.functional_group_atoms())

        # Depending on the number of functional groups, assigned a
        # building block to be either a linker or a building-block*.
        if n_fg1 < n_fg2:
            lk = bb1
            n_lk = n_fg1
            bb = bb2
            n_bb = n_fg2
        else:
            lk = bb2
            n_lk = n_fg2
            bb = bb1
            n_bb = n_fg1

        # Save the original orientations of the linker and building
        # block. This means that when orienation of molecules is done,
        # the starting position is always the same. Ensures
        # consistency.
        lk_pos = lk.position_matrix()
        bb_pos = bb.position_matrix()

        # This loop places all building-blocks* on the points at
        # `positions_A`. It then pairs all atoms which form a new bond
        # with the positions to which they will be bonding. It also
        # counts the nubmer of building-blocks* which make up the
        # structure.
        for i, position in enumerate(self.positions_A):

            # Get the id of atom which to be aligned with an edge.
            aligner = (0 if self.alignment is None else
                       self.alignment[i])
            # Position the molecule on the vertex.
            bb.set_position_from_matrix(bb_pos)
            bb_mol = position.place_mol(bb, aligner)
            macro_mol.mol = rdkit.CombineMols(macro_mol.mol, bb_mol)
            # Update the counter each time a building-block* is added.
            macro_mol.bb_counter.update([bb])

            # Get ids of atoms which form new bonds.
            bonder_ids = deque(maxlen=n_bb)
            for atom in macro_mol.mol.GetAtoms():
                if atom.HasProp('bonder'):
                    bonder_ids.append(atom.GetIdx())

            # Save the ids of atoms which form new bonds and pair them
            # up with positions.
            position.bonder_ids = sorted(bonder_ids)
            self.pair_bonders_with_positions(macro_mol, position)

        # This loop places all linkers on the points at `positions_B`.
        # It then saves all atoms which form a new bond to the position
        # they are found at. It also counts the number of linkers which
        # make up the structure.
        for position in self.positions_B:

            # Here `place_mol` can either belong to an Edge object
            # or a Vertex object. `self.alignment` can be either
            # ``None`` or a list of ints. If it is ``None`` it means
            # that a linker placed on the Edge should be flipped at
            # random. This means the second argument should be 1 to
            # turn on the flipping. If its ``None`` and `place_mol`
            # belongs to a Vertex it just means the second bonder atom
            # gets aligned. This has no practical effect. By default
            # the first bonder atom gets aligned, but this has no
            # objective benefit.

            # If `self.alignment` is a list of ints it means that
            # flipping should be turned off for edges and that the 0th
            # atom will get aligned if the position is a Vertex.
            lk.set_position_from_matrix(lk_pos)
            lk_mol = position.place_mol(lk, int(not self.alignment))
            macro_mol.mol = rdkit.CombineMols(macro_mol.mol, lk_mol)
            # Update the counter each time a linker is added.
            macro_mol.bb_counter.update([lk])

            # Get ids of atoms which form new bonds.
            bonder_ids = deque(maxlen=n_lk)
            for atom in macro_mol.mol.GetAtoms():
                if atom.HasProp('bonder'):
                    bonder_ids.append(atom.GetIdx())

            # Save the ids of atoms which form new bonds.
            position.bonder_ids = list(bonder_ids)

class _VertexOnlyCageTopology(_CageTopology):

    def __init__(self, random_placement=True):
        super().__init__()
        self.random_placement = random_placement
        self.connect()

    def place_mols(self, macro_mol):

        macro_mol.mol = rdkit.Mol()

        if self.random_placement:
            return self.place_mols_random(macro_mol)
        return self.place_mols_assigned(macro_mol)

    def place_mols_random(self, macro_mol):
        for position in self.positions_A:
            bb = np.random.choice(list(macro_mol.building_blocks))
            n_bb = len(bb.functional_group_atoms())

            macro_mol.mol = rdkit.CombineMols(
                                            macro_mol.mol,
                                            position.place_mol(bb))
            macro_mol.bb_counter.update([bb])

            bonder_ids = deque(maxlen=n_bb)
            for atom in macro_mol.mol.GetAtoms():
                if atom.HasProp('bonder'):
                    bonder_ids.append(atom.GetIdx())

            position.bonder_ids = sorted(bonder_ids)
            self.pair_bonders_with_positions(macro_mol, position)

    @classmethod
    def connect(cls):
        if getattr(cls, 'connected', False):
            return

        for v1, v2 in cls.connections:
            v1.connected.append(v2)
            v2.connected.append(v1)
        cls.connected = True

    def place_mols_assigned(self):
        pass
