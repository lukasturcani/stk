import itertools
from collections import deque
from scipy.spatial.distance import euclidean
import numpy as np
import rdkit.Chem.AllChem as rdkit

from ..base import Topology
from ....convenience_tools import (centroid, vector_theta,
                                   add_fragment_props,
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

    id : object, optional
        An id to identify the vertex. Used by `place_mols()`

    """

    def __init__(self, x, y, z, id_=None):
        self.coord = np.array([x,y,z])
        self.connected = []
        self.bonder_ids = []
        self.atom_position_pairs = []
        self.distances = []
        self.id = id_

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

    def place_mol(self, building_block, aligner=0, aligner_edge=0):
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

        aligner : int, optional
            The index of the atom within `bonder_ids` which is to be
            aligned with an edge.

        aligner_edge : int, optional
            The index of an edge in `connected`. It is the edge with
            which `aligner` is aligned.

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

        # Next, define the direction vector going from the edge
        # centroid to the edge with which the atom is aligned.
        building_block.set_bonder_centroid(self.coord)
        vector = (self.connected[aligner_edge].coord -
                  self.edge_centroid())
        # Get the id of the atom which is being aligned.
        atom = building_block.bonder_ids[aligner]
        # Minimize the angle between these things by rotating about the
        # normal of the edge plane.
        building_block.minimize_theta(atom,
                                      vector,
                                      self.edge_plane_normal())

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
            normal *= -1

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
        d = -np.sum(self.edge_plane_normal() * edge_coord)
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

    def __init__(self, v1, v2, id_=None):
        Vertex.__init__(self, *centroid(v1.coord, v2.coord), id_)
        self.direction = normalize_vector(v1.coord - v2.coord)
        self.connected.extend([v1, v2])
        v1.connected.append(self)
        v2.connected.append(self)

    def place_mol(self, linker, alignment):
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

        alignment : int
            1 for parallel alignment with `self.direction` and -1 for
            anti-parallel alignment with `self.direction`.


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
        linker.set_orientation2(self.direction * alignment)

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

    A_alignments : list of ints (default = None)
        The length of this list must be equal to the number of
        building blocks in the cage. When cages are built one of the
        bonder atoms of each building block is aligned with an edge
        during placement. The int indicates which bonder atom is
        aligned. The int corresponds to an index in `bonder_ids`.

        If ``None`` the first atom in `bonder_ids` is always aligned.

        For example,

            A_alignments = [0,2,1,2]

        In this case there must be 4 building blocks in the cage. The
        the first building block has its first (index 0) bonder atom
        aligned. The 2nd building block has the 3rd (index 2) atom
        aligned. The 3rd building block has the 2nd (index 1) atom
        aligned. The 4th building block has the 3rd (index 2) atom
        aligned.


    B_alignments : list of ints (default = None)
        The length of this list should be euqal to the number of
        linkers in the cage. The linkers of a cage can have either 2
        functional groups or 3 or more, depending on the topology.

        When the linkers have 2 functional groups the list should hold
        either 1 or -1. The value indicates that the linker is aligned
        parallel or antiparallel with the edge its placed on. For
        example in a tetrahedral topology,

            B_alignments = [-1, 1, 1, -1, 1, -1]

        It indicates that the first linker is aligned anti parallel and
        the second is aligned in parallel, and so on.

        If the linkers have 3 or more functional groups, the values
        in B_alignments have the same role as `A_alignments`. The only
        difference is that by default the second atom is aligned,
        rather than the first.

    edge_alignments : list of ints, optional
        The length of the list is equal to the number of building
        blocks in the cage. Each element is an int which holds the id
        of an edge. For example,

            edge_alignments = [1, 2, 3, 4]

        then the first building block is aligned with the edge with
        `id` of 1, the second building block is aligned with the edge
        with `id` 2 and so on. For this to work the edge defined by
        the class must have their `id` attributes defined.

    """

    def __init__(self, A_alignments=None, B_alignments=None,
                 edge_alignments=None):
        if A_alignments is None:
            A_alignments = np.zeros(len(self.positions_A))
        if B_alignments is None:
            B_alignments = np.ones(len(self.positions_B))
        if edge_alignments is None:
            edge_alignments = [None for x in
                                        range(len(self.positions_A))]

        super().__init__()
        self.A_alignments = A_alignments
        self.B_alignments = B_alignments
        self.edge_alignments = edge_alignments

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
            # Position the molecule on the vertex.
            bb.set_position_from_matrix(bb_pos)
            aligner_edge_id = self.edge_alignments[i]
            aligner_edge = next((position.connected.index(x) for x in
                                 position.connected if
                                 x.id == aligner_edge_id), 0)
            bb_mol = position.place_mol(bb, int(self.A_alignments[i]),
                                        aligner_edge)
            add_fragment_props(bb_mol,
                               macro_mol.building_blocks.index(bb),
                               i)

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
        for i, position in enumerate(self.positions_B):
            lk.set_position_from_matrix(lk_pos)
            lk_mol = position.place_mol(lk,  int(self.B_alignments[i]))
            add_fragment_props(lk_mol,
                               macro_mol.building_blocks.index(lk),
                               i)

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
    """
    Cage topolgies where all building block/linker mols have 3+ fgs.

    """

    def __init__(self, A_alignments=None, B_alignments=None,
                        edge_alignments=None):

        if A_alignments is None:
            A_alignments = np.zeros(len(self.positions_A))
        if B_alignments is None:
            B_alignments = np.zeros(len(self.positions_B))
        if edge_alignments is None:
            edge_alignments = [None for x in
                                        range(len(self.positions_A))]

        self.A_alignments = A_alignments
        self.B_alignments = B_alignments
        self.edge_alignments = edge_alignments

class _NoLinkerCageTopology(_CageTopology):
    """
    Cage topologies where all building units have equal number of fgs.

    Attributes
    ----------
    alignments : list of ints
        Same meaning as `A_alignments` in _CageTopology.

    placement : str
        The name of the placement type to be used. Valid options are

            'random' - For each vertex, a building block is picked
             randomly from `macro_mol.building_blocks` and placed on
             the vertex.

    """

    def __init__(self, alignments=None, placement='random'):
        if alignments is None:
            alignments = np.zeros(len(self.positions_A))

        self.alignments = alignments
        self.placement = placement
        self.connect()

    def place_mols(self, macro_mol):

        macro_mol.mol = rdkit.Mol()

        if self.placement == 'random':
            return self.place_mols_random(macro_mol)

    def place_mols_random(self, macro_mol):
        for position, orientation in zip(self.positions_A,
                                         self.alignments):
            bb = np.random.choice(list(macro_mol.building_blocks))
            ipos = bb.position_matrix()
            n_bb = len(bb.functional_group_atoms())

            macro_mol.mol = rdkit.CombineMols(macro_mol.mol,
                             position.place_mol(bb, int(orientation)))
            macro_mol.bb_counter.update([bb])

            bonder_ids = deque(maxlen=n_bb)
            for atom in macro_mol.mol.GetAtoms():
                if atom.HasProp('bonder'):
                    bonder_ids.append(atom.GetIdx())

            position.bonder_ids = sorted(bonder_ids)
            self.pair_bonders_with_positions(macro_mol, position)
            bb.set_position_from_matrix(ipos)

    @classmethod
    def connect(cls):
        """
        Updates each Vertex with a list of its neighbors.

        Returns
        -------
        None : NoneType

        """

        if getattr(cls, 'connected', False):
            return

        for v1, v2 in cls.connections:
            v1.connected.append(v2)
            v2.connected.append(v1)
        cls.connected = True
