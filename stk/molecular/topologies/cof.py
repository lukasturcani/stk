"""
Defines COF topologies.

.. _`cof assembly`:

How the assembly process works for 2D COF topologies.
-----------------------------------------------------

Linker topologies.
..................

Linker topologies are those which use a building block with two
functional groups and a building block with three or more functional
groups. These are implemented as subclasses to
:class:`LinkerCOFLattice`.

As for any assembly process, there are 3 stages: placement of
building blocks, joining of building blocks and deletion of extra
atoms.

A 2D COF topology can be defined by the vertices and edges which
connected them. The building block with 3 or more functional groups
sits on the vertices while the building block with 2 functional groups
sits on the edges. To assemble a COF topology, first the vertices
and edges are defined. They are represented by the classes
:class:`Vertex` and :class:`Edge`. Each of these objects defines their
positions in term of the fractional coordinates along the
``a``, ``b`` and ``c`` vectors of a unit cell. In addtion, each
of these objects defines a set of "atomic positiions".

For edges, there are two atomic positions, 0 and 1. 0 is the position
of the atom which gets connected to an atom on the first vertex
connected by the edge, while 1 is the
position of the atom which gets connected to the second vertex
connected by the edge. Identification of which bonder atom of a
di-functinalized building block sits on which atom position is done
in the following way. Check the distance of each bonder atom to the
first vertex. The closer one is on position 0 the further one is on
position 1.

For vertices there are multiple atomic positions. The number of
atomic positions on a vertex is equal to the number of functional
groups the building block placed on the vertex has. The atomic
positions are labelled from ``0`` to ``n`` where ``n`` is one less than
the number of functional groups in the building block. Going clockwise
around the vertex, starting at the 12 o clock position, the positions
are labelled starting at 0. After placing and aligning a building block
on a vertex, it is necessary to identify which bonder atom sits on
which atomic position. Alignment is done by rotating the builing block
around the z axis so that the distance of one of the bonder atoms
and one of the edges is minimized. Because each edge keeps a record
of which position it connects, the aligned bonder can be assigned to
a position. Then the bonder atoms are assigned to the next positions
in clockwise order.

To join molecules simply go through all the edges in a topology. Each
edge holds the vertices it is connected to and the positions on those
vertices it connects. Each edge and vertex defines
a dictionary which maps the position to the bonder atom which sits on
it. For each edge take the bonder at position ``0`` and the atomic
position on the vertex which is connceted to it. Use the atomic
position on the connected vertex to get the id of the bonder atom which
is to be connected. The same can be done for the bonder at position
``1`` on the edge.

Periodic bonds are not added to the rdkit molelecule, instead they
are registered in :attr:`.Periodic.periodic_bonds`. The method
:meth:`COFLattice.del_atoms` is slightly modified, so that the
coordinates of any deleted atoms are saved. This is necessary for
restoring those atoms when using :meth:`.Periodic.island`.

When creating periodic bonds during :meth:`LinkerCOFLattice.join_mols`
the ids of the bonder atoms involved are saved as the indices within
:attr:`.MacroMolecule.bonder_ids`. This is because deleting atoms
will change the atom ids but not the ordering in
:attr:`~.MacroMolecule.bonder_ids`. However, these are automatically
updated to bonder ids after assembly is completed by
:meth:`.Perdiodic.save_ids`. This is run automatically after assembly
and does not impact implementation but is noted here for completeness.

"""

import rdkit.Chem.AllChem as rdkit
import numpy as np
from scipy.spatial.distance import euclidean
from collections import deque, defaultdict

from .base import Topology
from ...convenience_tools import (PeriodicBond,
                                  add_fragment_props,
                                  normalize_vector)


class Vertex:
    """
    Represents a vertex of a 2D COF topology.

    Attributes
    ----------
    frac_coord : :class:`tuple` of :class:`float`
        The fractional positions along the ``a``, ``b`` and ``c``
        vectors of the vector.

    connected : :class:`list` of :class:`Edge`
        The edges connected to this vertex.

    bonder_map : :class:`dict`
        Maps the id of a bonder position on the vertex to the atom
        id of the bonder sitting on it.

    """

    def __init__(self, frac_coord):
        self.frac_coord = frac_coord
        self.connected = []

    def place_mol(self, cell_params, mol, aligner):
        """
        Places and aligns a molecule on the vertex.

        Parameters
        ----------
        cell_params : :class:`list` of :class:`numpy.array`
            The ``a``, ``b`` and ``c`` vectors of the unit cell.

        mol : :class:`.StructUnit3`
            The building block to be placed.

        aligner : :class:`int`
            The index of a bonder atom in
            :attr:`.StructUnit.bonder_ids` of `mol`. This is the
            bonder atom which gets aligned with an edge.

        Returns
        -------
        :class:`rdkit.Chem.rdchem.Mol`
            The rdkit instance of the placed `mol`.

        """

        coord = self.calc_coord(cell_params)
        original_position = mol.position_matrix()
        mol.set_orientation2([0, 0, 1])

        mol.set_bonder_centroid(coord)
        aligner_edge = next((e for e in self.connected if
                             all(b == 0 for b in e.bond)),
                            self.connected[0])
        vector = (aligner_edge.calc_coord(cell_params) - coord)
        atom = mol.bonder_ids[aligner]
        mol.minimize_theta2(atom, vector, [0, 0, 1])

        rdkit_mol = rdkit.Mol(mol.mol)
        mol.set_position_from_matrix(original_position)
        return rdkit_mol

    def calc_coord(self, cell_params, cell_position=[0, 0, 0]):
        """
        Calculate the position of the vertex within a unit cell.

        Parameters
        ----------
        cell_params : :class:`list` of :class:`numpy.array`
            The ``a``, ``b`` and ``c`` vectors of the unit cell.

        cell_position : :class:`list` of :class:`int`
            Indicates if the vertex is in a periodic cell or not.
            Analogous to the defintion of a periodic bond.

        Returns
        -------
        :class:`numpy.array`
            The position.

        """

        coord = np.zeros((3, ))
        for frac, dim, p in zip(self.frac_coord, cell_params, cell_position):
            coord += (frac+p) * dim

        return coord

    def create_bonder_map(self,
                          macro_mol,
                          cell_params,
                          nbonders,
                          aligned_bonder):
        """
        Creates the attribute :attr:`bonder_map`.

        Parameters
        ----------
        macro_mol : :class:`.MacroMolecule`
            The macromolecule being built.

        cell_params : :class:`list` of :class:`numpy.array`
            The ``a``, ``b`` and ``c`` vectors of the unit cell.

        nbonders : :class:`int`
            The number of bonder atoms the building block placed on the
            vertex has.

        aligned_bonder : :class:`int`
            The index of a bonder atom in
            :attr:`.StructUnit.bonder_ids` of `mol`. This is the
            bonder atom which gets aligned with an edge.

        Returns
        -------
        None : :class:`NoneType`

        """

        center = self.calc_coord(cell_params)

        # Get all the bonder ids.
        bonder_ids = deque(maxlen=nbonders)
        for atom in macro_mol.mol.GetAtoms():
            if atom.HasProp('bonder'):
                bonder_ids.append(atom.GetIdx())

        start = np.array([0, 1])
        angles = []
        for bonder in bonder_ids:
            bonder_coords = macro_mol.atom_coords(bonder) - center
            x, y, _ = normalize_vector(bonder_coords)
            angle = np.arccos(start@np.array([x, y]))
            if x < 0:
                angle = 2*np.pi - angle

            angles.append((angle, bonder))
        angles.sort()

        bonders = [bonder for angle, bonder in angles]
        aligned = bonder_ids[aligned_bonder]
        bonders = (bonders[bonders.index(aligned):] +
                   bonders[:bonders.index(aligned)])

        aligned_pos = self.aligned_position()
        positions = (list(range(aligned_pos, nbonders)) +
                     list(range(0, aligned_pos)))

        self.bonder_map = {}
        for bonder, position in zip(bonders, positions):
            self.bonder_map[position] = bonder

    def aligned_position(self):
        """
        Returns the position of the atom aligned by :meth:`place_mol`.

        Returns
        -------
        :class:`int`
            The id of a position of on the vertex. The bonder atom
            which gets aligned during :meth:`place_mol` sits on this
            position.

        """

        aligner_edge = next((e for e in self.connected if
                             all(b == 0 for b in e.bond)),
                            self.connected[0])
        vindex = 0 if self is aligner_edge.v1 else 1
        return aligner_edge.joint_positions[vindex]


class Edge:
    """
    Represents a bond in a COF unit cell.

    This class is used for positioning and aligning ditopic building
    blocks in 2D COF topologies.

    Attributes
    ----------
    v1 : :class:`Vertex`
        The first vertex connected by the edge.

    v2 : :class:Vertex`
        The second vertex connected by the edge.

    joint_positions : :class:`tuple` of :class:`int`
        Holds 2 numbers. The first is the id of the position on
        :attr:`v1` which gets connected by the edge. The second is the
        id of the position on :attr:`v2` which gets connected by the
        edge.

    bond : :class:`list` of :class:`int`:
        Indicates if the bond is periodic.

    connected : class:`list` of :class:`Vertex`
        The vertices conncetd by the edge.

    frac_coord : :class:`list` of :class`float`
        The position of the edge in terms of the fractional positions
        along the ``a``, ``b`` and ``c`` vectors of a unit cell.

    bonder_map : :class:`dict`
        Maps the id of a position to the id of the bonder atom which
        sits on it. The positions are ``0`` and ``1`` where ``0`` is
        the  position which is closer to :attr:`v1` while ``1``
        is the position closer to :attr:`v2`.

    """

    def __init__(self, v1, v2, joint_positions, bond=[0, 0, 0]):
        self.v1 = v1
        self.v2 = v2
        self.joint_positions = joint_positions
        self.bond = bond
        self.connected = [v1, v2]

        self.frac_coord = []
        for f1, f2, b in zip(v1.frac_coord, v2.frac_coord, bond):
            self.frac_coord.append((f1+f2+b) / 2)

        v1.connected.append(self)
        v2.connected.append(self)

    def place_mol(self, macro_mol, cell_params, mol, alignment):
        """
        Places and aligned a building block along the edge.

        Parameters
        ----------
        macro_mol : :class:`.MacroMolecule`
            The macromolecule being built.

        cell_params : :class:`list` of :class:`numpy.array`
            The ``a``, ``b`` and ``c`` vectors of the unit cell.

        mol : :class:`.StructUnit3`
            The building block to be placed.

        alignment : :class:`int`
            Can be ``1`` or ``-1`` to align `mol` either parallel or
            anti-parallel with the edge.

        Returns
        -------
        :class:`rdkit.Chem.rdchem.Mol`
            The rdkit instance of the placed `mol`.

        """

        coord = self.bonder_centroid(macro_mol, cell_params)
        original_position = mol.position_matrix()

        mol.set_bonder_centroid(coord)
        d = self.bonder_direction(macro_mol, cell_params)*alignment
        mol.set_orientation2(d)

        rdkit_mol = rdkit.Mol(mol.mol)
        mol.set_position_from_matrix(original_position)
        return rdkit_mol

    def bonder_direction(self, macro_mol, cell_params):
        """
        Calculates the direction vetor between the bonder atoms.

        Parameters
        ----------
        macro_mol : :class:`.MacroMolecule`
            The macromolecule being assembled.

        cell_params : :class:`list` of :class:`numpy.array`
            The ``a``, ``b`` and ``c`` vectors of the unit cell.

        Returns
        -------
        :class:`numpy.array`
            The direction vector.

        """

        coords = []
        for i, position in enumerate(self.joint_positions):
            vertex = self.connected[i]
            bonder = vertex.bonder_map[position]
            coords.append(macro_mol.atom_coords(bonder))

        for d, param in zip(self.bond, cell_params):
                coords[1] += d*param

        return normalize_vector(coords[1] - coords[0])

    def direction(self, cell_params):
        """
        The direction vector of the edge.

        Parameters
        ----------
        cell_params : :class:`list` of :class:`numpy.array`
            The ``a``, ``b`` and ``c`` vectors of the unit cell.

        Returns
        -------
        :class:`numpy.array`
            The direction vector.

        """

        v1coord = self.v1.calc_coord(cell_params)
        v2coord = self.v2.calc_coord(cell_params, self.bond)
        return normalize_vector((v1coord-v2coord)/2)

    def calc_coord(self, cell_params):
        """
        Calcualtes the position of the edge in a cell.

        Parameters
        ----------
        cell_params : :class:`list` of :class:`numpy.array`
            The ``a``, ``b`` and ``c`` vectors of the unit cell.

        Returns
        -------
        :class:`numpy.array`
            The position.
        """

        coord = np.zeros((3, ))
        for frac, dim in zip(self.frac_coord, cell_params):
            coord += frac * dim
        return coord

    def bonder_centroid(self, macro_mol, cell_params):
        """
        The centroid of the bonder molecules connected to the edge.

        So the bonders NOT sitting on the edge itself.

        Parameters
        ----------
        macro_mol : :class:`.MacroMolecule`
            The macromolecule being assembled.

        cell_params : :class:`list` of :class:`numpy.array`
            The ``a``, ``b`` and ``c`` vectors of the unit cell.

        Returns
        -------
        :class:`numpy.array`
            The centroid of the bonder attoms connected to those on the
            edge.

        """

        coord = np.zeros((3, ))
        for i, position in enumerate(self.joint_positions):
            vertex = self.connected[i]
            bonder = vertex.bonder_map[position]
            coord += macro_mol.atom_coords(bonder)

        for d, param in zip(self.bond, cell_params):
            coord += d*param

        return coord / (i+1)

    def create_bonder_map(self, macro_mol, cell_params):
        """
        Creates the attribute :attr:`bonder_map`.

        Parameters
        ----------
        macro_mol : :class:`.MacroMolecule`
            The macromolecule being built.

        cell_params : :class:`list` of :class:`numpy.array`
            The ``a``, ``b`` and ``c`` vectors of the unit cell.

        Returns
        -------
        None : :class:`NoneType`

        """

        # Get all the bonder ids.
        bonder_ids = deque(maxlen=2)
        for atom in macro_mol.mol.GetAtoms():
            if atom.HasProp('bonder'):
                bonder_ids.append(atom.GetIdx())

        v1coord = self.v1.calc_coord(cell_params)
        bonders = sorted(bonder_ids,
                         key=lambda x: euclidean(
                                             v1coord,
                                             macro_mol.atom_coords(x)))
        self.bonder_map = {0: bonders[0],
                           1: bonders[1]}


def is_bonder(macro_mol, atom_id):
    """
    ``True`` if atom has ``'bonder'`` property.

    Parameters
    ----------
    macro_mol : :class:`.MacroMolecule`
        The macromolecule to which the atom belongs.

    Returns
    -------
    :class:`bool`
        ``True`` if atom has ``'bonder'`` property else ``False``.

    """

    return (True if
            macro_mol.mol.GetAtomWithIdx(atom_id).HasProp('bonder')
            else False)


def bb_size(macro_mol):
    """
    Sums the diameters of the building blocks.

    Parameters
    ---------
    macro_mol : :class:`.MacroMolecule`
        The macromolecule being assembled.

    Returns
    -------
    :class:`float`
        The sum of the diameters of the building blocks. Used to
        scale the size of unit cells.

    """

    return sum(bb.max_diameter()[0] for bb in macro_mol.building_blocks)


def linker_cof_scale_func(macro_mol):
    """
    Calculates the scale factor for :class:`LinkerCOFLattice`.

    Parameters
    ---------
    macro_mol : :class:`.MacroMolecule`
        The macromolecule being assembled.

    Returns
    -------
    :class:`float`
        A product of the sum of the diameters of the building blocks
        and the number of vertices in the unit cell.

    """

    return bb_size(macro_mol)*((len(macro_mol.topology.vertices)+1)//1.5)


class COFLattice(Topology):
    """
    A base class for periodic topologies.

    This class behaves almost exactly like
    :class:`.Topology` with only minor additions to suit the
    representation of COF periodic lattices. The :meth:`del_atoms`
    method is extended to collect the coordinates of any deleted atoms.
    This is necessary for positioning terminating atoms on islands
    generated from the periodic structure by :meth:`.Periodic.island`.

    Attributes
    ----------
    scale_func : :class:`function`
        A function which takes a :class:`.MacroMolecule` as its only
        argument returns a number. The number is used
        to scale the size of the unit cell. See :func:`bb_size`
        for an example.

    """

    def __init__(self, scale_func=bb_size):
        self.scale_func = scale_func
        super().__init__()

    def del_atoms(self, macro_mol):
        """
        Deletes the atoms which are lost during assembly.

        Notes
        -----
        The parameter `macro_mol` has two attributes changed.
        :attr:`~.MacroMolecule.mol` has deleter atoms removed, while
        :attr:`~.Periodic.deleters` is updated with the
        coordinates, element type and bond type of every removed atom.

        Parameters
        ----------
        macro_mol : :class:`.Periodic`
            The periodic macromolecule being assembled.

        Returns
        -------
        None : :class:`NoneType`

        """

        # This loop checks if an atom is an bonder and if its is finds
        # a neighboring deleter atom. The deleter atom has its
        # coordinates relative to the bonder found and saved in
        # `terminator_coords`.
        macro_mol.deleters = defaultdict(list)
        for atom in macro_mol.mol.GetAtoms():
            if not atom.HasProp('bonder'):
                continue
            for neighbor in atom.GetNeighbors():
                if not neighbor.HasProp('del'):
                    continue
                bi = macro_mol.bonder_ids.index(atom.GetIdx())
                nid = neighbor.GetIdx()
                tcoords = (macro_mol.atom_coords(nid) -
                           macro_mol.atom_coords(atom.GetIdx()))

                bond = macro_mol.mol.GetBondBetweenAtoms(nid, atom.GetIdx())
                macro_mol.deleters[bi].append([tcoords,
                                               neighbor.GetAtomicNum(),
                                               bond.GetBondType()])

        super().del_atoms(macro_mol)
        macro_mol._ids_updated = False


class LinkerCOFLattice(COFLattice):
    """
    For 2D COF topologies which use a linker.

    These are topologies which consist of 2 building blocks. One
    which has two functional groups and one which has 3 or more
    functional groups.

    This class will have to be extended by sublcasses which define
    the vertices and edges which a particular topology consists of.

    Attributes
    ----------
    vertices : :class:`list` of :class:`Vertex`
        Class attribute added when defining a subclass. Holds the
        vertices which make up the topology.

    edges : :class:`list` of :class:`Edge`
        Class attribute added when defining a subclass. Holds the
        edges which make up the topology.

    ditopic_directions : :class:`list` of :class:`int`
        For each edge in the topology, holds ``1`` or ``-1`` to
        indicate if the building block is placed parallel or
        anti-parallel to the edge.

    multitopic_aligners : :class:`list` of :class:`int`
        For each vertex in the topology, holds the index of an
        atom in :attr:`.StructUnit.bonder_ids`. It specicifes which
        bonder atom is used for alignment of the building block on the
        vertex.

    """

    def __init__(self,
                 ditopic_directions=None,
                 multitopic_aligners=None,
                 scale_func=linker_cof_scale_func):
        super().__init__(scale_func=scale_func)

        if ditopic_directions is None:
            ditopic_directions = [1 for i in range(len(self.edges))]
        if multitopic_aligners is None:
            multitopic_aligners = [0 for i in range(len(self.vertices))]

        self.ditopic_directions = ditopic_directions
        self.multitopic_aligners = multitopic_aligners

    def place_mols(self, macro_mol):
        """
        Places the building blocks on the topology.

        Parameters
        ----------
        macro_mol : :class:`.MacroMolecule`
            The macromolecule being assembled.

        Returns
        -------
        None : :class:`NoneType`

        """

        # Create the rdkit molecule of the assembled molecule.
        macro_mol.mol = rdkit.Mol()

        # Identify which building block is ditopic and which is
        # tri or more topic.
        di = next(bb for bb in macro_mol.building_blocks if
                  len(bb.functional_group_atoms()) == 2)
        multi = next(bb for bb in macro_mol.building_blocks if
                     len(bb.functional_group_atoms()) >= 3)
        nbonders = len(multi.functional_group_atoms())

        # Calculate the size of the unit cell by scaling to the size of
        # building blocks.
        size = self.scale_func(macro_mol)
        cell_params = [size*p for p in self.cell_dimensions]
        macro_mol.cell_dimensions = cell_params

        # For each vertex in the topology, place a multitopic building
        # block on it. The Vertex object takes care of alignment.

        for i, v in enumerate(self.vertices):
            aligner = self.multitopic_aligners[i]
            mol = v.place_mol(cell_params, multi, aligner)
            add_fragment_props(mol,
                               macro_mol.building_blocks.index(multi),
                               i)
            macro_mol.mol = rdkit.CombineMols(macro_mol.mol, mol)
            macro_mol.bb_counter.update([multi])

            # Save the ids of the bonder atoms in the assembled molecule.
            # This is used when creating bonds later in the assembly
            # process.
            v.create_bonder_map(macro_mol, cell_params, nbonders, aligner)

        for i, e in enumerate(self.edges):
            mol = e.place_mol(macro_mol,
                              cell_params,
                              di,
                              self.ditopic_directions[i])
            add_fragment_props(mol,
                               macro_mol.building_blocks.index(di),
                               i)
            macro_mol.mol = rdkit.CombineMols(macro_mol.mol, mol)
            macro_mol.bb_counter.update([di])
            bonder_ids = deque(maxlen=2)
            for atom in macro_mol.mol.GetAtoms():
                if atom.HasProp('bonder'):
                    bonder_ids.append(atom.GetIdx())

            e.create_bonder_map(macro_mol, cell_params)

        super(macro_mol.__class__, macro_mol).save_ids()

    def join_mols(self, macro_mol):
        """
        Creates bonds between the building blocks.

        Parameters
        ----------
        macor_mol : :class:`.MacroMolecule`
            The macromolecule being assembled.

        Returns
        -------
        None : :class:`NoneType`

        """

        macro_mol.bonds_made = 0
        emol = rdkit.EditableMol(macro_mol.mol)
        for e in self.edges:
            bonder1 = e.bonder_map[0]
            bonder2 = e.v1.bonder_map[e.joint_positions[0]]
            bond_type = self.determine_bond_type(macro_mol,
                                                 bonder1,
                                                 bonder2)
            emol.AddBond(bonder1, bonder2, bond_type)

            bonder3 = e.bonder_map[1]
            bonder4 = e.v2.bonder_map[e.joint_positions[1]]
            if all(b == 0 for b in e.bond):
                bond_type = self.determine_bond_type(macro_mol,
                                                     bonder3,
                                                     bonder4)
                emol.AddBond(bonder3, bonder4, bond_type)
            else:
                macro_mol.periodic_bonds.append(
                    PeriodicBond(macro_mol.bonder_ids.index(bonder3),
                                 macro_mol.bonder_ids.index(bonder4),
                                 e.bond))

            macro_mol.bonds_made += 2

        macro_mol.mol = emol.GetMol()


class NoLinkerCOFLattice(COFLattice):
    ...


class NoLinkerHoneycomb(NoLinkerCOFLattice):
    """
    Represents a hexagonal lattice with 2 tritopic building blocks.

    """

    cell_dimensions = a, b, c = [np.array([1, 0, 0]),
                                 np.array([0.5, 0.866, 0]),
                                 np.array([0, 0, 5/1.7321])]

    vertices = [(a/3 + b/3 + c/2),
                (2*a/3 + 2*b/3 + c/2)]

    def join_mols(self, macro_mol):
        """
        Joins the building blocks in the unit cell.

        Notes
        -----
        The rdkit instance in the :attr:`~.MacroMolecule.mol` attribute
        has bonds added to it. The
        :attr:`~.MacroMolecule.periodic_bonds` attribute is also filled
        with :class:`PeriodicBond` instances.

        Parameters
        ----------
        macro_mol : :class:`.MacroMolecule`
            The unit cell being constructed.

        Returns
        -------
        None : :class:`NoneType`

        """

        # Get the fragments.
        frag1, frag2 = rdkit.GetMolFrags(macro_mol.mol,
                                         sanitizeFrags=False)
        # Get rid of any non bonder atoms.
        frag1 = [x for x in frag1 if is_bonder(macro_mol, x)]
        frag2 = [x for x in frag2 if is_bonder(macro_mol, x)]

        # The fragment with the larger bonder ids has higher x and y
        # values - due to place_mols() implmentation. It is the "top"
        # fragment.
        top = frag1 if frag1[0] > frag2[0] else frag2
        bottom = frag2 if top is frag1 else frag1

        # In the top fragment find the bonder atom with the
        # largest y value and connect it to the bonder atom in the
        # bottom fragment with the lowest y value. Note that the
        # connection must be registered as periodic, hence the
        # directions are 1/-1.
        top_atom = max(top, key=lambda x: macro_mol.atom_coords(x)[1])
        bottom_atom = min(bottom,
                          key=lambda x: macro_mol.atom_coords(x)[1])
        # The bonder atoms are registered by their index within
        # `bonder_ids`. This is because del_atoms will change the
        # atom ids but not the ordering of this list.
        top_atom = macro_mol.bonder_ids.index(top_atom)
        bottom_atom = macro_mol.bonder_ids.index(bottom_atom)
        macro_mol.periodic_bonds.append(
                    PeriodicBond(top_atom, bottom_atom, [0, 1, 0]))
        # Do the same for the x-axis periodic bonds.
        right_atom = max(top,
                         key=lambda x: macro_mol.atom_coords(x)[0])
        left_atom = min(bottom,
                        key=lambda x: macro_mol.atom_coords(x)[0])
        # The bonder atoms are registered by their index within
        # `bonder_ids`. This is because del_atoms will change the
        # atom ids.
        right_atom = macro_mol.bonder_ids.index(right_atom)
        left_atom = macro_mol.bonder_ids.index(left_atom)
        macro_mol.periodic_bonds.append(
                        PeriodicBond(right_atom, left_atom, [1, 0, 0]))

        # For the bond which gets created directly, find the bonder
        # atom in the bottom fragment closest to the position of the
        # top fragment. Create a bond between it and the bonder atom
        # in the top fragment closest to the position of the bottom
        # fragment.
        bottom_bonder = min(bottom, key=lambda x:
                            euclidean(self.vertices[1],
                                      macro_mol.atom_coords(x)))
        top_bonder = min(top, key=lambda x:
                         euclidean(self.vertices[0],
                                   macro_mol.atom_coords(x)))

        emol = rdkit.EditableMol(macro_mol.mol)
        bond_type = self.determine_bond_type(macro_mol,
                                             top_bonder,
                                             bottom_bonder)
        emol.AddBond(top_bonder, bottom_bonder, bond_type)
        macro_mol.mol = emol.GetMol()

    def place_mols(self, macro_mol):
        """
        Places and aligns building blocks in the unit cell.

        Notes
        -----
        This method modifies `macro_mol`. An rdkit molecule of the
        unit cell with the building blocks not joined up is placed in
        the :attr:`~.MacroMolecule.mol` attribute.

        Parameters
        ----------
        macro_mol : :class:`.MacroMolecule`
            The unit cell being constructed.

        Returns
        -------
        None : :class:`NoneType`

        """

        # Make the rdkit molecule.
        macro_mol.mol = rdkit.Mol()
        # Get the building blocks.
        bb1, bb2 = macro_mol.building_blocks
        cell_size = self.scale_func(macro_mol)
        macro_mol.cell_dimensions = [cell_size*x for x in
                                     self.cell_dimensions]
        self.vertices = [cell_size*x for x in self.vertices]
        # Place and set orientation of the first building block.
        bb1.set_bonder_centroid(self.vertices[0])
        bb1.set_orientation2([0, 0, 1])
        bb1.minimize_theta2(bb1.bonder_ids[0], [0, -1, 0], [0, 0, 1])
        # Add to the macromolecule.
        add_fragment_props(bb1.mol,
                           macro_mol.building_blocks.index(bb1),
                           0)
        macro_mol.mol = rdkit.CombineMols(macro_mol.mol, bb1.mol)
        # Place and set orientation of the second building block.
        bb2.set_bonder_centroid(self.vertices[1])
        bb2.set_orientation2([0, 0, 1])
        bb2.minimize_theta2(bb2.bonder_ids[0], [0, 1, 0], [0, 0, 1])
        # Add to the macromolecule.
        add_fragment_props(bb2.mol,
                           macro_mol.building_blocks.index(bb2),
                           0)
        macro_mol.mol = rdkit.CombineMols(macro_mol.mol, bb2.mol)
        # Add the bonder_ids prematurely for this topology. Needed for
        # making supercells - see join_mols(). Using the ``super``
        # version here because Periodic.save_ids() tries to update the
        # atom ids in periodic bonds, which is needed later in the
        # assembly process. Here only saving of the atom ids
        # is needed, which is done by the ``super`` version.
        super(macro_mol.__class__, macro_mol).save_ids()


class Honeycomb(LinkerCOFLattice):
    cell_dimensions = a, b, c = [np.array([1, 0, 0]),
                                 np.array([0.5, 0.866, 0]),
                                 np.array([0, 0, 5/1.7321])]

    vertices = v1, v2 = [Vertex((1/3, 1/3, 1/2)),
                         Vertex((2/3, 2/3, 1/2))]

    edges = [Edge(v1, v2, (0, 2)),
             Edge(v1, v2, (1, 0), [0, -1, 0]),
             Edge(v1, v2, (2, 1), [-1, 0, 0])]


class Hexagonal(LinkerCOFLattice):
    cell_dimensions = a, b, c = [np.array([1, 0, 0]),
                                 np.array([0.5, 0.866, 0]),
                                 np.array([0, 0, 5/1.7321])]

    vertices = v1, v2, v3, v4 = [Vertex((1/4, 1/4, 1/2)),
                                 Vertex((1/4, 3/4, 1/2)),
                                 Vertex((3/4, 1/4, 1/2)),
                                 Vertex((3/4, 3/4, 1/2))]

    edges = [Edge(v1, v2, (0, 3)),
             Edge(v1, v3, (1, 4)),
             Edge(v2, v3, (2, 5)),
             Edge(v2, v4, (1, 4)),
             Edge(v3, v4, (0, 3)),
             Edge(v1, v3, (4, 1), [-1, 0, 0]),
             Edge(v1, v2, (3, 0), [0, -1, 0]),
             Edge(v1, v4, (2, 5), [0, -1, 0]),
             Edge(v3, v2, (2, 5), [1, -1, 0]),
             Edge(v3, v4, (3, 0), [0, -1, 0]),
             Edge(v2, v4, (4, 1), [-1, 0, 0]),
             Edge(v4, v1, (2, 5), [1, 0, 0])]


class Square(LinkerCOFLattice):
    cell_dimensions = a, b, c = [np.array([1, 0, 0]),
                                 np.array([0, 1, 0]),
                                 np.array([0, 0, 1])]

    vertices = v1, = [Vertex((0.5, 0.5, 0.5))]
    edges = [Edge(v1, v1, (1, 3), [1, 0, 0]),
             Edge(v1, v1, (0, 2), [0, 1, 0])]


class Kagome(LinkerCOFLattice):
    cell_dimensions = a, b, c = [np.array([1, 0, 0]),
                                 np.array([0.5, 0.866, 0]),
                                 np.array([0, 0, 5/1.7321])]

    vertices = v1, v2, v3 = [Vertex((1/4, 3/4, 0.5)),
                             Vertex((3/4, 3/4, 1/2)),
                             Vertex((3/4, 1/4, 1/2))]

    edges = [Edge(v1, v2, (0, 3)),
             Edge(v1, v3, (1, 3)),
             Edge(v2, v3, (2, 0)),
             Edge(v1, v2, (2, 1), [-1, 0, 0]),
             Edge(v1, v3, (3, 1), [-1, 1, 0]),
             Edge(v2, v3, (0, 2), [0, 1, 0])]
