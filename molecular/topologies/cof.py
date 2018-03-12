import rdkit.Chem.AllChem as rdkit
import numpy as np
from scipy.spatial.distance import euclidean
from collections import deque, defaultdict

from .base import Topology
from ...convenience_tools import (PeriodicBond,
                                  add_fragment_props,
                                  normalize_vector)


class Vertex:
    def __init__(self, frac_coord):
        self.frac_coord = frac_coord
        self.connected = []

    def place_mol(self, cell_params, mol, aligner):
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
        coord = np.zeros((3, ))
        for frac, dim, p in zip(self.frac_coord, cell_params, cell_position):
            coord += (frac+p) * dim

        return coord

    def create_bonder_map(self,
                          macro_mol,
                          cell_params,
                          nbonders,
                          aligned_bonder):
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

        aligned_pos = self.aligned_position(aligned_bonder)
        positions = (list(range(aligned_pos, nbonders)) +
                     list(range(0, aligned_pos)))

        self.bonder_map = {}
        for bonder, position in zip(bonders, positions):
            self.bonder_map[position] = bonder

    def aligned_position(self, aligned_bonder):
        aligner_edge = next((e for e in self.connected if
                             all(b == 0 for b in e.bond)),
                            self.connected[0])
        vindex = 0 if self is aligner_edge.v1 else 1
        return aligner_edge.joint_positions[vindex]


class Edge:
    """
    Represents a bond in a COF unit cell.

    This class is used for positioning ditopic building blocks.

    Attributes
    ----------
    v1 : :class:`Vertex`


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

    def place_mol(self, cell_params, mol, alignment):
        coord = self.calc_coord(cell_params)
        original_position = mol.position_matrix()

        mol.set_bonder_centroid(coord)
        mol.set_orientation2(self.direction(cell_params)*alignment)

        rdkit_mol = rdkit.Mol(mol.mol)
        mol.set_position_from_matrix(original_position)
        return rdkit_mol

    def direction(self, cell_params):
        v1coord = self.v1.calc_coord(cell_params)
        v2coord = self.v2.calc_coord(cell_params, self.bond)
        return normalize_vector((v1coord-v2coord)/2)

    def calc_coord(self, cell_params):
        coord = np.zeros((3, ))
        for frac, dim in zip(self.frac_coord, cell_params):
            coord += frac * dim
        return coord

    def create_bonder_map(self, macro_mol, cell_params):
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


class COFLattice(Topology):
    """
    A base class for periodic topologies.

    This class behaves almost exactly like
    :class:`.Topology` with only minor additions to suit the
    representation of COF periodic lattices. The :meth:`del_atoms`
    method is extended to collect the coordinates of any deleted atoms.
    This is necessary for positioning terminating atoms on islands
    generated from the periodic structure by :meth:`.Periodic.island`.

    """

    def del_atoms(self, macro_mol):
        """
        Deletes the atoms which are lost during assembly.

        Notes
        -----
        The parameter `macro_mol` has two attributes changed.
        :attr:`~.MacroMolecule.mol` has deleter atoms removed, while
        :attr:`~.Periodic.terminator_coords` is updated with the
        coordinates of every removed atom.

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
        macro_mol.terminator_coords = defaultdict(list)
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
                macro_mol.terminator_coords[bi].append(tcoords)

        super().del_atoms(macro_mol)
        macro_mol._ids_updated = False


class LinkerCOFLattice(COFLattice):
    def __init__(self,
                 ditopic_directions=None,
                 multitopic_aligners=None):
        super().__init__()

        if ditopic_directions is None:
            ditopic_directions = [1 for i in range(len(self.edges))]
        if multitopic_aligners is None:
            multitopic_aligners = [0 for i in range(len(self.vertices))]

        self.ditopic_directions = ditopic_directions
        self.multitopic_aligners = multitopic_aligners

    def place_mols(self, macro_mol):

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
        size = di.max_diameter()[0] + multi.max_diameter()[0]
        size *= len(self.vertices)+1
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
            mol = e.place_mol(cell_params, di, self.ditopic_directions[i])
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
        cell_size = bb1.max_diameter()[0] + bb2.max_diameter()[0]
        self.cell_dimensions = [cell_size*x for x in
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
             Edge(v2, v4, (4, 1), [-1, 0, 0])]


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
