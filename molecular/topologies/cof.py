import rdkit.Chem.AllChem as rdkit
import numpy as np
from scipy.spatial.distance import euclidean
from collections import deque

from .base import Topology
from ...convenience_tools import PeriodicBond, add_fragment_props


class Edge:
    """
    Represents a bond in a COF unit cell.

    This class is used for positioning ditopic building blocks.

    Attributes
    ----------
    v1 : :class:`Vertex`


    """

    def __init__(self, v1, v2, bond=[0, 0, 0]):
        self.v1 = v1
        self.v2 = v2
        self.bond = bond


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
        macro_mol.terminator_coords = {}
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
                macro_mol.terminator_coords[bi] = tcoords

        super().del_atoms(macro_mol)
        macro_mol._ids_updated = False


class LinkerCOFLattice(COFLattice):
    def place_mols(self, macro_mol):

        # Create the rdkit molecule of the assembled molecule.
        macro_mol.mol = rdkit.Mol()

        # Identify which building block is ditopic and which is
        # tritopic.
        di = next(bb for bb in macro_mol.building_blocks if
                  len(bb.functional_group_atoms()) == 2)
        tri = next(bb for bb in macro_mol.building_blocks if
                   len(bb.functional_group_atoms() == 3))

        # For each vertex in the topology, place a tritopic building
        # block on it. The Vertex object takes care of alignment.

        for i, v in enumerate(self.vertices):
            mol = v.place_mol(tri,
                              self.tritopic_aligners[i],
                              self.edge_alignments[i])
            add_fragment_props(mol,
                               macro_mol.building_blocks.index(tri),
                               i)
            macro_mol.mol = rdkit.CombineMols(macro_mol.mol, mol)
            macro_mol.bb_counter.update([tri])

            # Save the ids of the bonder atoms in the assembled molecule.
            # This is used when creating bonds later in the assembly
            # process.
            bonder_ids = deque(maxlen=3)
            for atom in macro_mol.GetAtoms():
                if atom.HasProp('bonder'):
                    bonder_ids.append(atom.GetIdx())
            v.bonder_ids = sorted(bonder_ids)

        for i, e in enumerate(self.edges):
            mol = e.place_mol(di, self.ditopic_directions[i])
            add_fragment_props(mol,
                               macro_mol.building_blocks.index(di),
                               i)
            macro_mol.mol = rdkit.CombineMols(macro_mol.mol, mol)
            macro_mol.bb_counter.update([di])
            bonder_ids = deque(maxlen=2)
            for atom in macro_mol.mol.GetAtoms():
                if atom.HasProp('bonder'):
                    bonder_ids.append(atom.GetIdx())

            e.bonder_ids = sorted(bonder_ids)

    def join_mols(self, macro_mol):
        emol = rdkit.EditableMol(macro_mol.mol)
        macro_mol.bonds_made = 0
        

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

    vertices = v1, v2 = [(a/3 + b/3 + c/2),
                         (2*a/3 + 2*b/3 + c/2)]

    edges = [Edge(v1, v2),
             Edge(v1, v2, [0, -1, 0]),
             Edge(v1, v2, [-1, 0, 0])]


class Hexagonal(LinkerCOFLattice):
    cell_dimensions = a, b, c = [np.array([1, 0, 0]),
                                 np.array([0.5, 0.866, 0]),
                                 np.array([0, 0, 5/1.7321])]

    vertices = v1, v2 = [(a/3 + b/3 + c/2),
                         (2*a/3 + 2*b/3 + c/2)]

    edges = [Edge(v1, v2),
             Edge(v1, v1, [1, 0, 0]),
             Edge(v2, v2, [1, 0, 0]),
             Edge(v1, v2, [-1, -1, 0]),
             Edge(v1, v2, [-1, 0, 0]),
             Edge(v1, v2, [0, -1, 0])]


class Square(LinkerCOFLattice):
    cell_dimensions = a, b, c = [np.array([1, 0, 0]),
                                 np.array([0, 1, 0]),
                                 np.array([0, 0, 1])]

    vertices = v1, = [(a/2 + b/2 + c/2)]
    edges = [Edge(v1, v1, [1, 0, 0]),
             Edge(v1, v1, [0, 1, 0])]


class Kagome(LinkerCOFLattice):
    cell_dimensions = a, b, c = [np.array([1, 0, 0]),
                                 np.array([0.5, 0.866, 0]),
                                 np.array([0, 0, 5/1.7321])]

    vertices = v1, v2, v3 = [(a/3 + 2*b/3 + c/2),
                             (2*a/3 + 2*b/3 + c/2),
                             (2*a/3 + b/3 + c/2)]

    edges = [Edge(v1, v2),
             Edge(v1, v3),
             Edge(v2, v3),
             Edge(v1, v2, [-1, 0, 0]),
             Edge(v1, v3, [-1, 1, 0]),
             Edge(v2, v3, [0, 1, 0])]
