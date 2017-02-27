"""
Defines the base ``Topology`` type.

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


"""

import rdkit
import rdkit.Chem as chem
from collections import deque
import numpy as np

from ..fg_info import double_bond_combs
from ...convenience_tools import dedupe

class Topology:
    """
    Builds macromolecules.

    More accurately, child classes of ``Topology`` class take care of
    building macromolecules.

    This class directly defines any operations and attributes that are
    needed by any topology during assembly. However, this class is not
    used directly by MMEA. It is intended to be inherited from. All
    macromolecules within MMEA have a `topology` attribute which hold
    an instance of a Topology child class. Child classes of Topology
    define operations specific to that one topology. For example, each
    child class must define a `join_mols()` method which creates bonds
    between the building blocks of a macromolecule. The way in which
    this is done will depend on what kind of macromolecules are being
    built. In addition, each child class must define methods which
    place the building blocks in approriate positions.

    """

    def build(self, macro_mol):
        """
        Assembles rdkit instances of macromolecules.

        Parameters
        ----------
        macro_mol : MacroMolecule
            The MacroMolecule instance which needs to be built.

        Modifies
        --------
        macro_mol.mol : rdkit.Chem.rdchem.Mol
            An rdkit instance of the assembled macromolecule is placed
            in this attribute.

        macro_mol.bb_counter : Counter
            The counter is updated with the amounts of building blocks
            used to make the macromolecule.

        macro_mol.bonds_made : int
            This count is updated with the number of bonds made during
            assembly.

        Returns
        -------
        None : NoneType

        """

        # When running ``build()`` in parallel, the atom tags are
        # cleared by the multiprocessing module. Make sure to readd the
        # tags before running ``build()``.
        for bb in macro_mol.building_blocks:
            bb.tag_atoms()

        self.place_mols(macro_mol)
        self.join_mols(macro_mol)
        self.del_atoms(macro_mol)

    def del_atoms(self, macro_mol):
        """
        Deletes the atoms which are lost during assembly.

        Modifies
        --------
        macro_mol.mol : rdkit.Chem.rdchem.Mol
            Atoms are removed from this molecule.

        Returns
        -------
        None : NoneType

        """

        mol = chem.EditableMol(macro_mol.mol)
        # Delete atoms with largest id last, so that when deleting
        # later atoms their ids do not change.
        for atom in reversed(macro_mol.mol.GetAtoms()):
            if atom.HasProp('del'):
                mol.RemoveAtom(atom.GetIdx())

        macro_mol.mol = mol.GetMol()

    @staticmethod
    def determine_bond_type(macro_mol, atom1_id, atom2_id):
        """
        Returns the bond order to be formed between the atoms.

        Some atoms will need to have a double bond created between
        them. This is defined in the `double_bond_combs` list. If the
        atom ids provided as paramters belong to functional grups found
        in this list, the rdkit double bond type will be returned.
        If not the rdkit single bond type will be returned. These types
        are needed when adding bonds using ``EditableMol`` instances.

        Parameters
        ----------
        macro_mol : MacroMolecule
            The macromolecule being assembled.

        atom1_id : int
            The id number of the first atom.

        atom2_id : int
            The id number of the second atom.

        Returns
        -------
        rdkit.Chem.rdchem.BondType.SINGLE
            If the atoms don't belong to functional groups which form
            a double bond.

        rdkit.Chem.rdchem.BondType.DOUBLE
            If the atoms belong to functional groups which form a
            double bond.

        """

        # Get the functional groups of the of the atoms whose atom ids
        # were supplied as arguments. If the groups form a tuple in
        # `double_bond_combs` return a rdkit double bond type. If they
        # do not, return a rdkit single bond type.

        atom1 = macro_mol.mol.GetAtomWithIdx(atom1_id)
        atom1_grp = atom1.GetProp('fg')
        atom2 = macro_mol.mol.GetAtomWithIdx(atom2_id)
        atom2_grp= atom2.GetProp('fg')

        double_bond_present = ((atom1_grp, atom2_grp) == tup or
                               (atom2_grp, atom1_grp) == tup for
                               tup in double_bond_combs)

        if any(double_bond_present):
            return rdkit.Chem.rdchem.BondType.DOUBLE
        else:
            return rdkit.Chem.rdchem.BondType.SINGLE

    def __str__(self):
        return repr(self)

    def __repr__(self):
        c = ""
        for key, value in sorted(self.__dict__.items()):
            c += "{!s}={!r}".format(key, value)
        return "{}({})".format(self.__class__.__name__, c)

    def __eq__(self, other):
        return repr(self) == repr(other)

    def __hash__(self):
        return id(self)

class Linear(Topology):
    """
    A class represting linear polymers.

    Attributes
    ----------
    repeating_unit : str
        A string showing the repeating unit of the Polymer. For
        example, "AB" or "ABB". The building block with index 0 in
        `building-blocks` is labelled as "A" while index 1 as "B" and
        so on.

    orientation : list of ints
        For each character in the repeating unit, a value of -1, 0 or
        1 must be given as a list. It indicates the direction at
        which each monomer of the repeating unit is placed. 0 means
        that the direction is random.

    n : int
        The number of repeating units which are used to make the
        polymer.

    """

    def __init__(self, repeating_unit, orientation, n):
        self.repeating_unit = repeating_unit
        self.orientation = orientation
        self.n = n

    def place_mols(self, macro_mol):
        """
        Places monomers side by side and joins them.

        The monomers are placed along the x-axis, so that the vector
        running between the functional groups is placed on the axis.

        Parameters
        ----------
        macro_mol : Polymer
            The polymer being assembled.

        Modifies
        --------
        macro_mol.mol : rdkit.Chem.rdchem.Mol
            The monomers are placed.

        Returns
        -------
        None : NoneType

        """

        # Make a map from monomer label to object.
        mapping = {}
        # Assign every monomer a label ("A", "B", "C", etc.).
        for label, monomer in zip(dedupe(self.repeating_unit),
                                  macro_mol.building_blocks):
            mapping[label] = monomer

        # Go through the repeating unit. Place each monomer 50 A apart.
        # Also create a bond.
        self.bonders = deque(maxlen=2)
        macro_mol.mol = chem.Mol()
        polymer = self.repeating_unit*self.n
        dirs = ",".join(str(x) for x in self.orientation)
        dirs *= self.n
        dirs = [int(x) for x in orientations.split(',')]

        for i, (label, mdir) in enumerate(zip(polymer, dirs)):
            self.bonders.append([
                macro_mol.mol.GetNumAtoms() + id_ for
                           id_ in mapping[label].bonder_ids])

            mdir = np.random.choice([-1, 1]) if not mdir else mdir
            mapping[label].set_orientation2([mdir, 0, 0])
            macro_mol.mol = chem.CombineMols(macro_mol.mol,
                            mapping[label].set_position([i*50, 0, 0]))
            if i != 0:
                self.join(macro_mol)

    def join(self, macro_mol):
        """
        Joins 2 monomers.

        Parameters
        ----------
        macro_mol : Polymer
            The polymer being assembled.

        Modifies
        --------
        macro_mol.mol : rdkit.Chem.rdchem.Mol
            The polymer with the monomers connected.

        Returns
        -------
        None : NoneType

        """

        distances = []
        for bonder1 in self.bonders[0]:
            for bonder2 in self.bonders[1]:
                distances.append(
        (macro_mol.atom_distance(bonder1, bonder2), bonder1, bonder2))

        _, bonder1, bonder2 = min(distances)
        emol = chem.EditableMol(macro_mol.mol)
        emol.AddBond(bonder1, bonder2,
                self.determine_bond_type(macro_mol, bonder1, bonder2))

        macro_mol.mol = emol.GetMol()

    def join_mols(self, macro_mol):
        """
        Does nothing, place_mols() joins up the molecules too.

        Parameters
        ----------
        macro_mol : Polymer
            The polymer being assembled.

        Returns
        -------
        None : NoneType

        """

        return
