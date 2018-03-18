"""
Defines the base :class:`Topology` type.

.. _`adding topologies`:

Extending stk: Adding new topologies.
-------------------------------------

Cages
.....

To add a new cage topology a new class should be created, named
after the topology. This class should inherit :class:`._CageTopology`.
This will give access to various methods which are necessary
for dealing with any cage molecule. See the documenation of
:class:`._CageTopology` for more details.

The new class will only need to have five class attributes added:

    1. a :class:`list` called :attr:`vertices`
    2. a :class:`list` called :attr:`edges`
    3. :attr:`n_windows`, which holds the number of windows the cage
       topology has.
    4. :attr:`n_window_types`, which holds the number of different
       window types. For example, if :attr:`n_window_types` is ``2``,
       then the topology will have two kinds of windows, each with a
       different expected size. Windows of the same type are expected
       to be of the same size.

:attr:`vertices` holds instances of :class:`~.cage.base.Vertex`. Each
instance represents a vertex of a cage and needs to be initialized
with the coordinates of that vertex. Vertices of a cage are where
building blocks of cages are placed.

:attr:`edges` holds instances of the :class:`~.cage.base.Edge`. Each
instance represents an edge of a cage and needs to be initialized
with two instances of :class:`~.cage.base.Vertex`. The
:class:`~.cage.base.Vertex` instances
should be held in :attr:`vertices`, as mentioned above. The two
vertices are the ones which the edge connects. Linkers of cages are
placed on edges. The edge instances automatically derive their
positions from the vertices supplied during initialization.

The vertices need to be positioned such that the center of the
topology is at the origin.


"""

import rdkit.Chem.AllChem as rdkit
from collections import deque
import numpy as np
from itertools import chain
from inspect import signature

from ..fg_info import double_bond_combs
from ...convenience_tools import dedupe, flatten, add_fragment_props


def remove_confs(building_blocks, keep):
    """
    Removes all conformers from `building_blocks` except `keep`.

    All kept conformers have their id set to 0.

    Parameters
    ----------
    building_blocks : iterable of ``StructUnit`` instances
        A set of ``StructUnit`` instances which represent the
        building blocks forming a macromolecule.

    keep : :class:`list` of :class:`int`
        The ids of the building block conformers to be used for
        assembling the :class:`.MacroMolecule`. Must be equal in length
        to `building_blocks` and orders must correspond.

    Returns
    -------
    :class:`list`
        A :class:`list` of the form,

        .. code-block:: python

        returned = [[conf1, conf2, conf3],
                    [conf4, conf5],
                    [conf6, conf7, conf8, conf9]]

        where each sublist holds all the original conformers of a
        particular building block.

    """

    keep_ids = [bb.mol.GetConformer(id_).GetId() for
                bb, id_ in zip(building_blocks, keep)]

    original_confs = [[rdkit.Conformer(conf) for
                       conf in bb.mol.GetConformers()]
                      for bb in building_blocks]
    for bb, conf in zip(building_blocks, keep_ids):
        keep_conf = rdkit.Conformer(bb.mol.GetConformer(conf))
        keep_conf.SetId(0)
        bb.mol.RemoveAllConformers()
        bb.mol.AddConformer(keep_conf)
    return original_confs


class TopologyMeta(type):
    """
    Makes a repr of an instance, based initialization arguments used.

    """

    def __call__(self, *args, **kwargs):

        # Get the arguments, keyword arguments and defulat initialized
        # arguments used to make an instance of Topology.
        sig = signature(self.__init__).bind_partial(self,
                                                    *args, **kwargs)
        sig.apply_defaults()
        sig = dict(sig.arguments)
        sig.pop('self')
        # Create the Topology instance.
        obj = super().__call__(*args, **kwargs)
        # Use the arguments the object was initialized with to make
        # a repr of the object and place it in the `repr` attribute.
        # The __repr__() function in Topology will then just return
        # this attribute.
        c = ', '.join("{!s}={!r}".format(key, value) for key, value in
                      sorted(sig.items()))
        obj._repr = "{}({})".format(self.__name__, c)
        return obj


class Topology(metaclass=TopologyMeta):
    """
    Builds macromolecules.

    More accurately, child classes of :class:`Topology` take care of
    building macromolecules.

    This class directly defines any operations and attributes that are
    needed by any topology during assembly. However, this class is not
    used directly. It is intended to be inherited from. All
    :attr:`.MacroMolecule.topology` attributes hold an instance of a
    :class:`Topology` child class. Child classes of :class:`Topology`
    define operations specific to that one topology. For example, each
    child class must define a :meth:`join_mols`, which creates bonds
    between the building blocks of a macromolecule. The way in which
    this is done will depend on what kind of macromolecules are being
    built. In addition, each child class must define methods which
    place the building blocks in approriate positions.

    """

    def build(self, macro_mol, bb_conformers=None):
        """
        Assembles ``rdkit`` instances of macromolecules.

        This method places an ``rdkit`` molecule of the assembled
        macromolecule into the :attr:`~.Molecule.mol` attribute of
        :class:`.MacroMolecule`. It also updates the
        :attr:`.MacroMolecule.bb_counter` and
        :attr:`.MacroMolecule.bonds_made` attributes.

        Parameters
        ----------
        macro_mol : :class:`.MacroMolecule`
            The :class:`.MacroMolecule` instance which needs to be
            built.

        bb_conformers : :class:`list` of :class:`int`, optional
            The ids of the building block conformers to be used. Must
            be equal in length to `building_blocks` and orders must
            correspond. If ``None``, then ``-1`` is used for all
            building blocks.

        Returns
        -------
        None : :class:`NoneType`

        """

        if bb_conformers is None:
            bb_conformers = [-1 for _ in
                             range(len(macro_mol.building_blocks))]

        # When running ``build()`` in parallel, the atom tags are
        # cleared by the multiprocessing module. Make sure to reapply
        # the tags before running ``build()``.
        for bb in macro_mol.building_blocks:
            bb.tag_atoms()

        # When building, only a single conformer should exist per
        # building block. Otherwise, rdkit.CombineMols won't work. It
        # only combines conformers with the same id.
        original_confs = remove_confs(macro_mol.building_blocks,
                                      bb_conformers)

        self.place_mols(macro_mol)
        self.join_mols(macro_mol)
        self.del_atoms(macro_mol)

        # Make sure that the property cache of each atom is up to date.
        for atom in macro_mol.mol.GetAtoms():
            atom.UpdatePropertyCache()

        # Make sure that the tags showing which building block each
        # atom belongs to are saved in `fragment_assignments` and not
        # lost due to parallelism.
        macro_mol.update_fragments()

        # Restore the original conformers.
        for bb, confs in zip(macro_mol.building_blocks, original_confs):
            bb.mol.RemoveAllConformers()
            for conf in confs:
                bb.mol.AddConformer(conf)

    def del_atoms(self, macro_mol):
        """
        Deletes the atoms which are lost during assembly.

        Parameters
        ----------
        macro_mol : :class:`.MacroMolecule`
            The macromolecule being assembled.

        Returns
        -------
        None : :class:`NoneType`

        """

        mol = rdkit.EditableMol(macro_mol.mol)
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
        them. This is defined in the :data:`fg_info.double_bond_combs`.
        If the atom ids provided as parameters belong to functional
        groups found in this list, the ``rdkit`` double bond type will
        be returned. If not, the ``rdkit`` single bond type will be
        returned. These types are needed when adding bonds to an
        ``rdkit`` molecule.

        Parameters
        ----------
        macro_mol : :class:`.MacroMolecule`
            The macromolecule being assembled.

        atom1_id : :class:`int`
            The id number of the first atom.

        atom2_id : :class:`int`
            The id number of the second atom.

        Returns
        -------
        :class:`rdkit.Chem.rdchem.BondType.SINGLE`
            If the atoms don't belong to functional groups which form
            a double bond.

        :class:`rdkit.Chem.rdchem.BondType.DOUBLE`
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
        atom2_grp = atom2.GetProp('fg')

        double_bond_present = ((atom1_grp, atom2_grp) == tup or
                               (atom2_grp, atom1_grp) == tup for
                               tup in double_bond_combs)

        if any(double_bond_present):
            return rdkit.rdchem.BondType.DOUBLE
        else:
            return rdkit.rdchem.BondType.SINGLE

    def __str__(self):
        return repr(self)

    def __repr__(self):
        # The `_repr` attribute is made in the TopologyMeta __call__()
        # method, when the Topology object is instantiated.
        return self._repr

    def __eq__(self, other):
        return repr(self) == repr(other)

    def __hash__(self):
        return id(self)


class Linear(Topology):
    """
    A class represting linear polymers.

    Attributes
    ----------
    repeating_unit : :class:`str`
        A string showing the repeating unit of the :class:`.Polymer`.
        For example, ``"AB"`` or ``"ABB"``. The building block with
        index ``0`` in :attr:`.MacroMolecule.building_blocks` is
        labelled as ``"A"`` while index ``1`` as ``"B"`` and so on.

    orientation : :class:`tuple` of :class:`int`
        For each character in the repeating unit, a value of ``-1``,
        ``0`` or ``1`` must be given in a :class:`list`. It indicates
        the direction at which each monomer of the repeating unit is
        placed along the chain. ``0`` means that the direction is
        random.

    n : :class:`int`
        The number of repeating units which are used to make the
        polymer.

    ends : :class:`str`
        The string represents how the end groups of the polymer are
        treated. If ``'h'`` the functional groups at the end of the
        polymer are converted into hydrogem atoms. If ``'fg'`` they are
        kept as the original functional group.

    """

    def __init__(self, repeating_unit, orientation, n, ends='h'):
        """
        Initializes a :class:`Linear` instance.

        Parameters
        ----------
        repeating_unit : :class:`str`
            A string showing the repeating unit of the
            :class:`.Polymer`. For example, ``"AB"`` or ``"ABB"``. The
            building block with index ``0`` in
            :attr:`.MacroMolecule.building_blocks` is labelled as
            ``"A"`` while index ``1`` as ``"B"`` and so on.

        orientation : :class:`tuple` of :class:`float`
            For each character in the repeating unit, a value between
            ``0`` (inclusive) and ``1`` (inclusive) must be given.
            The values give the probability that each monomer is
            flipped by 180 degrees when being added to the chain. If
            ``0`` then the monomer is guaranteed to not flip. If ``1``
            it is guaranteed to flip. This allows the user to create
            head-to-head or head-to-tail chains, as well as chain with
            a preference for head-to-head or head-to-tail if a number
            between ``0`` and ``1`` is chosen.

        n : :class:`int`
            The number of repeating units which are used to make the
            polymer.

        ends : :class:`str`, optional
            The string represents how the end groups of the polymer are
            treated. If ``'h'`` the functional groups at the end of the
            polymer are converted into hydrogem atoms. If ``'fg'`` they
            are kept as the original functional group.

        """

        self.repeating_unit = repeating_unit
        self.orientation = tuple(orientation)
        self.n = n
        self.ends = ends

    def del_atoms(self, macro_mol):
        """
        Deletes the atoms which are lost during assembly.

        Parameters
        ----------
        macro_mol : :class:`.Polymer`
            The polymer being assembled.

        Returns
        -------
        None : :class:`NoneType`

        """

        if self.ends == 'h':
            self.hygrogen_ends(macro_mol)
        elif self.ends == 'fg':
            self.fg_ends(macro_mol)

    def fg_ends(self, macro_mol):
        """
        Removes almost all atoms tagged for deletion.

        In polymers, you don't want to delete the atoms at the ends of
        the chain.

        Parameters
        ----------
        macro_mol : :class:`.Polymer`
            The polymer being assembled.

        Returns
        -------
        None : :class:`NoneType`

        """

        fgs = set()
        # Get all atoms tagged for deletion, held in tuples
        # corresponding to individual functional groups.
        for bb in macro_mol.building_blocks:
            delmol = rdkit.MolFromSmarts(bb.func_grp.del_smarts)
            fgs = fgs.union(macro_mol.mol.GetSubstructMatches(delmol))

        # Get the functional groups which hold the atoms with the
        # smallest and largest values for the x coordinate need to have
        # deletion tags removed.
        maxid = max(flatten(fgs),
                    key=lambda x: macro_mol.atom_coords(x)[0])
        maxfg = next(fg for fg in fgs if maxid in fg)

        minid = min(flatten(fgs),
                    key=lambda x: macro_mol.atom_coords(x)[0])
        minfg = next(fg for fg in fgs if minid in fg)

        for atom_id in chain(flatten(minfg), flatten(maxfg)):
            atom = macro_mol.mol.GetAtomWithIdx(atom_id)
            atom.ClearProp('del')

        super().del_atoms(macro_mol)

    def hygrogen_ends(self, macro_mol):
        """
        Removes all atoms tagged for deletion and adds hydrogens.

        In polymers, you want to replace the functional groups at the
        ends with hydrogen atoms.

        Parameters
        ----------
        macro_mol : :class:`.Polymer`
            The polymer being assembled.

        Returns
        -------
        None : :class:`NoneType`

        """

        # Remove all extra atoms.
        super().del_atoms(macro_mol)
        # Add hydrogens.
        for atom in macro_mol.mol.GetAtoms():
            atom.UpdatePropertyCache()
        macro_mol.mol = rdkit.AddHs(macro_mol.mol, addCoords=True)

    def place_mols(self, macro_mol):
        """
        Places monomers side by side and joins them.

        The monomers are placed along the x-axis, so that the vector
        running between the functional groups is placed on the axis.

        Parameters
        ----------
        macro_mol : :class:`.Polymer`
            The polymer being assembled.

        Returns
        -------
        None : :class:`NoneType`

        """

        # Make a map from monomer label to object.
        mapping = {}
        # Assign every monomer a label ("A", "B", "C", etc.).
        for label, monomer in zip(dedupe(self.repeating_unit),
                                  macro_mol.building_blocks):
            mapping[label] = monomer

        # Make a que for holding bonder atom ids.
        self.bonders = deque(maxlen=2)
        # Make string representing the entire polymer, not just the
        # repeating unit.
        polymer = self.repeating_unit*self.n
        # Get the direction for each monomer along the entire chain,
        # not just the repeating unit.
        dirs = self.orientation*self.n

        # Go through the repeating unit. Place each monomer 50 A apart.
        # Also create a bond.
        macro_mol.mol = rdkit.Mol()
        for i, (label, mdir) in enumerate(zip(polymer, dirs)):
            self.bonders.append([
                macro_mol.mol.GetNumAtoms() + id_ for
                id_ in mapping[label].bonder_ids])

            # Flip or not flip the monomer as given by the probability
            # in `mdir`.
            mdir = np.random.choice([1, -1], p=[mdir, 1-mdir])
            mapping[label].set_orientation2([mdir, 0, 0])

            # The first building block should be placed at 0, the others
            # have positions calculated based on bb size.
            x_coord = self._x_position(macro_mol, mapping[label]) if i else 0
            monomer_mol = mapping[label].set_position([x_coord, 0, 0])

            bb_index = macro_mol.building_blocks.index(mapping[label])
            add_fragment_props(monomer_mol, bb_index, i)

            macro_mol.mol = rdkit.CombineMols(macro_mol.mol,
                                              monomer_mol)
            if i != 0:
                self.join(macro_mol)

    def join(self, macro_mol):
        """
        Joins 2 monomers.

        Parameters
        ----------
        macro_mol : :class:`.Polymer`
            The polymer being assembled.

        Returns
        -------
        None : :class:`NoneType`

        """

        distances = []
        for bonder1 in self.bonders[0]:
            for bonder2 in self.bonders[1]:
                distances.append(
                        (macro_mol.atom_distance(bonder1, bonder2),
                         bonder1,
                         bonder2))

        _, bonder1, bonder2 = min(distances)
        emol = rdkit.EditableMol(macro_mol.mol)
        emol.AddBond(bonder1, bonder2, self.determine_bond_type(
                                                            macro_mol,
                                                            bonder1,
                                                            bonder2))

        macro_mol.mol = emol.GetMol()

    def join_mols(self, macro_mol):
        """
        Does nothing, :meth:`place_mols` joins up the molecules too.

        Parameters
        ----------
        macro_mol : :class:`.Polymer`
            The polymer being assembled.

        Returns
        -------
        None : :class:`NoneType`

        """

        return

    def _x_position(self, macro_mol, bb):
        """
        Calculates the x coordinate on which to place `bb`.

        Does this by checking the most how for down the x axis
        `macro_mol` stretches and checking the distance between
        the minimum x position of `bb` and its centroid.
        It then tries to place `bb` about 3 A away from `macro_mol`.

        Parameters
        ----------
        macro_mol : :class:`.MacroMolecule`
            The macromolecule being assembled.

        bb : :class:`.StructUnit`
            The building block to be added to `macro_mol`.

        Returns
        -------
        :class:`float`
            The x coordinate on which to place `bb`.

        """

        mm_max_x = max(macro_mol.all_atom_coords(),
                       key=lambda x: x[1][0])[1][0]
        bb_min_x = min(bb.all_atom_coords(),
                       key=lambda x: x[1][0])[1][0]
        bb_len = bb.centroid()[0] - bb_min_x
        return mm_max_x + bb_len + 3
