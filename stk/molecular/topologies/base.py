"""
Defines the base :class:`Topology` type.

.. _`adding topologies`:

Extending stk: Adding new topologies.
-------------------------------------

General
.......

A new topology class must be defined. The class must inherit
:class:`Topology`. The new topology class will have to define
the methods :meth:`place_mols` and :meth:`bonded_fgs`. A description
of what these methods should do is given by :meth:`Topology.place_mols`
and :meth:`Topology.bonded_fgs`.

The new class may optionally define the methods :meth:`prepare` and
:meth:`cleanup`. The former perfroms operations on the molecule
before it is joined up and has atoms deleted via :func:`.react`.
The latter any final cleanup operations on the assembled molecule. For
example, converting the end functional groups of a polymer into
hydrogen atoms. See also :meth:`Topology.cleanup`.


Cages
.....

To add a new cage topology a new class should be created, named
after the topology. This class should inherit :class:`.CageTopology`.
This will give access to various methods which are necessary
for dealing with any cage molecule. See the documenation of
:class:`.CageTopology` for more details.

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
import numpy as np
from inspect import signature
from collections import Counter

from ..functional_groups import react
from ...utilities import dedupe, add_fragment_props, remake


def remove_confs(building_blocks, keep):
    """
    Removes all conformers from `building_blocks` except `keep`.

    All kept conformers have their id set to ``0``.

    Parameters
    ----------
    building_blocks : :class:`iterable` of :class:`.StructUnit`
        A set of :class:`.StructUnit` instances which represent the
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

    Attributes
    ----------
    react_del : :class:`bool`
        Toggles whether atoms with the ``'del'`` propety are deleted
        by :func:`.react`.

    """

    def __init__(self, react_del=True):
        self.react_del = react_del

    def build(self, macro_mol, bb_conformers=None):
        """
        Assembles ``rdkit`` instances of macromolecules.

        This method places an ``rdkit`` molecule of the assembled
        macromolecule into the :attr:`~.Molecule.mol` attribute of
        :class:`.MacroMolecule`. It also creates the
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

        macro_mol.bonds_made = 0
        macro_mol.mol = rdkit.Mol()
        macro_mol.bb_counter = Counter()

        self.place_mols(macro_mol)
        self.prepare(macro_mol)
        for fgs in self.bonded_fgs(macro_mol):
            macro_mol.mol, new_bonds = react(macro_mol.mol,
                                             self.react_del,
                                             *fgs)
            macro_mol.bonds_made += new_bonds
        self.cleanup(macro_mol)

        # Make sure that the property cache of each atom is up to date.
        for atom in macro_mol.mol.GetAtoms():
            atom.UpdatePropertyCache()

        # Restore the original conformers.
        for bb, confs in zip(macro_mol.building_blocks, original_confs):
            bb.mol.RemoveAllConformers()
            for conf in confs:
                bb.mol.AddConformer(conf)

    def place_mols(self, macro_mol):
        """
        Places building blocks.

        The ``rdkit`` molecules of the building blocks are
        combined into a single ``rdkit`` molecule and placed into
        `macro_mol.mol`. Beyond this, the function must give every
        functional group in the macromolecule a unique id. This is
        done by adding the property ``'fg_id'`` to every atom part
        of a functional group in `macro_mol`. Atoms in the same
        functional group will have the same value of ``'fg_id'`` while
        atoms in different functional groups will have a different
        ``'fg_id'``.

        The function is also reponsible for updating
        :attr:`~.MacroMolecule.bb_counter`.

        Finally, this function must also add the tags ``'bb_index'``
        and ``'mol_index'`` to every atom in the molecule. The
        ``'bb_index'`` tag identifies which building block the atom
        belongs to. The building block is identified by its index
        within :attr:`MacroMolecule.building_blocks`.
        The ``'mol_index'`` identifies which molecule of a specific
        building the atom belongs to. For example, if
        ``bb_index = 1`` and ``mol_index = 3`` the atom belongs to
        the 4th molecule of ``macro_mol.building_blocks[1]`` to
        be added to the macromolecule. The utility function
        :func:`.add_fragment_props` is provided to help with this.


        Parameters
        ----------
        macro_mol : :class:`MacroMolecule`
            The molecule being assembled.

        Raises
        ------
        :class:`NotImplementedError`

        """

        raise NotImplementedError()

    def bonded_fgs(self, macro_mol):
        """
        An iterator which yields ids of functional groups to be bonded.

        The ids of functional groups in `macro_mol` are assigned by
        :meth:`place_mols`. The ids are stored as atom properties
        under the name ``'fg_id'``, This method looks at the
        macromolecule and determines which functional groups should be
        bonded to create the final macromolecule. It then yields the
        ids functional groups as a :class:`tuple`.

        This :class:`tuple` gets passed to :func:`.react`.

        Parameters
        ----------
        macro_mol : :class:`MacroMolecule`
            The molecule being assembled.

        Raises
        ------
        :class:`NotImplementedError`

        """

        raise NotImplementedError()

    def cleanup(self, macro_mol):
        """
        Performs final clean up actions on an assembled molecule.

        Parameters
        ----------
        macro_mol : :class:`.MacroMolecule`
            The molecule being assembled.

        Returns
        -------
        None : :class:`NoneType`

        """

        return

    def prepare(self, macro_mol):
        """
        Performs ops between placing and reacting building blocks.

        Parameters
        ----------
        macro_mol : :class:`.MacroMolecule`
            The molecule being assembled.

        Returns
        -------
        None : :class:`NoneType`

        """

        return

    def update_fg_id(self, macro_mol, mol):
        """

        """

        mol = rdkit.Mol(mol)
        max_id = max((a.GetIntProp('fg_id') for
                      a in macro_mol.mol.GetAtoms() if
                      a.HasProp('fg_id')), default=-1) + 1
        for a in mol.GetAtoms():
            if a.HasProp('fg_id'):
                a.SetIntProp('fg_id', a.GetIntProp('fg_id')+max_id)
        return mol

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

    orientation : :class:`tuple` of :class:`float`
        For each character in the repeating unit, a value between ``0``
        and ``1`` (both inclusive) must be given in a :class:`list`. It
        indicates the probability that each monomer will have its
        orientation along the chain flipped.

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
        super().__init__()

    def cleanup(self, macro_mol):
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

        emol = rdkit.EditableMol(macro_mol.mol)
        # Remove all extra atoms.
        for atom in reversed(macro_mol.mol.GetAtoms()):
            if atom.HasProp('del'):
                emol.RemoveAtom(atom.GetIdx())

        macro_mol.mol = remake(emol.GetMol())
        macro_mol.mol = rdkit.AddHs(macro_mol.mol, addCoords=True)

    def place_mols(self, macro_mol):
        """
        Places monomers side by side.

        The monomers are placed along the x-axis, so that the vector
        running between the functional groups is placed along the axis.
        Functional groups are tagged with ``'fg_id'`` such that
        ``'fg_id'`` increases along the x-axis.

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

        # Make string representing the entire polymer, not just the
        # repeating unit.
        polymer = self.repeating_unit*self.n
        # Get the direction for each monomer along the entire chain,
        # not just the repeating unit.
        dirs = self.orientation*self.n

        # Go through the repeating unit and place each monomer.
        bb_counter = macro_mol.bb_counter = {}
        for i, (label, mdir) in enumerate(zip(polymer, dirs)):
            bb = mapping[label]
            bb_counter[bb] = bb_counter.get(bb, 0) + 1
            original_position = bb.position_matrix()

            # Flip or not flip the monomer as given by the probability
            # in `mdir`.
            mdir = np.random.choice([1, -1], p=[mdir, 1-mdir])
            bb.set_orientation2([mdir, 0, 0])

            # The first building block should be placed at 0, the
            # others have positions calculated based on bb size.
            x_coord = self._x_position(macro_mol, bb) if i else 0
            monomer_mol = bb.set_bonder_centroid([x_coord, 0, 0])
            monomer_mol = rdkit.Mol(monomer_mol)

            bb_index = macro_mol.building_blocks.index(bb)
            add_fragment_props(monomer_mol, bb_index, i)

            # Add fg_id tags.

            # Check which funcitonal group is at the back and which
            # one at the front.
            if len(bb.bonder_ids) == 1:
                c1, c2 = bb.bonder_centroid(), bb.centroid()
            else:
                c1, c2 = list(bb.bonder_centroids())
            front = 1 if c1[0] < c2[0] else 0
            back = 1 if front != 1 else 0

            for atom in monomer_mol.GetAtoms():
                if atom.HasProp('fg'):
                    fg_id = atom.GetIntProp('fg_id')
                    fg_id = 2*i if fg_id == back else 2*i+1
                    atom.SetIntProp('fg_id', fg_id)

            macro_mol.mol = rdkit.CombineMols(macro_mol.mol,
                                              monomer_mol)

            bb.set_position_from_matrix(original_position)

    def bonded_fgs(self, macro_mol):
        """
        Yields functional groups to react.

        Parameters
        ----------
        macro_mol : :class:`.Polymer`
            The polymer being assembled.

        Yields
        -------
        :class:`tuple` of :class:`int`
            Holds the ids of the functional groups set to react.

        """

        for i in range(1, 2*len(self.repeating_unit)*self.n-1, 2):
            yield i, i+1

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
