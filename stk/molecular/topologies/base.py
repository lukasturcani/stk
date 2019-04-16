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
:meth:`cleanup`. The former performs operations on the molecule
before it is joined up and has atoms deleted via
:meth:`.Reactor.react`. The latter any final cleanup operations on the
assembled molecule. For example, converting the end functional groups
of a polymer into hydrogen atoms. See also :meth:`Topology.cleanup`.

During the build process, every time a building block is placed
in the macromolecule, new :class:`FunctionalGroup` instances must be
made, which correspond to the functional groups added by virtue of
adding the building block. These must be added to the
:class:`Reactor` held in :attr:`.Topology.reactor`, specifically into
its :attr:`.Reactor.func_groups` attribute. This means that the
reactor will keep the atom ids in these functional groups up to date
when it deletes atoms. However, note that any functional groups
yielded by :meth:`.Topology.bonded_fgs` are automatically added, so
they do not have to be managed manually. If you do not wish to
automatically add the functional groups into
:attr:`.Reactor.func_groups` you can toggle it in
:attr:`Topology.track_fgs`.

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

from ..functional_groups import Reactor
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
        c = ', '.join(
            f'{key!s}={value!r}' for key, value in sorted(sig.items())
        )
        obj._repr = f'{self.__name__}({c})'
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
    del_atoms : :class:`bool`
        Toggles whether deleter atoms are deleted by
        :meth:`.Reactor.result`.

    reactor : :class:`.Reactor`
        The reactor which performs the reactions.

    track_fgs : :class:`bool`
        Toggles whether functional groups yielded by
        :meth:`bonded_fgs` are automatically added into
        :attr:`.Reactor.func_groups`.

    """

    def __init__(self, del_atoms=True, track_fgs=True):
        self.del_atoms = del_atoms
        self.track_fgs = track_fgs

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
            bb_conformers = [
                -1 for _ in range(len(macro_mol.building_blocks))
            ]

        # When building, only a single conformer should exist per
        # building block. Otherwise, rdkit.CombineMols won't work. It
        # only combines conformers with the same id.
        original_confs = remove_confs(macro_mol.building_blocks,
                                      bb_conformers)

        macro_mol.bonds_made = 0
        macro_mol.mol = rdkit.Mol()
        macro_mol.bb_counter = Counter()

        self.reactor = Reactor()
        self.place_mols(macro_mol)
        self.prepare(macro_mol)

        self.reactor.set_molecule(macro_mol.mol)
        macro_mol.func_groups = self.reactor.func_groups

        for fgs in self.bonded_fgs(macro_mol):
            self.reactor.react(*fgs, track_fgs=self.track_fgs)
        macro_mol.mol = self.reactor.result(self.del_atoms)
        macro_mol.bonds_made = self.reactor.bonds_made

        self.cleanup(macro_mol)

        # Make sure that the property cache of each atom is up to date.
        for atom in macro_mol.mol.GetAtoms():
            atom.UpdatePropertyCache()

        # Restore the original conformers.
        for bb, confs in zip(macro_mol.building_blocks,
                             original_confs):
            bb.mol.RemoveAllConformers()
            for conf in confs:
                bb.mol.AddConformer(conf)

        # Reactor can't be pickled because it contains an EditableMol,
        # which can't be pickled.
        self.reactor = None

    def place_mols(self, macro_mol):
        """
        Places building blocks.

        The ``rdkit`` molecules of the building blocks are
        combined into a single ``rdkit`` molecule and placed into
        `macro_mol.mol`.

        The function is also reponsible for updating
        :attr:`~.MacroMolecule.bb_counter`.

        This function must also add the tags ``'bb_index'``
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
        macro_mol : :class:`.MacroMolecule`
            The molecule being assembled.

        Raises
        ------
        :class:`NotImplementedError`

        """

        raise NotImplementedError()

    def bonded_fgs(self, macro_mol):
        """
        An iterator which yields functional groups to be bonded.

        This iterator must yield :class:`tuple`s of
        :class:`.FunctionalGroup` molecules. These are the functional
        groups in the macromolecule which need to be bonded.

        This :class:`tuple` gets passed to :meth:`Reactor.react`.

        Parameters
        ----------
        macro_mol : :class:`.MacroMolecule`
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
    A class representing linear polymers.

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

    def __init__(self, repeating_unit, orientation, n, ends='fg'):
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
        super().__init__(track_fgs=False)

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
        Removes all deleter atoms and adds hydrogens.

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

        deleters = []
        for func_group in self.reactor.func_groups:
            deleters.extend(func_group.deleter_ids)

        emol = rdkit.EditableMol(macro_mol.mol)
        for atom_id in sorted(deleters, reverse=True):
            emol.RemoveAtom(atom_id)
        macro_mol.mol = remake(emol.GetMol())
        macro_mol.mol = rdkit.AddHs(macro_mol.mol, addCoords=True)

    def place_mols(self, macro_mol):
        """
        Places monomers side by side.

        The monomers are placed along the x-axis, so that the vector
        running between the functional groups is placed along the axis.

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
        for i, (label, mdir) in enumerate(zip(polymer, dirs)):
            bb = mapping[label]
            macro_mol.bb_counter.update([bb])
            original_position = bb.mol.GetConformer().GetPositions().T

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

            # Check which functional group is at the back and which
            # one at the front.
            n_fgs = len(bb.func_groups)
            if n_fgs == 2:
                c1, c2 = list(bb.bonder_centroids())
                front = 1 if c1[0] < c2[0] else 0

            # Flag to see if the first functional group in
            # self.reactor.func_groups should be bonded.
            if i == 0:
                self.bond_first = True if n_fgs == 1 else False

            num_atoms = macro_mol.mol.GetNumAtoms()
            macro_mol.mol = rdkit.CombineMols(macro_mol.mol,
                                              monomer_mol)

            for fg in bb.func_groups:
                if n_fgs == 2:
                    id_ = 2*i + 1 if fg.id == front else 2*i
                elif len(self.reactor.func_groups) == 0:
                    id_ = 0
                else:
                    id_ = len(self.reactor.func_groups)

                func_group = fg.shifted_fg(id_, num_atoms)
                self.reactor.func_groups.append(func_group)

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

        fgs = sorted(self.reactor.func_groups, key=lambda fg: fg.id)

        start = 0 if self.bond_first else 1
        for i in range(start, len(self.reactor.func_groups)-1, 2):
            yield fgs[i], fgs[i+1]

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
