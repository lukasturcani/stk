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
:meth:`.Reactor.react`. The latter performs any final cleanup
operations on the constructed molecule. For example, converting the end
functional groups of a polymer into hydrogen atoms. See also
:meth:`Topology.cleanup`.

During the construction process, every time a building block is placed
in the :class:`.ConstructedMolecule`, new :class:`FunctionalGroup`
instances must be made, which correspond to the functional groups added
by virtue of adding the building block. These must be added to the
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
from inspect import signature
from collections import Counter

from ..functional_groups import Reactor


def remove_confs(building_blocks, keep):
    """
    Removes all conformers from `building_blocks` except `keep`.

    All kept conformers have their id set to ``0``.

    Parameters
    ----------
    building_blocks : :class:`iterable` of :class:`.StructUnit`
        A :class:`.set` of :class:`.StructUnit` instances which
        represent the building blocks forming a
        :class:`.ConstructedMolecule`.

    keep : :class:`list` of :class:`int`
        The ids of the building block conformers to be used for
        constructing the :class:`.ConstructedMolecule`. Must be equal
        in length to `building_blocks` and the orders must correspond.

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

    keep_ids = [
        bb.mol.GetConformer(id_).GetId()
        for bb, id_ in zip(building_blocks, keep)
    ]

    original_confs = [
        [rdkit.Conformer(conf) for conf in bb.mol.GetConformers()]
        for bb in building_blocks
    ]
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
    Constructs :class:`.ConstructedMolecule` from building blocks.

    More accurately, child classes of :class:`Topology` take care of
    constructing :class:`.ConstructedMolecule` from building blocks.

    This class directly defines any operations and attributes that are
    needed by any topology during construction. However, this class is
    not used directly. It is intended to be inherited. All
    :attr:`.ConstructedMolecule.topology` attributes hold an instance
    of a :class:`Topology` child class. Child classes of
    :class:`Topology` define operations specific to that one topology.
    For example, each child class must define a :meth:`join_mols`,
    which creates bonds between the building blocks of a
    :class:`.ConstructedMolecule`. The way in which this is done will
    depend on what kind of molecules are being constructed. In
    addition, each child class must define methods which place the
    building blocks in appropriate positions.

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

    def construct(self, mol, bb_conformers=None):
        """
        Constructs :mod:`rdkit` instances of molecules.

        This method places an :mod:`rdkit` molecule of the
        :class:`.ConstructedMolecule`
        into the :attr:`~.Molecule.mol` attribute of
        :class:`.ConstructedMolecule`. It also creates the
        :attr:`.ConstructedMolecule.bb_counter` and
        :attr:`.ConstructedMolecule.bonds_made` attributes.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The :class:`.ConstructedMolecule` instance which needs to
            be constructed.

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
                -1 for _ in range(len(mol.building_blocks))
            ]

        # During construction, only a single conformer should exist per
        # building block. Otherwise, rdkit.CombineMols won't work. It
        # only combines conformers with the same id.
        original_confs = remove_confs(
            mol.building_blocks,
            bb_conformers
        )

        mol.bonds_made = 0
        mol.mol = rdkit.Mol()
        mol.bb_counter = Counter()

        self.reactor = Reactor()
        self.place_mols(mol)
        self.prepare(mol)

        self.reactor.set_molecule(mol.mol)
        mol.func_groups = self.reactor.func_groups

        for fgs in self.bonded_fgs(mol):
            self.reactor.react(*fgs, track_fgs=self.track_fgs)
        mol.mol = self.reactor.result(self.del_atoms)
        mol.bonds_made = self.reactor.bonds_made

        self.cleanup(mol)

        # Make sure that the property cache of each atom is up to date.
        for atom in mol.mol.GetAtoms():
            atom.UpdatePropertyCache()

        # Restore the original conformers.
        for bb, confs in zip(mol.building_blocks, original_confs):
            bb.mol.RemoveAllConformers()
            for conf in confs:
                bb.mol.AddConformer(conf)

        # Reactor can't be pickled because it contains an EditableMol,
        # which can't be pickled.
        self.reactor = None

    def place_mols(self, mol):
        """
        Places building blocks.

        The :mod:`rdkit` molecules of the building blocks are
        combined into a single :mod:`rdkit` molecule and placed into
        `mol.mol`.

        The function is also reponsible for updating
        :attr:`~.ConstructedMolecule.bb_counter`.

        This function must also add the tags ``'bb_index'``
        and ``'mol_index'`` to every atom in the molecule. The
        ``'bb_index'`` tag identifies which building block the atom
        belongs to. The building block is identified by its index
        within :attr:`ConstructedMolecule.building_blocks`.
        The ``'mol_index'`` identifies which molecule of a specific
        building block the atom belongs to. For example, if
        ``bb_index = 1`` and ``mol_index = 3`` the atom belongs to
        the 4th molecule of ``mol.building_blocks[1]`` to
        be added to the :class:`.ConstructedMolecule`. The utility
        function :func:`.add_fragment_props` is provided to help with
        this.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The molecule being constructed.

        Raises
        ------
        :class:`NotImplementedError`

        """

        raise NotImplementedError()

    def bonded_fgs(self, mol):
        """
        An iterator which yields functional groups to be bonded.

        This iterator must yield :class:`tuple`s of
        :class:`.FunctionalGroup` molecules. These are the functional
        groups in the molecule being constructed which need to be
        bonded.

        This :class:`tuple` gets passed to :meth:`Reactor.react`.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The molecule being constructed.

        Raises
        ------
        :class:`NotImplementedError`

        """

        raise NotImplementedError()

    def cleanup(self, mol):
        """
        Performs final clean up actions on a constructed molecule.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The molecule being constructed.

        Returns
        -------
        None : :class:`NoneType`

        """

        return

    def prepare(self, mol):
        """
        Performs ops between placing and reacting building blocks.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The molecule being constructed.

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
