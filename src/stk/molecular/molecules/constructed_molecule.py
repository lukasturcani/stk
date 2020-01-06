"""
Constructed Molecule
====================

"""

import logging
import numpy as np
from collections import Counter

from .. import elements, bonds, topology_graphs
from .molecule import Molecule
from ..functional_groups import FunctionalGroup, fg_types

logger = logging.getLogger(__name__)


class ConstructionError(Exception):
    ...


class ConstructedMolecule(Molecule):
    """
    Represents constructed molecules.

    A :class:`ConstructedMolecule` requires at least 2 basic pieces of
    information: which building block molecules are used to construct
    the molecule and what the :class:`.TopologyGraph` of the
    constructed molecule is. The construction of the molecular
    structure is performed by :meth:`.TopologyGraph.construct`. This
    method does not have to be called explicitly by the user, it will
    be called automatically during initialization.

    The building block molecules used for construction can be either
    :class:`.BuildingBlock` instances or other
    :class:`.ConstructedMolecule` instances, or a combination both.

    Each :class:`.TopologyGraph` subclass may add additional attributes
    to the :class:`ConstructedMolecule`, which will be described within
    its documentation.

    Attributes
    ----------
    atoms : :class:`tuple` of :class:`.Atom`
        The atoms of the molecule. Each :class:`.Atom`
        instance is guaranteed to have two attributes. The
        first is :attr:`building_block`, which holds the building
        block :class:`.Molecule` from which that
        :class:`.Atom` came. If the :class:`.Atom` did not come from a
        building block, but was added by a reaction, the value
        of this attribute will be ``None``.

        The second attribute is :attr:`building_block_id`.
        This will be the same value on all atoms that came from the
        building block. Note that if a building block is used multiple
        times during construction, the :attr:`building_block_id` will
        be different for each time it is used.

    bonds : :class:`tuple` of :class:`.Bond`
        The bonds of the molecule.

    building_block_vertices : :class:`dict`
        Maps the :class:`.Molecule` instances used for construction,
        which can be either :class:`.BuildingBlock` or
        :class:`.ConstructedMolecule`, to the
        :class:`~.topologies.base.Vertex` objects they are placed on
        during construction. The :class:`dict` has the form

        .. code-block:: python

            building_block_vertices = {
                BuildingBlock(...): [Vertex(...), Vertex(...)],
                BuildingBlock(...): [
                    Vertex(...),
                    Vertex(...),
                    Vertex(...),
                ]
                ConstructedMolecule(...): [Vertex(...)]
            }

    building_block_counter : :class:`collections.Counter`
        A counter keeping track of how many times each building block
        molecule appears in the :class:`ConstructedMolecule`.

    topology_graph : :class:`.TopologyGraph`
        Defines the topology graph of :class:`ConstructedMolecule` and
        is responsible for constructing it.

    construction_bonds : :class:`tuple` of :class:`.Bond`
        Holds the bonds in :attr:`bonds`, which were added by the
        construction process.

    func_groups : :class:`tuple` of :class:`.FunctionalGroup`
        The remnants of building block functional groups present in the
        molecule. They track which atoms belonged to functional groups
        in the building block molecules. The id of each
        :class:`.FunctionalGroup` should match its index in
        :attr:`func_groups`.

    Examples
    --------
    *Initialization*

    A :class:`ConstructedMolecule` can be created from a set of
    building blocks and a :class:`.TopologyGraph`

    .. code-block:: python

        import stk

        bb1 = stk.BuildingBlock('NCCCN', ['amine'])
        bb2 = stk.BuildingBlock('O=CC(C=O)CC=O', ['aldehyde'])
        tetrahedron = stk.cage.FourPlusSix()
        cage1 = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2],
            topology_graph=tetrahedron
        )

    A :class:`ConstructedMolecule` can be used to construct other
    :class:`ConstructedMolecule` instances

    .. code-block:: python

        benzene = stk.BuildingBlock('c1ccccc1')
        cage_complex = stk.ConstructedMolecule(
            building_blocks=[cage1, benzene],
            topology_graph=stk.host_guest.Complex()
        )

    During initialization it is possible to force building blocks to
    be placed on specific :attr:`~.TopologyGraph.vertices` of the
    :class:`.TopologyGraph`

    .. code-block:: python

        bb3 = stk.BuildingBlock('NCOCN', ['amine'])
        bb4 = stk.BuildingBlock('NCOCCCOCN', ['amine'])
        cage2 = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2, bb3, bb4],
            topology_graph=tetrahedron
            building_block_vertices={
                bb1: tetrahedron.vertices[4:6]
                bb3: tetrahedron.vertices[6:7]
                bb4: tetrahedron.vertices[7:8]
            }
        )

    *Building blocks with the wrong number of functional groups.*

    If the building block has too many functional groups, you can
    remove some in order to use it

    .. code-block:: python

        chain = stk.polymer.Linear('AB', [0, 0],  3)

        # This won't work, bb2 has 3 functional groups but 2 are needed
        # for monomers in a linear polymer chain.
        failed = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2],
            topology_graph=chain
        )

        # Remove one of the functional groups and you will be able to
        # construct the chain.
        bb2.func_groups = (bb2.func_groups[0], bb2.func_groups[2])
        failed = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2],
            topology_graph=chain
        )

    """

    @classmethod
    def _construct(
        cls,
        building_blocks,
        topology_graph,
        building_block_vertices=None,
        use_cache=False
    ):
        if building_block_vertices is None:
            building_block_vertices = (
                topology_graph.assign_building_blocks_to_vertices(
                    building_blocks=building_blocks
                )
            )
        identity_key = cls._get_identity_key_from_components(
            building_blocks=building_blocks,
            topology_graph=topology_graph,
            building_block_vertices=building_block_vertices
        )
        if use_cache and identity_key in cls._cache:
            return cls._cache[identity_key]

        obj = cls._init_from_components(
            building_blocks=building_blocks,
            topology_graph=topology_graph,
            building_block_vertices=building_block_vertices,
            identity_key=identity_key
        )
        if use_cache:
            cls._cache[identity_key] = obj
        return obj

    def __init__(
        self,
        building_blocks,
        topology_graph,
        building_block_vertices=None,
        use_cache=False
    ):
        """
        Initialize a :class:`ConstructedMolecule`.

        Parameters
        ----------
        building_blocks : :class:`list` of :class:`.Molecule`
            The :class:`.BuildingBlock` and
            :class:`ConstructedMolecule` instances which
            represent the building block molecules used for
            construction. Only one instance is present per building
            block molecule, even if multiples of that building block
            join up to form the :class:`ConstructedMolecule`.

        topology_graph : :class:`.TopologyGraph`
            Defines the topology graph of the
            :class:`ConstructedMolecule` and constructs it.

        building_block_vertices : :class:`dict`, optional
            Maps the :class:`.Molecule` in  `building_blocks` to the
            :class:`~.topologies.base.Vertex` instances in
            `topology_graph` it is placed on. Each
            :class:`.BuildingBlock` and :class:`ConstructedMolecule`
            can be mapped to multiple :class:`~.topologies.base.Vertex`
            objects. See the examples section in the
            :class:`.ConstructedMolecule` class docstring to help
            understand how this parameter is used. If ``None``,
            building block molecules will be assigned to vertices at
            random.

        use_cache : :class:`bool`, optional
            If ``True``, a new :class:`.ConstructedMolecule` will
            not be made if a cached and identical one already exists,
            the one which already exists will be returned. If ``True``
            and a cached, identical :class:`ConstructedMolecule` does
            not yet exist the created one will be added to the cache.

        """

        # This method does not get called, See _construct().
        raise RuntimeError('This method should not be getting called.')

    @classmethod
    def _init_from_components(
        cls,
        building_blocks,
        topology_graph,
        building_block_vertices,
        identity_key
    ):
        """
        Initialize a :class:`ConstructedMolecule`.

        Parameters
        ----------
        building_blocks : :class:`list` of :class:`.Molecule`
            The :class:`.BuildingBlock` and
            :class:`ConstructedMolecule` instances which
            represent the building block molecules used for
            construction. Only one instance is present per building
            block molecule, even if multiples of that building block
            join up to form the :class:`ConstructedMolecule`.

        topology_graph : :class:`.TopologyGraph`
            Defines the topology graph of the
            :class:`ConstructedMolecule` and constructs it.

        building_block_vertices : :class:`dict`
            Maps the :class:`.Molecule` in  `building_blocks` to the
            :class:`~.topologies.base.Vertex` in `topology_graph`.
            Each :class:`.BuildingBlock` and
            :class:`ConstructedMolecule` can be mapped to multiple
            :class:`~.topologies.base.Vertex` objects. See the
            examples section in the :class:`.ConstructedMolecule`
            class docstring to help understand how this parameter
            is used. If ``None``, building block molecules will be
            assigned to vertices at random.

        identity_key : :class:`tuple`
            The identity key of the molecule. The identity key wil be
            equal for two molecules which ``stk`` sees as identical.

        Returns
        -------
        :class:`.ConstructedMolecule`
            The molecule.

        """

        obj = cls.__new__(cls)
        obj._identity_key = identity_key
        obj.building_block_vertices = building_block_vertices
        obj.topology_graph = topology_graph
        obj.atoms = []
        obj.bonds = []
        obj.construction_bonds = []
        obj.func_groups = []
        obj.building_block_counter = Counter()
        # A (3, n) numpy.ndarray holding the position of every atom in
        # the molecule.
        obj._position_matrix = []

        try:
            topology_graph.construct(obj)

        except Exception as ex:
            errormsg = (
                'Construction failure.\n'
                '\n'
                'topology_graph\n'
                '--------------\n'
                f'{topology_graph}\n'
                '\n'
                'building blocks\n'
                '---------------\n'
            )

            bb_blocks = []
            for i, bb in enumerate(building_blocks):
                bb_blocks.append(
                    f'{bb}\n\n'
                    'MDL MOL BLOCK\n'
                    '-------------\n'
                    f'{bb._to_mdl_mol_block()}'
                )

            errormsg += '\n'.join(bb_blocks)
            raise ConstructionError(errormsg) from ex

        obj.atoms = tuple(obj.atoms)
        obj.bonds = tuple(obj.bonds)
        obj.construction_bonds = tuple(obj.construction_bonds)
        obj.func_groups = tuple(obj.func_groups)

        # Ensure that functional group ids are set correctly.
        for id_, func_group in enumerate(obj.func_groups):
            func_group.id = id_

        return obj

    def clone(self):
        """
        Return a clone.

        Returns
        -------
        :class:`.ConstructedMolecule`
            The clone.

        """

        clone = super().clone()
        clone.building_block_vertices = dict(
            self.building_block_vertices
        )
        clone.building_block_counter = Counter(
            self.building_block_counter
        )
        clone.topology_graph = self.topology_graph
        construction_bonds = set(self.construction_bonds)
        clone.construction_bonds = tuple(
            clone.bonds[i] for i, bond in enumerate(self.bonds)
            if bond in construction_bonds
        )
        atom_map = {
            original: clone
            for original, clone in zip(self.atoms, clone.atoms)
        }
        clone.func_groups = tuple(
            fg.clone(atom_map) for fg in self.func_groups
        )
        return clone

    def get_building_blocks(self):
        """
        Yield the building blocks.

        Yields
        ------
        :class:`.Molecule`
            A building block of the :class:`ConstructedMolecule`.

        """

        yield from self.building_block_vertices.keys()

    def to_dict(self, include_attrs=None, ignore_missing_attrs=False):
        """
        Return a :class:`dict` representation.

        Parameters
        ----------
        include_attrs : :class:`list` of :class:`str`, optional
            The names of additional attributes of the molecule to be
            added to the :class:`dict`. Each attribute is saved as a
            string using :func:`repr`. These attributes are also
            passed down recursively to the building block molecules.

        ignore_missing_attrs : :class:`bool`, optional
            If ``False`` and an attribute in `include_attrs` is not
            held by the :class:`ConstructedMolecule`, an error will be
            raised.

        Returns
        -------
        :class:`dict`
            A :class:`dict` which represents the molecule.

        """

        if include_attrs is None:
            include_attrs = []

        bb_counter = []
        building_blocks = []
        building_block_vertices = {}
        bb_index = {}
        for i, bb in enumerate(self.building_block_vertices):
            bb_index[bb] = i
            bb_counter.append(self.building_block_counter[bb])
            building_blocks.append(bb.to_dict(include_attrs, True))
            building_block_vertices[i] = [
                self.topology_graph.vertices.index(v)
                for v in self.building_block_vertices[bb]
            ]

        # For atoms, bond sand func groups, need to change all atom
        # references with the id of the atom. This is because when
        # these objects a reconstructed you want to create pointers
        # back to the atoms found in the atoms attribute.
        atoms = []
        for atom in self.atoms:
            clone = atom.clone()
            # This attribute is a pointer that will need to get
            # re-created by load().
            clone.building_block = bb_index[clone.building_block]
            atoms.append(clone)

        bonds = []
        bond_indices = {}
        for i, bond in enumerate(self.bonds):
            clone = bond.clone()
            clone.atom1 = clone.atom1.id
            clone.atom2 = clone.atom2.id
            bonds.append(clone)
            bond_indices[bond] = i

        construction_bonds = [
            bond_indices[bond] for bond in self.construction_bonds
        ]

        func_groups = []
        for fg in self.func_groups:
            clone = fg.clone()
            clone.atoms = tuple(clone.get_atom_ids())
            clone.bonders = tuple(clone.get_bonder_ids())
            clone.deleters = tuple(clone.get_deleter_ids())
            clone.fg_type = f'{clone.fg_type.name!r}'
            func_groups.append(clone)

        d = {
            'building_blocks': building_blocks,
            'building_block_vertices': building_block_vertices,
            'building_block_counter': bb_counter,
            'construction_bonds': construction_bonds,
            'class': self.__class__.__name__,
            'position_matrix': self.get_position_matrix().tolist(),
            'topology_graph': repr(self.topology_graph),
            'func_groups': repr(tuple(func_groups)),
            'atoms': repr(tuple(atoms)),
            'bonds': repr(tuple(bonds)),
            'identity_key': repr(self._identity_key),
        }

        if ignore_missing_attrs:
            d.update({
                attr: repr(getattr(self, attr))
                for attr in include_attrs
                if hasattr(self, attr)
            })
        else:
            d.update({
                attr: repr(getattr(self, attr))
                for attr in include_attrs
            })

        return d

    @classmethod
    def _init_from_dict(cls, mol_dict, use_cache):
        """
        Intialize from a :class:`dict` representation.

        Parameters
        ----------
        mol_dict : :class:`dict`
            A :class:`dict` representation of a molecule generated
            by :meth:`to_dict`.

        use_cache : :class:`bool`
            If ``True``, a new instance will not be made if a cached
            and identical one already exists, the one which already
            exists will be returned. If ``True`` and a cached,
            identical instance does not yet exist the created one will
            be added to the cache.

        Returns
        -------
        :class:`ConstructedMolecule`
            The molecule described by `mol_dict`.

        """

        d = dict(mol_dict)
        identity_key = eval(d.pop('identity_key'))
        if use_cache and identity_key in cls._cache:
            return cls._cache[identity_key]

        tops = vars(topology_graphs)
        topology_graph = eval(d.pop('topology_graph'), tops)
        bbs = [
            Molecule.init_from_dict(bb_dict, use_cache)
            for bb_dict in d.pop('building_blocks')
        ]

        obj = cls.__new__(cls)
        obj._identity_key = identity_key
        obj.building_block_counter = Counter()
        obj.building_block_vertices = {}

        counter = d.pop('building_block_counter')
        vertices = d.pop('building_block_vertices')
        for i, bb in enumerate(bbs):
            obj.building_block_counter[bb] = counter[i]
            obj.building_block_vertices[bb] = [
                topology_graph.vertices[i] for i in vertices[i]
            ]

        obj.topology_graph = topology_graph
        obj._position_matrix = np.array(d.pop('position_matrix')).T

        obj.atoms = eval(d.pop('atoms'), vars(elements))
        for atom in obj.atoms:
            atom.building_block = bbs[atom.building_block]

        obj.bonds = eval(d.pop('bonds'), vars(bonds))
        for bond in obj.bonds:
            bond.atom1 = obj.atoms[bond.atom1]
            bond.atom2 = obj.atoms[bond.atom2]

        obj.construction_bonds = [
            obj.bonds[i] for i in d.pop('construction_bonds')
        ]

        g = {'FunctionalGroup': FunctionalGroup}
        obj.func_groups = tuple(eval(d.pop('func_groups'), g))
        for func_group in obj.func_groups:
            func_group.atoms = tuple(
                obj.atoms[i] for i in func_group.atoms
            )
            func_group.bonders = tuple(
                obj.atoms[i] for i in func_group.bonders
            )
            func_group.deleters = tuple(
                obj.atoms[i] for i in func_group.deleters
            )
            func_group.fg_type = fg_types[func_group.fg_type]

        d.pop('class')
        for attr, val in d.items():
            setattr(obj, attr, eval(val))

        if use_cache:
            cls._cache[identity_key] = obj
        return obj

    @classmethod
    def _get_identity_key_from_components(
        cls,
        building_blocks,
        topology_graph,
        building_block_vertices,
    ):
        """
        Return the identity key.

        Parameters
        ----------
        building_blocks : :class:`list` of :class:`.Molecule`
            The :class:`.BuildingBlock` and
            :class:`ConstructedMolecule` instances which
            represent the building block molecules used for
            construction. Only one instance is present per building
            block molecule, even if multiples of that building block
            join up to form the :class:`ConstructedMolecule`.

        topology_graph : :class:`.TopologyGraph`
            Defines the topology graph of the
            :class:`ConstructedMolecule` and constructs it.

        building_block_vertices : :class:`dict`
            Maps the :class:`.Molecule` in  `building_blocks` to the
            :class:`~.topologies.base.Vertex` in `topology_graph`.
            Each :class:`.BuildingBlock` and
            :class:`ConstructedMolecule` can be mapped to multiple
            :class:`~.topologies.base.Vertex` objects. See the
            examples section in the :class:`.ConstructedMolecule`
            class docstring to help understand how this parameter
            is used. If ``None``, building block molecules will be
            assigned to vertices at random.

        Returns
        -------
        :class:`tuple`
            The identity key.

        """

        building_block_types = sorted(
            {bb.__class__ for bb in building_blocks},
            key=lambda cls: cls.__name__
        )
        bb_keys = cls._get_sorted_building_block_keys(
            building_blocks=building_blocks,
            building_block_types=building_block_types
        )
        vertices = cls._get_sorted_building_block_vertices(
            building_blocks=building_blocks,
            building_block_vertices=building_block_vertices,
            building_block_types=building_block_types
        )
        return tuple(bb_keys), repr(topology_graph), tuple(vertices)

    @classmethod
    def _get_sorted_building_block_keys(
        cls,
        building_blocks,
        building_block_types
    ):
        for bb_type in building_block_types:
            bbs = filter(
                lambda bb: bb.__class__ is bb_type,
                building_blocks
            )
            yield from sorted(x.get_identity_key() for x in bbs)

    @staticmethod
    def _get_building_block_vertices_of_type(
        building_block_vertices,
        building_block_type
    ):

        bbs = filter(
            lambda bb: bb.__class__ is building_block_type,
            building_block_vertices
        )
        for bb in bbs:
            vertices = tuple(sorted(
                v.id for v in building_block_vertices[bb]
            ))
            yield bb.get_identity_key(), vertices

    @classmethod
    def _get_sorted_building_block_vertices(
        cls,
        building_blocks,
        building_block_vertices,
        building_block_types
    ):
        # The identity key needs to account for which building block
        # is placed on which vertex. This means that a tuple matching
        # the identity key of each building block with each vertex
        # needs to be included. This tuple needs to be sorted so that
        # the identity key of a constructed molecule is always the
        # same. However, the constructed molecule can be constructed
        # with both BuildingBlocks and other ConstructedMolecules.
        # They each have differnt types as identity keys which means
        # they connot be sorted together. As a result, first
        # sort the identity keys of all buliding blocks of the same
        # type and then join these sorted tuples togther.
        # The sorted tuples are joined together in order based on the
        # names of the classes.
        for building_block_type in building_block_types:
            yield from sorted(
                cls._get_building_block_vertices_of_type(
                    building_block_vertices=building_block_vertices,
                    building_block_type=building_block_type
                )
            )

    def __str__(self):
        return (
            f'{self.__class__.__name__}'
            f'(building_blocks={list(self.get_building_blocks())}, '
            f'topology_graph={self.topology_graph!r})'
        )

    def __repr__(self):
        return str(self)
