"""
Constructed Molecule
====================

"""

import logging
import numpy as np
from collections import Counter

from .molecule import Molecule
from ..functional_groups import FunctionalGroup
from .. import atoms, topology_graphs

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

    Examples
    --------
    *Initialization*

    A :class:`ConstructedMolecule` can be created from a set of
    building blocks and a :class:`.TopologyGraph`

    .. code-block:: python

        import stk

        bb1 = stk.BuildingBlock('NCCCN', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock(
            smiles='O=CC(C=O)CC=O',
            functional_groups=[stk.AldehydeFactory()],
        )
        tetrahedron = stk.cage.FourPlusSix()
        cage1 = stk.ConstructedMolecule(
            building_blocks=(bb1, bb2),
            topology_graph=tetrahedron,
        )

    A :class:`ConstructedMolecule` can be used to construct other
    :class:`ConstructedMolecule` instances, but you have to convert
    them into a :class:`.BuildingBlock` first

    .. code-block:: python

        benzene = stk.BuildingBlock('c1ccccc1')
        cage_complex = stk.ConstructedMolecule(
            building_blocks=(
                stk.BuildingBlock.init_from_molecule(cage1),
                benzene,
            ),
            topology_graph=stk.host_guest.Complex(),
        )

    During initialization, it is possible to force building blocks to
    be placed on a specific :class:`.Vertex` of the
    :class:`.TopologyGraph` by specifying the vertex

    .. code-block:: python

        bb3 = stk.BuildingBlock('NCOCN', [stk.PrimaryAminoFactory()])
        bb4 = stk.BuildingBlock(
            smiles='NCOCCCOCN',
            functional_groups=[stk.PrimaryAminoFactory()],
        )
        cage2 = stk.ConstructedMolecule(
            building_blocks=(bb1, bb2, bb3, bb4),
            topology_graph=tetrahedron,
            building_block_vertices={
                bb1: tetrahedron.get_vertices(vertex_ids=(4, 5)),
                bb2: tetrahedron.get_vertices(vertex_ids=range(4)),
                bb3: tetrahedron.get_vertices(vertex_ids=6),
                bb4: tetrahedron.get_vertices(vertex_ids=range(7, 10)),
            },
        )

    """

    def __init__(
        self,
        building_blocks,
        topology_graph,
        building_block_vertices=None,
    ):
        """
        Initialize a :class:`.ConstructedMolecule`.

        Parameters
        ----------
        building_blocks : :class:`tuple` of :class:`.BuildingBlock`
            The :class:`.BuildingBlock` instances which
            represent the building block molecules used for
            construction. Only one instance is present per building
            block molecule, even if multiples of that building block
            join up to form the :class:`.ConstructedMolecule`.

        topology_graph : :class:`.TopologyGraph`
            Defines the topology graph of the
            :class:`.ConstructedMolecule` and constructs it.

        building_block_vertices : :class:`dict`, optional
            Maps the :class:`.BuildingBlock` in  `building_blocks` to
            the :class:`~.topologies.base.Vertex` instances in
            `topology_graph` it is placed on. Each
            :class:`.BuildingBlock` can be mapped to multiple
            :class:`~.topologies.base.Vertex`
            objects. See the examples section in the
            :class:`.ConstructedMolecule` class docstring to help
            understand how this parameter is used. If ``None``,
            a building block will be placed on an vertex with a
            degree equal to its number of functional groups.
            If multiple building blocks with the same number of
            functional groups are present in `building_blocks`, this
            parameter is required, as otherwise placement would be
            ambiguous.

        Raises
        ------
        :class:`.ConstructionError`
            If multiple building blocks with the same number of
            functional groups are present, and
            `building_block_vertices` is not provided.

        """

        if building_block_vertices is None:
            if self._is_placement_ambiguous(building_blocks):
                raise ConstructionError(
                    'Multiple building blocks have the same number '
                    'of functional groups and building_block_vertices '
                    'is None. Desired placement of building '
                    'blocks on vertices is therefore unclear. '
                    'Please use the building_block_vertices parameter '
                    'to fix this error.'
                )

            building_block_vertices = (
                topology_graph.assign_building_blocks_to_vertices(
                    building_blocks=building_blocks
                )
            )
        else:
            building_block_vertices = dict(building_block_vertices)

        try:
            construction_result = topology_graph.construct(
                building_block=building_block_vertices,
            )
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

        super().__init__(
            atoms=construction_result.atoms,
            bonds=construction_result.bonds,
            position_matrix=construction_result.position_matrix,
        )
        self._topology_graph = topology_graph
        self._building_block_vertices = building_block_vertices
        self._atom_infos = construction_result.atom_infos
        self._reaction_infos = construction_result.reaction_infos

    @staticmethod
    def _is_placement_ambiguous(building_blocks):
        """
        Return ``True``, if desired placement is unclear.

        The desired placement of `building_blocks` is unclear if
        multiple building blocks have the same number of functional
        groups. In cases like this, it is not clear which of these
        building blocks the user wants to place on a given vertex,
        since they are both valid candidates.

        Returns
        -------
        :class:`bool`
            ``True`` if placement is ambiguous and ``False``
            otherwise.

        """

        seen = set()
        for building_block in building_blocks:
            num_functional_groups = (
                building_block.get_num_functional_groups()
            )
            if num_functional_groups in seen:
                return True
            seen.add(num_functional_groups)
        return False

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
        :class:`.BuildingBlock`
            A building block of the :class:`ConstructedMolecule`.

        """

        yield from self._building_block_vertices.keys()

    def get_topology_graph(self):
        """
        Get the :class:`.TopologyGraph` used for construction.

        Returns
        -------
        :class:`.TopologyGraph`
            The topology graph used for construction.

        """

        return self._topology_graph

    def get_atom_infos(self, atom_ids=None):
        """
        Yield data about atoms in the molecule.

        Parameters
        ----------
        atom_ids : :class:`iterable` of :class:`int`, optional
            The ids of atoms whose data is desired. If ``None``,
            data on all atoms will be yielded. Can be a single
            :class:`int`, if data on a single atom is desired.

        Yields
        ------
        :class:`.AtomInfo`
            Data about an atom.

        """

        if atom_ids is None:
            atom_ids = range(len(self._atoms))
        elif isinstance(atom_ids, int):
            atom_ids = (atom_ids, )

        for atom_id in atom_ids:
            yield self._atom_infos[atom_id]

    def get_reaction_infos(self):
        """
        Yield data about reactions performed during construction.

        Yields
        ------
        :class:`.ReactionInfo`
            Data about a :class:`.Reaction` performed during
            construction.

        """

        yield from self._reaction_infos

    def get_building_block_vertices(self):
        """
        Get the

        Returns
        -------
        :class:`dict`

        """

        return dict(self._building_block_vertices)

    def to_dict(self):
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

    def init_from_dict(cls, molecule_dict):
        """
        Intialize from a :class:`dict` representation.

        Parameters
        ----------
        molecule_dict : :class:`dict`
            A :class:`dict` representation of a molecule generated
            by :meth:`to_dict`.

        Returns
        -------
        :class:`ConstructedMolecule`
            The molecule described by `mol_dict`.

        """

        d = dict(mol_dict)

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
                topology_graph.vertices[i] for i in vertices[str(i)]
            ]

        obj.topology_graph = topology_graph
        obj._position_matrix = np.array(d.pop('position_matrix')).T

        obj.atoms = eval(d.pop('atoms'), vars(atoms))
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
