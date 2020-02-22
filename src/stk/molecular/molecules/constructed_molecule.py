"""
Constructed Molecule
====================

"""

import logging

from .molecule import Molecule

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
                    'position matrix\n'
                    '---------------\n'
                    f'{bb.get_position_matrix()}'
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
        self._bond_infos = construction_result.bond_infos
        self._building_block_counts = (
            construction_result.building_block_counts
        )

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
        clone._building_block_vertices = dict(
            self._building_block_vertices
        )
        clone._building_block_counts = dict(
            self._building_block_counts
        )
        clone._topology_graph = self._topology_graph
        clone._atom_infos = self._atom_infos
        clone._bond_infos = self._bond_infos
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

    def get_building_block_counts(self):
        """
        Get the count of each building block.

        Returns
        -------
        :class:`dict`
            Maps each :class:`.BuildingBlock` used during construction,
            to the number of times it is present in the
            :class:`.ConstructedMolecule`.

        """

        return dict(self._building_block_counts)

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

    def get_bond_infos(self):
        """
        Yield data about bonds in the molecule.

        Yields
        ------
        :class:`.BondInfo`
            Data about a bond.

        """

        yield from self._bond_infos

    def get_building_block_vertices(self):
        """
        Get vertices on which each building block was placed.

        Returns
        -------
        :class:`dict`
            Maps each :class:`.BuildingBlock` used during construction,
            to a :class:`tuple` of the :class:`.Vertex` instances on
            which it was placed.

        """

        return dict(self._building_block_vertices)

    def __str__(self):
        return (
            f'{self.__class__.__name__}'
            f'(building_blocks={list(self.get_building_blocks())}, '
            f'topology_graph={self.topology_graph!r})'
        )

    def __repr__(self):
        return str(self)
