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

    def __init__(self, topology_graph):
        """
        Initialize a :class:`.ConstructedMolecule`.

        Parameters
        ----------
        topology_graph : :class:`.TopologyGraph`
            Defines the topology graph of the
            :class:`.ConstructedMolecule` and constructs it.

        """

        construction_result = topology_graph.construct()
        super().__init__(
            atoms=construction_result.atoms,
            bonds=construction_result.bonds,
            position_matrix=construction_result.position_matrix,
        )
        self._topology_graph = topology_graph
        self._atom_infos = construction_result.atom_infos
        self._bond_infos = construction_result.bond_infos

    def clone(self):
        """
        Return a clone.

        Returns
        -------
        :class:`.ConstructedMolecule`
            The clone.

        """

        clone = super().clone()
        clone._topology_graph = self._topology_graph
        clone._atom_infos = self._atom_infos
        clone._bond_infos = self._bond_infos
        return clone

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

    def __str__(self):
        return f'{self.__class__.__name__}({self._topology_graph!r})'

    def __repr__(self):
        return str(self)
