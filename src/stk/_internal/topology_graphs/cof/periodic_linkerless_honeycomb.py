"""
Periodic Linkerless Honeycomb
=============================

"""


import numpy as np

from stk._internal.optimizers.null import NullOptimizer
from stk._internal.reaction_factories.generic_reaction_factory import (
    GenericReactionFactory,
)
from stk._internal.topology_graphs.edge import Edge

from .cof import Cof
from .vertices import NonLinearVertex


class PeriodicLinkerlessHoneycomb(Cof):
    """
    Represents a periodic honeycomb COF topology graph.

    Unoptimized construction

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        cof = stk.ConstructedMolecule(
            topology_graph=stk.cof.PeriodicLinkerlessHoneycomb(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrCC(CBr)CBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                lattice_size=(3, 3, 1),
            ),
        )
        moldoc_display_molecule = molecule.Molecule(
            atoms=(
                molecule.Atom(
                    atomic_number=atom.get_atomic_number(),
                    position=position,
                ) for atom, position in zip(
                    cof.get_atoms(),
                    cof.get_position_matrix(),
                )
            ),
            bonds=(
                molecule.Bond(
                    atom1_id=bond.get_atom1().get_id(),
                    atom2_id=bond.get_atom2().get_id(),
                    order=bond.get_order(),
                ) for bond in cof.get_bonds()
                if all(p == 0 for p in bond.get_periodicity())
            ),
        )

    ``Collapser(scale_steps=False)`` optimized construction

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        cof = stk.ConstructedMolecule(
            topology_graph=stk.cof.PeriodicLinkerlessHoneycomb(
                building_blocks=(
                    stk.BuildingBlock(
                        smiles='BrCC(CBr)CBr',
                        functional_groups=[stk.BromoFactory()],
                    ),
                ),
                lattice_size=(3, 3, 1),
                optimizer=stk.Collapser(scale_steps=False),
            ),
        )
        moldoc_display_molecule = molecule.Molecule(
            atoms=(
                molecule.Atom(
                    atomic_number=atom.get_atomic_number(),
                    position=position,
                ) for atom, position in zip(
                    cof.get_atoms(),
                    cof.get_position_matrix(),
                )
            ),
            bonds=(
                molecule.Bond(
                    atom1_id=bond.get_atom1().get_id(),
                    atom2_id=bond.get_atom2().get_id(),
                    order=bond.get_order(),
                ) for bond in cof.get_bonds()
                if all(p == 0 for p in bond.get_periodicity())
            ),
        )

    Building blocks with three functional groups are required
    for this topology graph.

    When using a :class:`dict` for the `building_blocks` parameter,
    as in :ref:`cof-topology-graph-examples`:
    *Multi-Building Block COF Construction*, a
    :class:`.BuildingBlock`, with the following number of functional
    groups, needs to be assigned to each of the following vertex ids:

        | 3-functional groups: 0 to 1

    Note that optimizers may not optimize the :class:`.PeriodicInfo`.
    The documentation of the optimizer will state if it does.

    See :class:`.Cof` for more details and examples.

    """

    def __init__(
        self,
        building_blocks,
        lattice_size,
        vertex_alignments=None,
        reaction_factory=GenericReactionFactory(),
        num_processes=1,
        optimizer=NullOptimizer(),
    ):
        """
        Initialize a :class:`.PeriodicLinkerlessHoneycomb` instance.

        Parameters
        ----------
        building_blocks : :class:`tuple` or :class:`dict`
            Can be a :class:`tuple` of :class:`.BuildingBlock`
            instances, which should be placed on the topology graph.

            Can also be a :class:`dict` which maps the
            :class:`.BuildingBlock` instances to the ids of the
            vertices it should be placed on. A :class:`dict` is
            required when there are multiple building blocks with the
            same number of functional groups, because in this case
            the desired placement is ambiguous.

        lattice_size : :class:`tuple` of :class:`int`
            The size of the lattice in the x, y and z directions.

        vertex_alignments : :class:`dict`, optional
            A mapping from the id of a :class:`.Vertex`
            to an :class:`.Edge` connected to it.
            The :class:`.Edge` is used to align the first
            :class:`.FunctionalGroup` of a :class:`.BuildingBlock`
            placed on that vertex. Only vertices which need to have
            their default edge changed need to be present in the
            :class:`dict`. If ``None`` then the default edge is used
            for each vertex. Changing which :class:`.Edge` is used will
            mean that the topology graph represents different
            structural isomers. The edge is referred to by a number
            between ``0`` (inclusive) and the number of edges the
            vertex is connected to (exclusive).

        reaction_factory : :class:`.ReactionFactory`, optional
            The reaction factory to use for creating bonds between
            building blocks.

        num_processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        optimizer : :class:`.Optimizer`, optional
            Used to optimize the structure of the constructed
            molecule.

        Raises
        ------
        :class:`AssertionError`
            If the any building block does not have a
            valid number of functional groups.

        :class:`ValueError`
            If the there are multiple building blocks with the
            same number of functional_groups in `building_blocks`,
            and they are not explicitly assigned to vertices. The
            desired placement of building blocks is ambiguous in
            this case.

        :class:`~.cof.UnoccupiedVertexError`
            If a vertex of the COF topology graph does not have a
            building block placed on it.

        :class:`~.cof.OverlyOccupiedVertexError`
            If a vertex of the COF topology graph has more than one
            building block placed on it.

        """

        super().__init__(
            building_blocks=building_blocks,
            lattice_size=lattice_size,
            periodic=True,
            vertex_alignments=vertex_alignments,
            reaction_factory=reaction_factory,
            num_processes=num_processes,
            optimizer=optimizer,
        )

    _lattice_constants = _a, _b, _c = (
        np.array([1.0, 0.0, 0.0]),
        np.array([0.5, 0.866, 0.0]),
        np.array([0.0, 0.0, 5 / 1.7321]),
    )

    _vertex_prototypes = (
        NonLinearVertex(0, (1 / 3) * _a + (1 / 3) * _b + (1 / 2) * _c),
        NonLinearVertex(1, (2 / 3) * _a + (2 / 3) * _b + (1 / 2) * _c),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[1]),
        Edge(
            id=1,
            vertex1=_vertex_prototypes[0],
            vertex2=_vertex_prototypes[1],
            periodicity=(-1, 0, 0),
        ),
        Edge(
            id=2,
            vertex1=_vertex_prototypes[0],
            vertex2=_vertex_prototypes[1],
            periodicity=(0, -1, 0),
        ),
    )
