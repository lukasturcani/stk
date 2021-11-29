"""
Metal Complex
=============

.. toctree::
    :maxdepth: 2

    Paddlewheel <\
stk.molecular.topology_graphs.metal_complex.paddlewheel.paddlewheel\
>
    Porphyrin <\
stk.molecular.topology_graphs.metal_complex.porphyrin.porphyrin\
>
    Octahedral <\
stk.molecular.topology_graphs.metal_complex.octahedral.octahedral\
>
    Octahedral Lambda <\
stk.molecular.topology_graphs.metal_complex.octahedral\
.octahedral_lambda\
>
    Octahedral Delta <\
stk.molecular.topology_graphs.metal_complex.octahedral\
.octahedral_delta\
>
    Square Planar <\
stk.molecular.topology_graphs.metal_complex.square_planar\
.square_planar\
>
    Bidentate Square Planar <\
stk.molecular.topology_graphs.metal_complex.square_planar\
.bidentate_square_planar\
>
    Cis Protected Square Planar <\
stk.molecular.topology_graphs.metal_complex.square_planar\
.cis_protected_square_planar\
>

"""

from collections import Counter, defaultdict
from itertools import product

from ...reactions import DativeReactionFactory, GenericReactionFactory
from ..topology_graph import NullOptimizer, TopologyGraph


class MetalComplex(TopologyGraph):
    """
    Represents a metal complex topology graph.

    Notes
    -----
    *Subclass Implementation*

    Each subclass needs to define the attributes,
    :attr:`_metal_vertex_prototypes` and
    :attr:`_ligand_vertex_prototypes`, which are :class:`tuple` of
    :class:`.Vertex` instances.

    Examples
    --------
    *Subclass Implementation*

    The source code of the subclasses, listed in
    :mod:`~.metal_complex.metal_complex`, can serve as good examples.

    *Basic Construction*

    For most :class:`.MetalComplex` topology graphs, we first
    need to define a metal :class:`.BuildingBlock`, consisting of
    1 atom and multiple functional groups

    .. testcode:: basic-construction

        import stk

        metal = stk.BuildingBlock(
            smiles='[Fe+2]',
            functional_groups=(
                stk.SingleAtom(stk.Fe(0, charge=2))
                for i in range(6)
            ),
            position_matrix=[[0, 0, 0]],
        )

    We also need to define an organic ligand :class:`.BuildingBlock`

    .. testcode:: basic-construction

        # Define an organic linker with two functional groups.
        bidentate = stk.BuildingBlock(
            smiles='C=NC/C=N/Br',
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#35]',
                    bonders=(1, ),
                    deleters=(),
                ),
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#6]',
                    bonders=(1, ),
                    deleters=(),
                ),
            ],
        )

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bidentate = stk.BuildingBlock(
            smiles='C=NC/C=N/Br',
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#35]',
                    bonders=(1, ),
                    deleters=(),
                ),
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#6]',
                    bonders=(1, ),
                    deleters=(),
                ),
            ],
        )

        moldoc_display_molecule = molecule.Molecule(
            atoms=(
                molecule.Atom(
                    atomic_number=atom.get_atomic_number(),
                    position=position,
                ) for atom, position in zip(
                    bidentate.get_atoms(),
                    bidentate.get_position_matrix(),
                )
            ),
            bonds=(
                molecule.Bond(
                    atom1_id=bond.get_atom1().get_id(),
                    atom2_id=bond.get_atom2().get_id(),
                    order=bond.get_order(),
                ) for bond in bidentate.get_bonds()
            ),
        )


    Finally, we can create the :class:`.MetalComplex`.

    .. testcode:: basic-construction

        complex = stk.ConstructedMolecule(
            topology_graph=stk.metal_complex.OctahedralLambda(
                metals=metal,
                ligands=bidentate,
            )
        )

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        metal = stk.BuildingBlock(
            smiles='[Fe+2]',
            functional_groups=(
                stk.SingleAtom(stk.Fe(0, charge=2))
                for i in range(6)
            ),
            position_matrix=[[0, 0, 0]],
        )

        bidentate = stk.BuildingBlock(
            smiles='C=NC/C=N/Br',
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#35]',
                    bonders=(1, ),
                    deleters=(),
                ),
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#6]',
                    bonders=(1, ),
                    deleters=(),
                ),
            ],
        )

        complex = stk.ConstructedMolecule(
            topology_graph=stk.metal_complex.OctahedralLambda(
                metals=metal,
                ligands=bidentate,
            )
        )

        moldoc_display_molecule = molecule.Molecule(
            atoms=(
                molecule.Atom(
                    atomic_number=atom.get_atomic_number(),
                    position=position,
                ) for atom, position in zip(
                    complex.get_atoms(),
                    complex.get_position_matrix(),
                )
            ),
            bonds=(
                molecule.Bond(
                    atom1_id=bond.get_atom1().get_id(),
                    atom2_id=bond.get_atom2().get_id(),
                    order=(
                        bond.get_order()
                        if bond.get_order() != 9
                        else 1
                    ),
                ) for bond in complex.get_bonds()
            ),
        )

    *Suggested Optimization*

    For :class:`.MetalComplex` topologies, it is recommend to use the
    :class:`.MCHammer` optimizer.

    .. testcode:: suggested-optimization

        import stk

        metal = stk.BuildingBlock(
            smiles='[Fe+2]',
            functional_groups=(
                stk.SingleAtom(stk.Fe(0, charge=2))
                for i in range(6)
            ),
            position_matrix=[[0, 0, 0]],
        )

        bidentate = stk.BuildingBlock(
            smiles='C=NC/C=N/Br',
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#35]',
                    bonders=(1, ),
                    deleters=(),
                ),
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#6]',
                    bonders=(1, ),
                    deleters=(),
                ),
            ]
        )

        complex = stk.ConstructedMolecule(
            topology_graph=stk.metal_complex.OctahedralLambda(
                metals=metal,
                ligands=bidentate,
                optimizer=stk.MCHammer(),
            ),
        )

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        metal = stk.BuildingBlock(
            smiles='[Fe+2]',
            functional_groups=(
                stk.SingleAtom(stk.Fe(0, charge=2))
                for i in range(6)
            ),
            position_matrix=[[0, 0, 0]],
        )

        bidentate = stk.BuildingBlock(
            smiles='C=NC/C=N/Br',
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#35]',
                    bonders=(1, ),
                    deleters=(),
                ),
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#6]',
                    bonders=(1, ),
                    deleters=(),
                ),
            ]
        )

        complex = stk.ConstructedMolecule(
            topology_graph=stk.metal_complex.OctahedralLambda(
                metals=metal,
                ligands=bidentate,
                optimizer=stk.MCHammer(),
            ),
        )

        moldoc_display_molecule = molecule.Molecule(
            atoms=(
                molecule.Atom(
                    atomic_number=atom.get_atomic_number(),
                    position=position,
                ) for atom, position in zip(
                    complex.get_atoms(),
                    complex.get_position_matrix(),
                )
            ),
            bonds=(
                molecule.Bond(
                    atom1_id=bond.get_atom1().get_id(),
                    atom2_id=bond.get_atom2().get_id(),
                    order=(
                        bond.get_order()
                        if bond.get_order() != 9
                        else 1
                    ),
                ) for bond in complex.get_bonds()
            ),
        )

    *Construction with Multiple Metals & Ligands*

    When multiple metals or ligands are used, the `metals` and
    `ligands` parameters accept values of type :class:`dict`, which
    specify the exact vertex each metal or ligand needs to be placed
    on.

    .. testcode:: construction-with-multiple-metals-and-ligands

        import stk

        metal = stk.BuildingBlock(
            smiles='[Fe+2]',
            functional_groups=(
                stk.SingleAtom(stk.Fe(0, charge=2))
                for i in range(6)
            ),
            position_matrix=[[0, 0, 0]],
        )

        bidentate1 = stk.BuildingBlock(
            smiles='C=NC/C=N/Br',
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#35]',
                    bonders=(1, ),
                    deleters=(),
                ),
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#6]',
                    bonders=(1, ),
                    deleters=(),
                ),
            ],
        )

        # Define a second organic linker with two functional groups.
        bidentate2 = stk.BuildingBlock(
            smiles='C=NC(C)(C)/C(C)=N/Br',
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#35]',
                    bonders=(1, ),
                    deleters=(),
                ),
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#6]',
                    bonders=(1, ),
                    deleters=(),
                ),
            ],
        )

        # Build heteroleptic complex.
        complex = stk.ConstructedMolecule(
            topology_graph=stk.metal_complex.OctahedralLambda(
                metals=metal,
                ligands={
                    bidentate1: (0, 1),
                    bidentate2: (2, ),
                },
            ),
        )

    However, if each ligand is has a different number of
    functional groups, they can be provided together in a
    :class:`tuple`.

    Note that the valid vertex identifiers depend on the exact
    metal complex you are using. These are detailed in the docstring
    for that specific metal vertex topology graph.

    *Unsubstituted Metal Complexes*

    Some metal complex topologies represent metal complexes with
    unsubstituted metal sites. For example,
    :class:`.BidentateSquarePlanar` has all sites substituted and
    :class:`.CisProtectedSquarePlanar` is the equivalent metal complex
    with some unsubstituted sites

    .. testcode:: leaving-unsubstituted-sites

        import stk

        pd = stk.BuildingBlock(
            smiles='[Pd+2]',
            functional_groups=(
                stk.SingleAtom(stk.Pd(0, charge=2))
                for i in range(4)
            ),
            position_matrix=[[0, 0, 0]],
        )

        # Define a bidentate ligand with two functional groups.
        bidentate_ligand = stk.BuildingBlock(
            smiles='NCCN',
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#7]~[#6]',
                    bonders=(0, ),
                    deleters=(),
                ),
            ],
        )

        # Construct a cis-protected square planar metal complex.
        complex = stk.ConstructedMolecule(
            topology_graph=stk.metal_complex.CisProtectedSquarePlanar(
                metals=pd,
                ligands=bidentate_ligand,
            ),
        )

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        pd = stk.BuildingBlock(
            smiles='[Pd+2]',
            functional_groups=(
                stk.SingleAtom(stk.Pd(0, charge=2))
                for i in range(4)
            ),
            position_matrix=[[0, 0, 0]],
        )

        # Define a bidentate ligand with two functional groups.
        bidentate_ligand = stk.BuildingBlock(
            smiles='NCCN',
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#7]~[#6]',
                    bonders=(0, ),
                    deleters=(),
                ),
            ],
        )

        # Construct a cis-protected square planar metal complex.
        complex = stk.ConstructedMolecule(
            topology_graph=stk.metal_complex.CisProtectedSquarePlanar(
                metals=pd,
                ligands=bidentate_ligand,
            ),
        )

        moldoc_display_molecule = molecule.Molecule(
            atoms=(
                molecule.Atom(
                    atomic_number=atom.get_atomic_number(),
                    position=position,
                ) for atom, position in zip(
                    complex.get_atoms(),
                    complex.get_position_matrix(),
                )
            ),
            bonds=(
                molecule.Bond(
                    atom1_id=bond.get_atom1().get_id(),
                    atom2_id=bond.get_atom2().get_id(),
                    order=(
                        bond.get_order()
                        if bond.get_order() != 9
                        else 1
                    ),
                ) for bond in complex.get_bonds()
            ),
        )

    """

    def __init_subclass__(cls, **kwargs):
        cls._vertex_degrees = Counter(
            vertex_id
            for edge in cls._edge_prototypes
            for vertex_id in edge.get_vertex_ids()
        )
        cls._metal_vertices_of_degree = defaultdict(set)
        for vertex in cls._metal_vertex_prototypes:
            degree = cls._vertex_degrees[vertex.get_id()]
            cls._metal_vertices_of_degree[degree].add(vertex.get_id())

        cls._ligand_vertices_of_degree = defaultdict(set)
        for vertex in cls._ligand_vertex_prototypes:
            degree = cls._vertex_degrees[vertex.get_id()]
            cls._ligand_vertices_of_degree[degree].add(vertex.get_id())

    def __init__(
        self,
        metals,
        ligands,
        reaction_factory=None,
        num_processes=1,
        optimizer=NullOptimizer(),
    ):
        """
        Initialize a :class:`.MetalComplex`.

        Parameters
        ----------
        metals : :class:`dict` or :class:`.BuildingBlock` or \
                :class:`tuple`
            Can be a :class:`dict` which maps the
            :class:`.BuildingBlock` instances to the indices of the
            vertices in :attr:`_metal_vertex_prototypes` it should
            be placed on.

            If each :class:`.BuildingBlock` has a different number
            of functional groups, they can be supplied together in
            a :class:`tuple`.

            Can also be a :class:`.BuildingBlock` instance, which
            should be placed on all :attr:`_metal_vertex_prototypes`
            on the topology graph.

        ligands : :class:`dict` or :class:`.BuildingBlock` or \
                :class:`tuple`
            Can be a :class:`dict` which maps the
            :class:`.BuildingBlock` instances to the indices of the
            vertices in :attr:`_ligand_vertex_prototypes` it should be
            placed on.

            If each :class:`.BuildingBlock` has a different number
            of functional groups, they can be supplied together in
            a :class:`tuple`.

            Can also be a :class:`.BuildingBlock` instance, which
            should be placed on all :attr:`_ligand_vertex_prototypes`
            on the topology graph.

        reaction_factory : :class:`.ReactionFactory`, optional
            The reaction factory to use for creating bonds between
            building blocks.
            By default, a :class:`.DativeReactionFactory` is used,
            which produces only dative bonds in any reactions done by
            this topology construction.

        num_processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        optimizer : :class:`.Optimizer`, optional
            Used to optimize the structure of the constructed
            molecule.

        """

        building_block_vertices = self._normalize_metals(metals)
        building_block_vertices.update(
            (building_block, vertices)
            for building_block, vertices
            in self._normalize_ligands(ligands).items()
        )

        # By default, assign a dative bond order to available
        # functional groups.
        if reaction_factory is None:
            metal_functional_groups = set(
                type(functional_group)
                for metal in self._normalize_metals(metals)
                for functional_group in metal.get_functional_groups()
            )
            ligand_functional_groups = set(
                type(functional_group)
                for ligand in self._normalize_ligands(ligands)
                for functional_group in ligand.get_functional_groups()
            )
            functional_group_pairs = product(
                metal_functional_groups,
                ligand_functional_groups,
            )
            reaction_factory = DativeReactionFactory(
                GenericReactionFactory(
                    bond_orders={
                        frozenset(pair): 9
                        for pair in functional_group_pairs
                    }
                )
            )

        super().__init__(
            building_block_vertices=building_block_vertices,
            edges=self._edge_prototypes,
            reaction_factory=reaction_factory,
            construction_stages=(),
            num_processes=num_processes,
            optimizer=optimizer,
            edge_groups=None,
        )

    def _normalize_metals(self, metals):
        """
        Return a map between metals and vertices.

        Parameters
        ----------
        metals : :class:`dict` or :class:`.BuildingBlock`
            The metal-based building blocks.

        Returns
        -------
        :class:`dict`
            Map of :class:`.BuildingBlock` to a :class:`tuple` of
            :class:`.Vertex`

        """

        if isinstance(metals, dict):
            metals_dict = {
                metal: tuple(self._get_metal_vertices(ids))
                for metal, ids in metals.items()
            }
        elif isinstance(metals, tuple):
            functional_group_counter = Counter(
                metal.get_num_functional_groups()
                for metal in metals
            )
            assert (
                all(
                    count == 1
                    for count
                    in functional_group_counter.values()
                )
            ), (
                'Cannot use a tuple when multiple metals '
                'have the same number of functional groups. '
                'Use a dictionary instead.'
            )
            metals_dict = {
                metal: tuple(
                    self._get_metal_vertices(
                        self._metal_vertices_of_degree[
                            metal.get_num_functional_groups()
                        ]
                    )
                )
                for metal in metals
            }
        else:
            ids = range(len(self._metal_vertex_prototypes))
            metals_dict = {
                metals: tuple(self._get_metal_vertices(ids))
            }

        return metals_dict

    def _normalize_ligands(self, ligands):
        """
        Return a map ligands and vertices.

        Parameters
        ----------
        ligands : :class:`dict` or :class:`.BuildingBlock` or \
                :class:`tuple`
            The organic-based building blocks.

        Returns
        -------
        :class:`dict`
            Map of :class:`.BuildingBlock` to a :class:`tuple` of
            :class:`.Vertex`

        Raises
        ------
        :class:`AssertionError`
            If a :class:`tuple` is provided for ligands but there is
            ambiguity on ligand-vertex assignment because two ligands
            have the same number of functional groups.

        """

        if isinstance(ligands, dict):
            ligands_dict = {
                ligand: tuple(self._get_ligand_vertices(ids))
                for ligand, ids in ligands.items()
            }
        elif isinstance(ligands, tuple):
            functional_group_counter = Counter(
                ligand.get_num_functional_groups()
                for ligand in ligands
            )
            assert (
                all(
                    count == 1
                    for count
                    in functional_group_counter.values()
                )
            ), (
                'Cannot use a tuple when multiple ligands '
                'have the same number of functional groups. '
                'Use a dictionary instead.'
            )
            ligands_dict = {
                ligand: tuple(
                    self._get_ligand_vertices(
                        self._ligand_vertices_of_degree[
                            ligand.get_num_functional_groups()
                        ]
                    )
                )
                for ligand in ligands
            }

        else:
            ids = range(len(self._ligand_vertex_prototypes))
            ligands_dict = {
                ligands: tuple(self._get_ligand_vertices(ids))
            }

        return ligands_dict

    def _get_metal_vertices(self, vertex_ids):
        """
        Yield vertex prototypes.

        Parameters
        ----------
        vertex_ids : :class:`iterable` of :class:`int`
            The ids of the vertices to yield.

        Yields
        ------
        :class:`.Vertex`
            A vertex prototype of the topology graph.

        """

        if isinstance(vertex_ids, int):
            vertex_ids = (vertex_ids, )

        for vertex_id in vertex_ids:
            yield self._metal_vertex_prototypes[vertex_id]

    def _get_ligand_vertices(self, vertex_ids):
        """
        Yield vertex prototypes.

        Parameters
        ----------
        vertex_ids : :class:`iterable` of :class:`int`
            The ids of the vertices to yield.

        Yields
        ------
        :class:`.Vertex`
            A vertex prototype of the topology graph.

        """

        if isinstance(vertex_ids, int):
            vertex_ids = (vertex_ids, )

        for vertex_id in vertex_ids:
            yield self._ligand_vertex_prototypes[vertex_id]

    def _get_scale(self, building_block_vertices):
        return 1

    def __repr__(self):
        return (
            f'metal_complex.{self.__class__.__name__}'
            f'()'
        )
