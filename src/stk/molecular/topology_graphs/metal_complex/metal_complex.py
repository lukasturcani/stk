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
stk.molecular.topology_graphs.metal_complex.octahedral.octahedral_lambda\
>
    Octahedral Delta <\
stk.molecular.topology_graphs.metal_complex.octahedral.octahedral_delta\
>
    Square Planar <\
stk.molecular.topology_graphs.metal_complex.square_planar.square_planar\
>
    Bidentate Square Planar <\
stk.molecular.topology_graphs.metal_complex.square_planar.bidentate_square_planar\
>
    Cis Protected Square Planar <\
stk.molecular.topology_graphs.metal_complex.square_planar.cis_protected_square_planar\
>

"""

from ..topology_graph import TopologyGraph
from ...reactions import GenericReactionFactory


class MetalComplex(TopologyGraph):
    """
    Represents a metal complex topology graph.

    Examples
    --------
    *Subclass Implementation*

    The source code of the subclasses, listed in
    :mod:`~.metal_complex.metal_complex`, can serve as good examples.

    Each subclass needs to define the iterable attributes,
    :attr:`_metal_vertex_prototypes` and
    :attr:`_ligand_vertex_prototypes`, which define the positions for
    metal-based and ligand-based vertices, respectively.

    *Construction*

    In some cases, constructing metal-containing archiectures requires
    more explicit input than other :class:`.TopologyGraph`. For
    example, to build a six-coordinate iron(II) complex with an
    octahedral geometry and lambda stereochemistry, we firstly need to
    define a metal atom (with six functional groups) and the organic
    ligand :class:`.BuildingBlock` (with two functional groups).

    .. code-block:: python
        import stk

        # Define a single atom of desired smiles with coordinates
        # (0, 0, 0) using RDKit.
        atom = rdkit.MolFromSmiles('[Fe+2]')
        atom.AddConformer(rdkit.Conformer(atom.GetNumAtoms()))
        stk_atom = stk.BuildingBlock.init_from_rdkit_mol(atom)
        atom_0, = stk_atom.get_atoms(0)

        # Assign the atom with 6 functional groups.
        stk_atom = stk_atom.with_functional_groups(
            (stk.SingleAtom(atom_0) for i in range(6))
        )

        # Define an organic linker with two functional groups.
        # In this example, the ordering of the functional groups is set
        # such that the ligand orientation is consistent and the target
        # stereochemistry is obtained.
        bidentate_ligand = stk.BuildingBlock(
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

    Next we perform construction of the :class:`.TopologyGraph`. Note
    that :class:`.MetalComplex` topologies have a different
    initializer, that clearly defines the placement of metal and ligand
    building blocks. Each :class:`.MetalComplex` provides example
    placements in their docstring.
    It is crucial that metal-ligand bonds are correctly created to be
    dative (i.e. the valence of the ligand `bonder` is not affected by
    the bond). This is handled by providing the
    :class:`DativeReactionFactory` a :class:`GenericReactionFactory`
    with the bond order of the expected functional groups defined to be
    9.

    .. code-block:: python

        complex = stk.ConstructedMolecule(
            stk.metal_complex.OctahedralLambda(
                metals={stk_atom: 0},
                ligands={bidentate_ligand: (0, 1, 2)},
                reaction_factory=stk.DativeReactionFactory(
                    stk.GenericReactionFactory(
                        bond_orders={
                            frozenset({
                                stk.GenericFunctionalGroup,
                                stk.SingleAtom
                            }): 9
                        }
                    )
                )
            )
        )

    *Leaving Unsubstitued Sites*

    Sometimes when building metal complexes to be used as
    :class:`.BuildingBlock` for a :class:`.ConstructedMolecule`, it is
    necessary to build a metal with open metal sites. Here, we provide
    the :class:`.CisProtectedSquarePlanar` topology to allow the user
    to build a square planar palladium(II) complex with two open metal
    sites. Notice in this example, that the metal and ligand placement
    is not specified upon initialisation. This is possible for all
    topologies.

    .. code-block:: python

        import stk

        # Define single metal atom with four functional groups.
        pd_metal = build_single_atom(
            smiles='[Pd+2]',
            no_fgs=4
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
            ]
        )

        # Construct a cis-protected square planar metal complex.
        complex = stk.ConstructedMolecule(
            stk.metal_complex.CisProtectedSquarePlanar(
                metals=pd_metal,
                ligands=bidentate_ligand,
                reaction_factory=stk.DativeReactionFactory(
                    stk.GenericReactionFactory(
                        bond_orders={
                            frozenset({
                                stk.GenericFunctionalGroup,
                                stk.SingleAtom
                            }): 9
                        }
                    )
                )
            )
        )
    Notes
    -----
    *Subclass Implementation*
    
    <Fill me in.>
    
    Examples
    --------
    *Subclass Implmentation*
    
    <Fill me in. You can just use the standard 
    disclaimer here, see Cage documentation.>
    
    *Construction*
    
    <Fill me in.>

    """

    def __init__(
        self,
        metals,
        ligands,
        reaction_factory=GenericReactionFactory(),
        num_processes=1,
    ):
        """
        Initialize a :class:`.MetalComplex`.

        Parameters
        ----------
        metals : :class:`dict` or :class:`.BuildingBlock`
            Can be a :class:`dict` which maps the
            :class:`.BuildingBlock` instances to the ids of the
            vertices in :attr:`_metal_vertex_prototypes` it should
            be placed on.

            Can also be a :class:`.BuildingBlock` instance, which
            should be placed at all :attr:`_metal_vertex_prototypes`
            on the topology graph.

        ligands : :class:`dict` or :class:`.BuildingBlock`
            Can be a :class:`dict` which maps the
            :class:`.BuildingBlock` instances to the ids of the
            vertices in :attr:`_ligand_vertex_prototypes` it should be
            placed on.

            Can also be a :class:`.BuildingBlock` instance, which
            should be placed at all :attr:`_ligand_vertex_prototypes`
            on the topology graph.

        reaction_factory : :class:`.ReactionFactory`, optional
            The reaction factory to use for creating bonds between
            building blocks.

        num_processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        """

        if isinstance(metals, dict):
            metals = {
                metal: self._get_metal_vertices(ids)
                for metal, ids in metals.items()
            }
        else:
            metals = {
                metals: self._get_metal_vertices(ids)
                for ids in range(len(self._metal_vertex_prototypes))
            }

        if isinstance(metals, dict):
            ligands = {
                ligand: self._get_ligand_vertices(ids)
                for ligand, ids in ligands.items()
            }
        else:
            ligands = {
                ligands: self._get_metal_vertices(ids)
                for ids in range(len(self._ligand_vertex_prototypes))
            }

        building_blocks = {**metals, **ligands}

        super().__init__(
            building_block_vertices=building_block_vertices,
            edges=self._edge_prototypes,
            reaction_factory=reaction_factory,
            construction_stages=(),
            num_processes=num_processes,
            edge_groups=None,
        )

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
