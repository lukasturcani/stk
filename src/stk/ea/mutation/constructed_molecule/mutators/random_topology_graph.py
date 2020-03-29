from .mutator import Mutator


class RandomTopologyGraph(Mutator):
    """
    Changes topology graphs at random.

    Examples
    --------
    .. code-block:: python

        import stk

        # Create a molecule which is to be mutated.
        bb1 = stk.BuildingBlock('NCCN', ['amine'])
        bb2 = stk.BuildingBlock('O=CCC(=O)CC=O', ['aldehyde'])
        cage = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2],
            topology_graph=stk.cage.FourPlusSix()
        )

        # Create topologies used for substition.
        topology_graphs = [
            stk.cage.TwoPlusThree(),
            stk.cage.EightPlusTwelve(),
            stk.cage.TwentyPlusThirty()
        ]

        # Create the mutator.
        random_topology = stk.RandomTopologyGraph(topology_graphs)

        # Mutate a molecule.
        mutant1 = random_topology.mutate(cage)

        # Mutate the molecule a second time.
        mutant2 = random_topology.mutate(cage)

        # Mutate a mutant.
        mutant3 = random_topology.mutate(mutant1)

    """

    def __init__(
        self,
        topology_graphs,
        random_seed=None,
        use_cache=False
    ):
        """
        Initialize a :class:`RandomTopology` instance.

        Parameters
        ----------
        topology_graphs : :class:`list` of :class:`.TopologyGraph`
            This lists holds the topology instances from which one is
            selected at random to form a new molecule.

        random_seed : :class:`bool`, optional
            The random seed to use.

        use_cache : :class:`bool`, optional
            Toggles use of the molecular cache.

        """

        self._topology_graphs = topology_graphs
        self._generator = np.random.RandomState(random_seed)
        super().__init__(use_cache=use_cache)

    def _mutate(self, mol):
        """
        Return a mutant of `mol`.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The molecule to be mutated.

        Returns
        -------
        mol : :class:`.ConstructedMolecule`
            The mutant.

        """

        tops = [
            x for x in self._topology_graphs
            if repr(x) != repr(mol.topology_graph)
        ]
        topology_graph = self._generator.choice(tops)
        return mol.__class__(
            building_blocks=list(mol.building_block_vertices.keys()),
            topology_graph=topology_graph,
            use_cache=self._use_cache
        )
