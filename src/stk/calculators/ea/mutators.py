"""
Mutation
========

#. :class:`RandomBuildingBlock`
#. :class:`.SimilarBuildingBlock`
#. :class:`.RandomTopology`
#. :class:`.RandomMutation`

Mutation is carried by :class:`.Mutator` objects. They inherit
:class:`.Mutator` and define a method :meth:`~Mutator.mutate`. This
method must take a single molecule and return a mutant.

Examples of how mutators work can be seen the documentation of
the various :class:`Mutator` classes, for example
:class:`RandomBuildingBlock`, :class:`SimilarBuildingBlock`
or :class:`RandomMutation`.


.. _`adding mutators`:

Making New Mutators
-------------------

Mutators must simple inherit the :class:`.Mutator` class and define a
method called :meth:`~Mutatator.mutate`, which take a single molecule
and returns a mutant molecule. There are no requirements besides this.

"""

import logging
import numpy as np


from ...utilities import dice_similarity


logger = logging.getLogger(__name__)


class MutationError(Exception):
    """
    Used for errors which occuring during mutation operations.

    """

    ...


class Mutator:
    """
    Creates mutants.

    """

    def __init__(self, use_cache):
        """
        Initialize a :class:`Mutator`.

        Parameters
        ----------
        use_cache : :class:`bool`
            Toggles use of the molecular cache.

        """

        self._use_cache = use_cache

    def set_cache_use(self, use_cache):
        """
        Set cache use on or off.

        Parameters
        ----------
        use_cache : :class:`bool`
            ``True`` if the cache is to be used.

        Returns
        -------
        None : :class:`NoneType`

        """

        self._use_cache = use_cache

    def mutate(self, mol):
        """
        Return a mutant of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be mutated.

        Returns
        -------
        mol : :class:`.Molecule`
            The mutant.

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method which must be implemented by
            a subclass.

        """

        raise NotImplementedError()


class RandomMutation(Mutator):
    """
    Uses a random :class:`Mutator` to carry out mutations.

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


        # Create the first mutator.
        topology_graphs = [
            stk.cage.TwoPlusThree(),
            stk.cage.EightPlusTwelve(),
            stk.cage.TwentyPlusThirty()
        ]
        random_topology = stk.RandomTopologyGraph(topology_graphs)

        # Create the second and third mutator.
        building_blocks = [
            stk.BuildingBlock('NC[Si]CCN', ['amine']),
            stk.BuildingBlock('NCCCCCCCN', ['amine']),
            stk.BuildingBlock('NC1CCCCC1N', ['amine'])
        ]

        random_bb = stk.RandomBuildingBlock(
            building_blocks=building_blocks,
            key=lambda mol: mol.func_groups[0].fg_type.name == 'amine'
        )

        similar_bb = stk.SimilarBuildingBlock(
            building_blocks=building_blocks,
            key=lambda mol: mol.func_groups[0].fg_type.name == 'amine'
        )

        # Create the mutator used to carry out the mutations.
        random_mutator = stk.RandomMutation(
            random_topology,
            random_bb,
            similar_bb
        )

        # Mutate a molecule, one of random_topology,
        # random_bb or similar_bb will be used.
        mutant1 = random_mutator.mutate(cage)

        # Mutate the molecule a second time, one of random_topology,
        # random_bb or similar_bb will be used.
        mutant2 = random_mutator.mutate(cage)

        # Mutate a mutant, one of random_topology,
        # random_bb or similar_bb will be used.
        mutant3 = random_mutator.mutate(mutant1)

    """

    def __init__(self, *mutators, weights=None, random_seed=None):
        """
        Initialize a :class:`RandomMutation` instance.

        Parameters
        ----------
        *mutators : :class:`Mutator`
            :class:`Mutator` objects which are used to carry out the
            mutations.

        weights : :class:`list` of :class:`float`, optional
            The probability that each :class:`Mutator` will be selected
            to carry out a mutation.

        random_seed : :class:`bool`, optional
            The random seed to use.

        """

        self._mutators = mutators
        self._weights = weights
        self._generator = np.random.RandomState(random_seed)

    def set_cache_use(self, use_cache):
        """
        Set cache use on or off.

        Parameters
        ----------
        use_cache : :class:`bool`
            ``True`` if the cache is to be used.

        Returns
        -------
        None : :class:`NoneType`

        """

        for mutator in self._mutators:
            mutator.set_cache_use(use_cache)
        super().set_cache_use(use_cache)

    def mutate(self, mol):
        """
        Return a mutant of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be mutated.

        Returns
        -------
        mol : :class:`.Molecule`
            The mutant.

        """

        mutator = self._generator.choice(
            a=self._mutators,
            p=self._weights
        )
        return mutator.mutate(mol)


class RandomBuildingBlock(Mutator):
    """
    Substitutes random building blocks.

    This mutator takes a :class:`.ConstructedMolecule` and substitutes
    the building blocks with one chosen at random from a given set.

    Examples
    --------
    .. code-block:: python

        import stk

        # Create a molecule which is to be mutated.
        bb1 = stk.BuildingBlock('NCCN', ['amine'])
        bb2 = stk.BuildingBlock('O=CCC=O', ['aldehyde'])
        polymer = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2],
            topology_graphs=stk.polymer.Linear('AB', [0, 0], n=3)
        )

        # Create molecules used to substitute building blocks.
        building_blocks = [
            stk.BuildingBlock('NC[Si]CCN', ['amine']),
            stk.BuildingBlock('NCCCCCCCN', ['amine']),
            stk.BuildingBlock('NC1CCCCC1N', ['amine'])
        ]

        # Create the mutator.
        random_bb = stk.RandomBuildingBlock(
            building_blocks=building_blocks,
            key=lambda mol: mol.func_groups[0].fg_type.name == 'amine'
        )

        # Mutate a molecule.
        mutant1 = random_bb.mutate(polymer)

        # Mutate the molecule a second time.
        mutant2 = random_bb.mutate(polymer)

        # Mutate a mutant.
        mutant3 = random_bb.mutate(mutant1)

    """

    def __init__(
        self,
        building_blocks,
        key,
        duplicate_building_blocks=False,
        random_seed=None,
        use_cache=False
    ):
        """
        Initialize a :class:`RandomBuildingBlock` instance.

        Parameters
        ----------
        building_blocks : :class:`tuple` of :class:`.Molecule`
            A group of molecules which are used to replace building
            blocks in molecules being mutated.

        key : :class:`function`
            A function which takes a :class:`.Molecule` and
            returns ``True`` or ``False``. This function is applied to
            every building block in the molecule being mutated.
            Building blocks which returned ``True`` are liable for
            substition by one of the molecules in `building_blocks`.

        duplicate_building_blocks : :class:`bool`, optional
            Indicates whether the building blocks used to construct the
            mutant must all be unique.

        random_seed : :class:`bool`, optional
            The random seed to use.

        use_cache : :class:`bool`, optional
            Toggles use of the molecular cache.

        """

        self._building_blocks = building_blocks
        self._key = key
        self._duplicate_building_blocks = duplicate_building_blocks
        self._generator = np.random.RandomState(random_seed)
        super().__init__(use_cache=use_cache)

    def mutate(self, mol):
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

        # Choose the building block which undergoes mutation.
        valid_bbs = [
            bb for bb in mol.building_block_vertices if self._key(bb)
        ]
        chosen_bb = self._generator.choice(valid_bbs)

        # If the mutant can have more than one of the same building
        # block, prevent only the building block getting replaced
        # from being used as a replacement.
        if self._duplicate_building_blocks:
            excluded_bbs = {chosen_bb}
        # If the mutant is to be composed of unique building blocks
        # only, prevent any building block already present in the
        # mol from being used as a replacement.
        else:
            excluded_bbs = set(mol.building_block_vertices.keys())
        # Make sure that the excluded_bbs will not be picked.
        mols = [
            mol for mol in self._building_blocks
            if mol not in excluded_bbs
        ]

        # Choose a replacement building block.
        replacement = self._generator.choice(mols)

        # Build the new ConstructedMolecule.
        new_bbs = [
            bb for bb in mol.building_block_vertices
            if bb is not chosen_bb
        ]
        new_bbs.append(replacement)
        return mol.__class__(
            building_blocks=new_bbs,
            topology_graph=mol.topology_graph,
            use_cache=self._use_cache
        )


class SimilarBuildingBlock(Mutator):
    """
    Substitutes similar building blocks.

    This mutator takes a :class:`.ConstructedMolecule` and substitutes
    the building blocks with the most similar one from a given set.

    Examples
    --------
    .. code-block:: python

        import stk

        # Create a molecule which is to be mutated.
        bb1 = stk.BuildingBlock('NCCN', ['amine'])
        bb2 = stk.BuildingBlock('O=CCC=O', ['aldehyde'])
        polymer = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2],
            topology_graph=stk.polymer.Linear('AB', [0, 0], n=3)
        )

        # Create molecules used to substitute building blocks.
        building_blocks = [
            stk.BuildingBlock('NC[Si]CCN', ['amine']),
            stk.BuildingBlock('NCCCCCCCN', ['amine']),
            stk.BuildingBlock('NC1CCCCC1N', ['amine'])
        ]

        # Create the mutator.
        similar_bb = stk.SimilarBuildingBlock(
            building_blocks=building_blocks,
            key=lambda mol: mol.func_groups[0].fg_type.name == 'amine'
        )

        # Mutate a molecule.
        mutant1 = random_bb.mutate(polymer)

        # Mutate the molecule a second time.
        mutant2 = random_bb.mutate(polymer)

        # Mutate a mutant.
        mutant3 = random_bb.mutate(mutant1)

    """

    def __init__(
        self,
        building_blocks,
        key,
        duplicate_building_blocks,
        random_seed=None,
        use_cache=False
    ):
        """
        Initialize a :class:`RandomBuildingBlock` instance.

        Parameters
        ----------
        building_blocks : :class:`tuple` of :class:`.Molecule`
            A group of molecules which are used to replace building
            blocks in molecules being mutated.

        key : :class:`function`
            A function which takes a :class:`.Molecule` and returns
            ``True`` or ``False``. This function is applied to every
            building block of the molecule being mutated. Building
            blocks which returned ``True`` are liable for substition by
            one of the molecules in :attr:`building_blocks`.

        duplicate_building_blocks : :class:`bool`, optional
            Indicates whether the building blocks used to construct the
            mutant must all be unique.

        random_seed : :class:`bool`, optional
            The random seed to use.

        use_cache : :class:`bool`, optional
            Toggles use of the molecular cache.

        """

        self._building_blocks = building_blocks
        self._key = key
        self._duplicate_building_blocks = duplicate_building_blocks
        self._similar_bbs = {}
        self._generator = np.random.RandomState(random_seed)
        super().__init__(use_cache=use_cache)

    def mutate(self, mol):
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

        if mol not in self._similar_bbs:
            # Maps the mol to a dict. The dict maps each
            # building block of the mol to an iterator.
            # The iterators yield the next most similar molecules in
            # `building_blocks` to the building block.
            self._similar_bbs[mol] = {}

        # Choose the building block which undergoes mutation.
        valid_bbs = [
            bb for bb in mol.building_block_vertices if self._key(bb)
        ]
        chosen_bb = self._generator.choice(valid_bbs)

        # If the building block has not been chosen before, create an
        # iterator yielding similar molecules from `building_blocks`
        # for it.
        if chosen_bb not in self._similar_bbs[mol]:
            self._similar_bbs[mol][chosen_bb] = iter(sorted(
                self._building_blocks,
                key=lambda m: dice_similarity(m, chosen_bb),
                reverse=True
            ))

        try:
            new_bb = next(self._similar_bbs[mol][chosen_bb])
        except StopIteration:
            self._similar_bbs[mol][chosen_bb] = iter(sorted(
                self._building_blocks,
                key=lambda m: dice_similarity(m, chosen_bb),
                reverse=True
            ))
            new_bb = next(self._similar_bbs[mol][chosen_bb])

        # If the most similar molecule in `mols` is itself, then take
        # the next most similar one.
        if new_bb is chosen_bb:
            try:
                new_bb = next(self._similar_bbs[mol][chosen_bb])
            except StopIteration:
                self._similar_bbs[mol][chosen_bb] = iter(sorted(
                    self._building_blocks,
                    key=lambda m: dice_similarity(m, chosen_bb),
                    reverse=True
                ))
                new_bb = next(self._similar_bbs[mol][chosen_bb])

        # Build the new ConstructedMolecule.
        new_bbs = [
            bb for bb in mol.building_block_vertices
            if bb is not chosen_bb
        ]
        new_bbs.append(new_bb)
        return mol.__class__(
            building_blocks=new_bbs,
            topology_graph=mol.topology_graph,
            use_cache=self._use_cache
        )


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

    def mutate(self, mol):
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
