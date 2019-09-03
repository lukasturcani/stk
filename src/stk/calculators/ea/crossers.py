"""
Crossover
=========

#. :class:`.RandomCrossover`
#. :class:`.GeneticRecombination`
#. :class:`.Jumble`

Crossover is implement through :class:`Crosser` objects. Crossers take
a group of molecules and recombine
them to produce offspring molecules. How crossers are used can be
seen in the documentation of the various :class:`Crosser` classes,
for example :class:`GeneticRecombination`,
:class:`Jumble` or :class:`RandomCrossover`.

.. _`adding crossers`:

Making New Crossers
-------------------

To add a new :class:`Crosser`, make a new class which inherits
:class:`Crosser` and defines a method called
:meth:`~Crosser.cross`. The method is a generator, which can take
any number of molecules and yields the offspring molecules.

"""

import logging
import numpy as np
import itertools as it
from collections import defaultdict

from ...utilities import dedupe


logger = logging.getLogger(__name__)


class Crosser:
    """
    Carries out crossover on molecules.

    """

    def __init__(self, use_cache):
        """
        Initialize a :class:`Crosser`.

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

    def cross(self, *mols):
        """
        Cross `mols`.

        Parameters
        ----------
        *mols : :class:`.Molecule`
            The molecules on which a crossover operation is performed.

        Yields
        -------
        :class:`.Molecule`
            The generated offspring.

        """

        return NotImplementedError()


class RandomCrossover(Crosser):
    """
    Uses a random :class:`Crosser` to carry out crossovers.

    Examples
    --------
    Note that the crossover operations used in this example accept
    any number of parents.

    .. code-block:: python

        import stk

        # Create the molecules which will be crossed.
        bb1 = stk.BuildingBlock('NCCN', ['amine'])
        bb2 = stk.BuildingBlock('O=CCCCC=O', ['aldehyde'])
        polymer1  = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2],
            topology_graph=stk.polymer.Linear('AB', [0, 0], n=2)
        )

        bb3 = stk.BuildingBlock('NCCCN', ['amine'])
        bb4 = stk.BuildingBlock('O=C[Si]CCC=O', ['aldehyde'])
        polymer2  = stk.ConstructedMolecule(
            building_blocks=[bb3, bb4],
            topology_graph=stk.polymer.Linear('AB', [0, 0], n=2)
        )

        bb5 = stk.BuildingBlock('NC[Si]CN', ['amine'])
        bb6 = stk.BuildingBlock('O=CCNNCCC=O', ['aldehyde'])
        polymer3  = stk.ConstructedMolecule(
            building_blocks=[bb5, bb6],
            topology_graph=stk.polymer.Linear('AB', [0, 0], n=2)
        )

        # Create the crossers.
        recombination = stk.GeneticRecombination(
            key=lambda mol: mol.func_groups[0].fg_type.name
        )
        jumble = stk.Jumble(num_offspring_building_blocks=2)
        random_crossover = stk.RandomCrossover(recombination, jumble)

        # Get the offspring molecules, either recombination or jumble
        # will be used to make them.
        cohort1 = list(
            random_crossover.cross(polymer1, polymer2, polymer3)
        )

        # Get a second set of offspring molecules, either recombination
        # or jumble will be used to make them.
        cohort2 = list(
            random_crossover.cross(polymer1, polymer2, polymer3)
        )

        # Make a third set of offspring molecules by crossing two of
        # the offspring molecules. Either recombination or jumble will
        # be used to make them.
        offspring1, offspring2, *rest = cohort1
        cohort3 = list(
            random_crossover.cross(offspring1, offspring2)
        )

    """

    def __init__(self, *crossers, weights=None, random_seed=None):
        """
        Initialize a :class:`RandomCrossover` instance.

        Parameters
        ----------
        *crossers : :class:`Crosser`
            :class:`Crosser` objects which are used to carry out the
            crossovers.

        weights : :class:`list` of :class:`float`, optional
            The probability that each :class:`Crosser` will be selected
            to carry out a crossover.

        random_seed : :class:`int`, optional
            The random seed to use.

        """

        self._crossers = crossers
        self._weights = weights
        self._generator = np.random.RandomState(random_seed)
        super().__init__(use_cache=False)

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

        for crosser in self._crossers:
            crosser.set_cache_use(use_cache)
        super().set_cache_use(use_cache)

    def cross(self, *mols):
        """
        Cross `mols`.

        Parameters
        ----------
        *mols : :class:`.Molecule`
            The molecules to be crossed.

        Yields
        -------
        :class:`.Molecule`
            A generated offspring

        """

        crosser = self._generator.choice(
            a=self._crossers,
            p=self._weights
        )
        yield from crosser.cross(*mols)


class GeneticRecombination(Crosser):
    """
    Recombine building blocks using biological systems as a model.

    Overall, this crosser mimics how animals and plants inherit
    DNA from their parents, except generalized to work with any
    number of parents. First it is worth discussing some
    terminology. A gene is a the smallest packet of genetic
    information. In animals, each gene can have multiple alleles.
    For example, there is a gene for hair color, and individual
    alleles for black, red, brown, etc. hair. This means that every
    person has a gene for hair color, but a person with black hair
    will have the black hair allele and a person with red hair will
    have the red hair allele. When two parents produce an
    offspring, the offspring will have a hair color gene and will
    inherit the allele of one of the parents at random. Therefore,
    if you have two parents, one with black hair and one with red
    hair, the offspring will either have black or red hair,
    depending on which allele they inherit.

    In ``stk`` molecules, each building block represents an allele.
    The question is, which gene is each building block an allele
    of? To answer that, let's first construct a couple of
    building block molecules

    .. code-block:: python

        bb1 = stk.BuildingBlock('NCC(N)CN', ['amine'])
        bb2 = stk.BuildingBlock('O=CCC=O', ['aldehyde'])
        bb3 = stk.BuildingBlock('O=CCNC(C=O)C=O', ['aldehyde'])
        bb4 = stk.BuildingBlock('NCOCN', ['amine'])

    We can define a function which analyzes a building block
    molecule and returns the gene it belongs to, for example

    .. code-block:: python

        def determine_gene(building_block):
            return building_block.func_groups[0].fg_type.name

    Here, we can see that the gene to which each building block
    molecule belongs is given by the functional group name.
    Therefore there is an ``'amine'`` gene which has two alleles
    ``bb1`` and ``bb4`` and there is an ``'aldehyde'`` gene which
    has two alleles ``bb2`` and ``bb3``.

    Alternatively, we could have defined a function such as

    .. code-block:: python

        def determine_gene(building_block):
            return len(building_block.func_groups)

    Now we can see that we end up with the gene called
    ``3`` which has two alleles ``bb1`` and ``bb3``
    and a second gene called ``2`` which has the alleles ``bb2`` and
    ``bb4``.

    To produce offspring molecules, this class categorizes
    each building block of the parent molecules into genes using
    the `key` parameter. Then, to generate a single offspring, it
    picks a random building block for every gene. The picked
    building blocks are used to construct the offspring. The
    topoogy graph of the offspring is one of the parent's.
    For obvious reasons, this approach works with any number of
    parents.

    Examples
    --------
    Note that any number of parents can be used for the crossover

    .. code-block:: python

        import stk

        # Create the molecules which will be crossed.
        bb1 = stk.BuildingBlock('NCCN', ['amine'])
        bb2 = stk.BuildingBlock('O=CCCCC=O', ['aldehyde'])
        polymer1  = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2],
            topolgy_graph=stk.polymer.Linear('AB', [0, 0], n=2)
        )

        bb3 = stk.BuildingBlock('NCCCN', ['amine'])
        bb4 = stk.BuildingBlock('O=C[Si]CCC=O', ['aldehyde'])
        polymer2  = stk.ConstructedMolecule(
            building_blocks=[bb3, bb4],
            topolog_graph=stk.polymer.Linear('AB', [0, 0], n=2)
        )

        bb5 = stk.BuildingBlock('NC[Si]CN', ['amine'])
        bb6 = stk.BuildingBlock('O=CCNNCCC=O', ['aldehyde'])
        polymer3  = stk.ConstructedMolecule(
            building_blocks=[bb5, bb6],
            topology_graph=stk.polymer.Linear('AB', [0, 0], n=2)
        )

        # Create the crosser.
        recombination = stk.GeneticRecombination(
            key=lambda mol: mol.func_groups[0].fg_type.name
        )

        # Get the offspring molecules.
        cohort1 = list(
            recombination.cross(polymer1, polymer2, polymer3)
        )

        # Get a second set of offspring molecules.
        cohort2 = list(
            recombination.cross(polymer1, polymer2, polymer3)
        )

        # Make a third set of offspring molecules by crossing two of
        # the offspring molecules.
        offspring1, offspring2, *rest = cohort1
        cohort3 = list(
            recombination.cross(offspring1, offspring2)
        )

    """

    def __init__(
        self,
        key,
        random_yield_order=True,
        random_seed=None,
        use_cache=False
    ):
        """
        Initialize a :class:`GeneticRecombination` instance.

        Parameters
        ----------
        key : :class:`callable`
            A :class:`callable`, which takes a :class:`.Molecule`
            object and returns its gene or category. To produce an
            offspring, one of the building blocks from each category is
            picked at random.

        random_seed : :class:`int`, optional
            The random seed to use.

        use_cache : :class:`bool`, optional
            Toggles use of the molecular cache.

        """

        self._key = key
        self._generator = np.random.RandomState(random_seed)
        super().__init__(use_cache=use_cache)

    def cross(self, *mols):
        """
        Cross `mols`.

        Parameters
        ----------
        *mols : :class:`.ConstructedMolecule`
            The molecules to crossed.

        Yields
        ------
        :class:`.ConstructedMolecule`
            The generated offspring.

        """

        cls = mols[0].__class__

        # Use a dict here to prevent the same structure for being
        # used in the same gene twice.
        genes = defaultdict(dict)
        for mol in mols:
            for allele in mol.building_block_vertices:
                gene = genes[self._key(allele)]
                if allele.get_identity_key() not in gene:
                    gene[allele.get_identity_key()] = allele

        genes = {
            gene: self._generator.permutation(list(alleles.values()))
            for gene, alleles in genes.items()
        }
        tops = dedupe((mol.topology_graph for mol in mols), key=repr)
        product = list(it.product(*genes.values(), tops))
        self._generator.shuffle(product)

        parents = {
            (
                *tuple(sorted(
                    bb.get_identity_key()
                    for bb in mol.building_block_vertices
                )),
                repr(mol.topology_graph)
            )
            for mol in mols
        }
        for *building_blocks, top in product:
            # Do not yield the parents.
            mol = (
                *tuple(sorted(
                    bb.get_identity_key() for bb in building_blocks
                )),
                repr(top)
            )
            if mol in parents:
                continue

            yield cls(
                building_blocks=building_blocks,
                topology_graph=top,
                use_cache=self._use_cache
            )


class Jumble(Crosser):
    """
    Distributes all building blocks among offspring.

    Puts all the building blocks from each parent into one big pot
    and building blocks are drawn from the pot to generate the
    offspring. The offspring inherit the topology  graph of one of the
    parents.

    Examples
    --------
    Note that any number of parents can be used for the crossover

    .. code-block:: python

        import stk

        # Create the molecules which will be crossed.
        bb1 = stk.BuildingBlock('NCCN', ['amine'])
        bb2 = stk.BuildingBlock('O=CCCCC=O', ['aldehyde'])
        polymer1  = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2],
            topology_graph=stk.polymer.Linear('AB', [0, 0], n=2)
        )

        bb3 = stk.BuildingBlock('NCCCN', ['amine'])
        bb4 = stk.BuildingBlock('O=C[Si]CCC=O', ['aldehyde'])
        polymer2  = stk.ConstructedMolecule(
            building_blocks=[bb3, bb4],
            topology_graph=stk.polymer.Linear('AB', [0, 0], n=2)
        )

        bb5 = stk.BuildingBlock('NC[Si]CN', ['amine'])
        bb6 = stk.BuildingBlock('O=CCNNCCC=O', ['aldehyde'])
        polymer3  = stk.ConstructedMolecule(
            building_blocks=[bb5, bb6],
            topology_graph=stk.polymer.Linear('AB', [0, 0], n=2)
        )

        # Create the crosser.
        jumble = stk.Jumble(num_offspring_building_blocks=2)

        # Get the offspring molecules.
        cohort1 = list(jumble.cross(polymer1, polymer2, polymer3))

        # Get a second set of offspring molecules.
        cohort2 = list(jumble.cross(polymer1, polymer2, polymer3))

        # Make a third set of offspring molecules by crossing two of
        # the offspring molecules.
        offspring1, offspring2, *rest = cohort1
        cohort3 = list(jumble.cross(offspring1, offspring2))

    """

    def __init__(
        self,
        num_offspring_building_blocks,
        duplicate_building_blocks=False,
        random_yield_order=True,
        random_seed=None,
        use_cache=False
    ):
        """
        Initialize a :class:`Jumble` instance.

        Parameters
        ----------
        num_offspring_building_blocks : :class:`int`
            The number of building blocks each offspring is made from.

        duplicate_building_blocks : :class:`bool`, optional
            Indicates whether the building blocks used to construct the
            offspring must all be unique.

        random_seed : :class:`int`, optional
            The random seed to use.

        use_cache : :class:`bool`, optional
            Toggles use of the molecular cache.

        """

        n = num_offspring_building_blocks
        self._num_offspring_building_blocks = n
        self._duplicate_building_blocks = duplicate_building_blocks
        self._generator = np.random.RandomState(random_seed)
        super().__init__(use_cache=use_cache)

    def cross(self, *mols):
        """
        Cross `mols`.

        Parameters
        ----------
        *mols : :class:`.ConstructedMolecule`
            The molecules to cross.

        Yields
        ------
        :class:`.ConstructedMolecule`
            An offspring molecule.

        """

        cls = mols[0].__class__
        building_blocks = dedupe(
            (bb for mol in mols for bb in mol.building_block_vertices),
            key=lambda bb: bb.get_identity_key()
        )

        if self._duplicate_building_blocks:
            combinations = it.combinations_with_replacement
        else:
            combinations = it.combinations
        building_block_groups = combinations(
            iterable=building_blocks,
            r=self._num_offspring_building_blocks
        )
        topologies = dedupe(
            (mol.topology_graph for mol in mols),
            key=repr
        )
        product = list(it.product(building_block_groups, topologies))
        self._generator.shuffle(product)

        parents = {
            (
                *tuple(sorted(
                    bb.get_identity_key()
                    for bb in mol.building_block_vertices
                )),
                repr(mol.topology_graph)
            )
            for mol in mols
        }

        for bbs, top in product:
            # Do not yield the parents.
            mol = (
                *tuple(sorted(bb.get_identity_key() for bb in bbs)),
                repr(top)
            )
            if mol in parents:
                continue

            yield cls(
                building_blocks=bbs,
                topology_graph=top,
                use_cache=self._use_cache
            )
