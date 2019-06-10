"""
Defines crossers.

Crossers are objects which take a group of molecules and recombine
them to produce offspring molecules. How crossers are used can be
seen in the documentation of the various :class:`Crosser` classes,
for example :class:`GeneticRecombination`,
:class:`Jumble` or :class:`RandomCrossover`.

Available crossers.
-------------------

#. :class:`.RandomCrossover`
#. :class:`.GeneticRecombination`
#. :class:`.Jumble`

.. _`adding crossers`:

Extending stk: Making new crossers.
-----------------------------------

To add a new :class:`Crosser`, make a new class which inherits
:class:`Crosser` and defines a method called
:meth:`~Crosser.crossover`. The method is a generator, which can take
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

    def crossover(self, *mols):
        """
        Applies a crossover operation on some molecules.

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

    Attributes
    ----------
    crossers : :class:`tuple` of :class:`Crosser`
        :class:`Crosser` objects which are used to carry out the
        crossovers.

    weights : :class:`list` of :class:`float`
        The probability that each :class:`Crosser` will be selected to
        carry out a crossover.

    Examples
    --------
    Note that the crossover operations used in this example accept
    any number of parents.

    .. code-block:: python

        # Create the molecules which will be crossed.
        bb1 = StructUnit2.smiles_init('NCCN', ['amine'])
        bb2 = StructUnit2.smiles_init('O=CCCCC=O', ['aldehyde'])
        polymer1  = Polymer([bb1, bb2], Linear('AB', [0, 0], n=2))

        bb3 = StructUnit2.smiles_init('NCCCN', ['amine'])
        bb4 = StructUnit2.smiles_init('O=C[Si]CCC=O', ['aldehyde'])
        polymer2  = Polymer([bb3, bb4], Linear('AB', [0, 0], n=2))

        bb5 = StructUnit2.smiles_init('NC[Si]CN', ['amine'])
        bb6 = StructUnit2.smiles_init('O=CCNNCCC=O', ['aldehyde'])
        polymer3  = Polymer([bb5, bb6], Linear('AB', [0, 0], n=2))

        # Create the crossers.
        recombination = GeneticRecombination(
                         key=lambda mol: mol.func_group_infos[0].name)
        jumble = Jumble(num_offspring_building_blocks=2)
        random_crossover = RandomCrossover(recombination, jumble)

        # Get the offspring molecules, either recombination or jumble
        # will be used to make them.
        cohort1 = list(
            random_crossover.crossover(polymer1, polymer2, polymer3)
        )

        # Get a second set of offspring molecules, either recombination
        # or jumble will be used to make them.
        cohort2 = list(
            random_crossover.crossover(polymer1, polymer2, polymer3)
        )

        # Make a third set of offspring molecules by crossing two of
        # the offspring molecules. Either recombination or jumble will
        # be used to make them.
        offspring1, offspring2, *rest = cohort1
        cohort3 = list(
            random_crossover.crossover(offspring1, offspring2)
        )

    """

    def __init__(self, *crossers, weights=None):
        """
        Initializes a :class:`RandomCrossover` instance.

        Parameters
        ----------
        crossers : :class:`tuple` of :class:`Crosser`
            :class:`Crosser` objects which are used to carry out the
            crossovers.

        weights : :class:`list` of :class:`float`, optional
            The probability that each :class:`Crosser` will be selected
            to carry out a crossover.

        """

        self.crossers = crossers
        self.weights = weights

    def crossover(self, *args, **kwargs):
        """
        Create offspring with a :class:`Crosser` in :attr:`crossers`.

        Parameters
        ----------
        *args : :class:`object`
            Arguments passed to the chosen :class:`Crosser`'s
            :meth:`~Crosser.crossover` method.

        **kwargs : :class:`object`
            Keyword arguments passed to the chosen :class:`Crosser`'s
            :meth:`~Crosser.crossover` method.

        Yields
        -------
        :class:`.Molecule`
            The generated offspring

        """

        crosser = np.random.choice(self.crossers, p=self.weights)
        yield from crosser.crossover(*args, **kwargs)


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

        bb1 = StructUnit2('filename1.mol', ['amine'])
        bb2 = StructUnit3('filename2.mol', ['aldehyde'])
        bb3 = StructUnit2('filename3.mol', ['aldehyde'])
        bb4 = StructUnit3('filename4.mol', ['amine'])

    We can define a function which analyzes a building block
    molecule and returns the gene it belongs to, for example

    .. code-block:: python

        def determine_gene(building_block):
            return building_block.func_group_infos[0].name

    Here, we can see that the gene to which each building block
    molecule belongs is given by the functional group name.
    Therefore there is an ``'amine'`` gene which has two alleles
    ``bb1`` and ``bb4`` and there is an ``'aldehyde'`` gene which
    has two alleles ``bb2`` and ``bb3``.

    Alternatively, we could have defined a function such as

    .. code-block:: python

        def determine_gene(building_block):
            return building_block.__class__.__name__

    Now we can see that we end up with the gene called
    ``'StructUnit2'`` which has two alleles ``bb1`` and ``bb3``
    and a second gene called ``'StructUnit3'`` which has the
    alleles ``bb2`` and ``bb4``.

    To produce offspring molecules, this class categorizes
    each building block of the parent molecules into genes using
    the `key` parameter. Then, to generate a single offspring, it
    picks a random building block for every gene. The picked
    building blocks are used to construct the offspring. The
    topoogy of the offspring is one of the parent's topologies.
    For obvious reasons, this approach works with any number of
    parents.

    Attributes
    ----------
    key : :class:`function`
        A function, which takes a :class:`.StructUnit` object
        and returns its gene or category. To produce an offspring,
        one of the building blocks from each category is picked
        at random.

    random_yield_order : :class:`bool`
        Toggles if the order in which the building blocks are made
        should be random.

    Examples
    --------
    Note that any number of parents can be used for the crossover.

    .. code-block:: python

        # Create the molecules which will be crossed.
        bb1 = StructUnit2.smiles_init('NCCN', ['amine'])
        bb2 = StructUnit2.smiles_init('O=CCCCC=O', ['aldehyde'])
        polymer1  = Polymer([bb1, bb2], Linear('AB', [0, 0], n=2))

        bb3 = StructUnit2.smiles_init('NCCCN', ['amine'])
        bb4 = StructUnit2.smiles_init('O=C[Si]CCC=O', ['aldehyde'])
        polymer2  = Polymer([bb3, bb4], Linear('AB', [0, 0], n=2))

        bb5 = StructUnit2.smiles_init('NC[Si]CN', ['amine'])
        bb6 = StructUnit2.smiles_init('O=CCNNCCC=O', ['aldehyde'])
        polymer3  = Polymer([bb5, bb6], Linear('AB', [0, 0], n=2))

        # Create the crosser.
        recombination = GeneticRecombination(
                         key=lambda mol: mol.func_group_infos[0].name)

        # Get the offspring molecules.
        cohort1 = list(
            recombination.crossover(polymer1, polymer2, polymer3)
        )

        # Get a second set of offspring molecules.
        cohort2 = list(
            recombination.crossover(polymer1, polymer2, polymer3)
        )

        # Make a third set of offspring molecules by crossing two of
        # the offspring molecules.
        offspring1, offspring2, *rest = cohort1
        cohort3 = list(
            recombination.crossover(offspring1, offspring2)
        )

    """

    def __init__(self, key, random_yield_order=True):
        """
        Initializes a :class:`GeneticRecombination` instance.

        Parameters
        ----------
        key : :class:`function`
            A function, which takes a :class:`.StructUnit` object
            and returns its gene or category. To produce an offspring,
            one of the building blocks from each category is picked
            at random.

        random_yield_order : :class:`bool`, optional
            Toggles if the order in which the building blocks are
            made should be random.

        """

        self.key = key
        self.random_yield_order = random_yield_order
        super().__init__()

    def crossover(self, *mols):
        """
        Applies a crossover operation on some molecules.

        Parameters
        ----------
        *mols : :class:`.MacroMolecule`
            The molecules to be crossed.

        Yields
        ------
        :class:`.MacroMolecule`
            The generated offspring.

        """

        cls = mols[0].__class__

        genes = defaultdict(set)
        for mol in mols:
            for allele in mol.building_blocks:
                genes[self.key(allele)].add(allele)

        genes = {gene: np.random.permutation(list(alleles))
                 for gene, alleles in genes.items()}

        tops = dedupe((mol.topology for mol in mols), key=repr)

        product = it.product(*genes.values(), tops)

        if self.random_yield_order:
            product = list(product)
            np.random.shuffle(product)

        for *building_blocks, top in product:
            yield cls(building_blocks, top)


class Jumble(Crosser):
    """
    Distributes all building blocks among offspring.

    Puts all the building blocks from each parent into one big pot
    and building blocks are drawn from the pot to generate the
    offspring. The offspring inherit the topology of one of the
    parents.

    Attributes
    ----------
    num_offspring_building_blocks : :class:`int`
        The number of building blocks each offspring is made from.

    duplicate_building_blocks : :class:`bool`
        Indicates whether the building blocks used to construct the
        offspring must all be unique.

    random_yield_order : :class:`bool`
        Toggles if the order in which the building blocks are made
        should be random.

    Examples
    --------
    Note that any number of parents can be used for the crossover.

    .. code-block:: python

        # Create the molecules which will be crossed.
        bb1 = StructUnit2.smiles_init('NCCN', ['amine'])
        bb2 = StructUnit2.smiles_init('O=CCCCC=O', ['aldehyde'])
        polymer1  = Polymer([bb1, bb2], Linear('AB', [0, 0], n=2))

        bb3 = StructUnit2.smiles_init('NCCCN', ['amine'])
        bb4 = StructUnit2.smiles_init('O=C[Si]CCC=O', ['aldehyde'])
        polymer2  = Polymer([bb3, bb4], Linear('AB', [0, 0], n=2))

        bb5 = StructUnit2.smiles_init('NC[Si]CN', ['amine'])
        bb6 = StructUnit2.smiles_init('O=CCNNCCC=O', ['aldehyde'])
        polymer3  = Polymer([bb5, bb6], Linear('AB', [0, 0], n=2))

        # Create the crosser.
        jumble = Jumble(num_offspring_building_blocks=2)

        # Get the offspring molecules.
        cohort1 = list(jumble.crossover(polymer1, polymer2, polymer3))

        # Get a second set of offspring molecules.
        cohort2 = list(jumble.crossover(polymer1, polymer2, polymer3))

        # Make a third set of offspring molecules by crossing two of
        # the offspring molecules.
        offspring1, offspring2, *rest = cohort1
        cohort3 = list(jumble.crossover(offspring1, offspring2))

    """

    def __init__(self,
                 num_offspring_building_blocks,
                 duplicate_building_blocks=False,
                 random_yield_order=True):
        """
        Initializes a :class:`Jumble` instance.

        Parameters
        ----------
        num_offspring_building_blocks : :class:`int`
            The number of building blocks each offspring is made from.

        duplicate_building_blocks : :class:`bool`, optional
            Indicates whether the building blocks used to construct the
            offspring must all be unique.

        random_yield_order : :class:`bool`, optional
            Toggles if the order in which the building blocks are
            made should be random.

        """

        n = num_offspring_building_blocks
        self.num_offspring_building_blocks = n
        self.duplicate_building_blocks = duplicate_building_blocks
        self.random_yield_order = random_yield_order
        super().__init__()

    def crossover(self, *mols):
        """
        Carries out crossover on some molecules.

        Parameters
        ----------
        *mols : :class:`.MacroMolecule`
            The parent molecules.

        Yields
        ------
        :class:`.MacroMolecule`
            The offspring molecule.

        """

        cls = mols[0].__class__
        building_blocks = dedupe(
            bb for mol in mols for bb in mol.building_blocks
        )

        if self.duplicate_building_blocks:
            combinations = it.combinations_with_replacement
        else:
            combinations = it.combinations
        building_block_groups = combinations(
                                  iterable=building_blocks,
                                  r=self.num_offspring_building_blocks)
        topologies = dedupe((mol.topology for mol in mols), key=repr)
        product = it.product(building_block_groups, topologies)

        if self.random_yield_order:
            product = list(product)
            np.random.shuffle(product)

        for bbs, top in product:
            yield cls(bbs, top)
