"""
Genetic Recombination
=====================

"""

from .constructed_molecule import ConstructedMoleculeCrosser


class GeneticRecombination(ConstructedMoleculeCrosser):
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
    topology graph of the offspring is one of the parent's.
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

    def _cross(self, *mols):
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
