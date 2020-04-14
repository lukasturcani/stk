"""
Genetic Recombination
=====================

"""

import itertools as it
from collections import defaultdict

from .crosser import MoleculeCrosser
from ...records import CrossoverRecord
from ....molecule_records import MoleculeRecord


class GeneticRecombination(MoleculeCrosser):
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

    In an :mod:`.stk` :class:`.ConstructedMolecule`, each building
    block represents an allele. The question is, which gene is each
    building block an allele of? To answer that, let's first construct
    a couple of building block molecules

    .. code-block:: python

        bb1 = stk.BuildingBlock(
            smiles='NCC(N)CN',
            functional_groups=[stk.PrimaryAminoFactory()],
        )
        bb2 = stk.BuildingBlock('O=CCC=O', [stk.AldehydeFactory()])
        bb3 = stk.BuildingBlock(
            smiles='O=CCNC(C=O)C=O',
            functional_groups=[stk.AldehydeFactory()],
        )
        bb4 = stk.BuildingBlock(
            smiles='NCOCN',
            functional_groups=[stk.PrimaryAminoFactory()],
        )

    We can define a function which analyzes a building block
    molecule and returns the gene it belongs to, for example

    .. code-block:: python

        def get_gene(building_block):
            fg, = building_block.get_functional_groups(0)
            return fg.__class__

    Here, we can see that the gene, to which each building block
    molecule belongs, is given by the class of its first functional
    group. Therefore there is an :class:`.PrimaryAmino` gene,
    which has two alleles ``bb1`` and ``bb4``, and there is an
    :class:`.Aldehyde` gene, which has two alleles ``bb2`` and ``bb3``.

    Alternatively, we could have defined a function such as

    .. code-block:: python

        def get_gene(building_block):
            return building_block.get_num_functional_groups()

    Now we can see that we end up with the gene called
    ``3``, which has two alleles ``bb1`` and ``bb3``,
    and a second gene called ``2``, which has the alleles ``bb2`` and
    ``bb4``.

    To produce offspring molecules, this class categorizes
    each building block of the parent molecules into genes using
    the `get_gene` parameter. Then, to generate a single offspring, it
    picks a building block for every gene. The picked
    building blocks are used to construct the offspring. The
    topology graph of the offspring is one of the parent's.
    For obvious reasons, this approach works with any number of
    parents.

    Examples
    --------
    *Crossing Constructed Molecules*

    Note that any number of parents can be used for the crossover

    .. code-block:: python

        import stk

        # Create the molecule records which will be crossed.

        bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock('O=CCCCC=O', [stk.AldehydeFactory()])
        graph1 = stk.polymer.Linear((bb1, bb2), 'AB', 2)
        polymer1  = stk.ConstructedMolecule(graph1)
        record1 = stk.MoleculeRecord(polymer1, graph1)

        bb3 = stk.BuildingBlock('NCCCN', [stk.PrimaryAminoFactory()])
        bb4 = stk.BuildingBlock(
            smiles='O=C[Si]CCC=O',
            functional_groups=[stk.AldehydeFactory()],
        )
        graph2 = stk.polymer.Linear((bb3, bb4), 'AB', 2)
        polymer2  = stk.ConstructedMolecule(graph2)
        record2 = stk.MoleculeRecord(polymer2, graph2)

        # Create the crosser.

        def get_functional_group_type(building_block):
            fg, = building_block.get_functional_groups(0)
            return fg.__class__

        recombination = stk.GeneticRecombination(
            get_gene=get_functional_group_type,
        )

        # Get the offspring molecules.

        cohort1 = tuple(recombination.cross(
            records=(record1, record2),
        ))

    """

    def __init__(
        self,
        get_gene,
        name='GeneticRecombination',
    ):
        """
        Initialize a :class:`GeneticRecombination` instance.

        Parameters
        ----------
        get_gene : :class:`callable`
            A :class:`callable`, which takes a :class:`.BuildingBlock`
            object and returns its gene. To produce an offspring, one
            of the building blocks from each gene is picked.

        name : :class:`str`, optional
            A name to identify the crosser instance.

        """

        self._get_gene = get_gene
        self._name = name

    def cross(self, records):
        topology_graphs = (
            record.get_topology_graph() for record in records
        )
        for topology_graph, alleles in it.product(
            topology_graphs,
            self._get_alleles(records),
        ):

            def get_replacement(building_block):
                gene = self._get_gene(building_block)
                return next(
                    allele for allele in alleles
                    if self._get_gene(allele) == gene
                )

            topology_graph = topology_graph.with_building_blocks(
                building_block_map={
                    building_block: get_replacement(building_block)
                    for building_block
                    in topology_graph.get_building_blocks()
                },
            )
            yield CrossoverRecord(
                molecule_record=MoleculeRecord(
                    topology_graph=topology_graph,
                ),
                crosser_name=self._name,
            )

    def _get_alleles(self, records):
        """
        Yield every possible combination of alleles.

        """

        genes = defaultdict(list)
        topology_graphs = (
            record.get_topology_graph() for record in records
        )
        for topology_graph in topology_graphs:
            for allele in topology_graph.get_building_blocks():
                genes[self._get_gene(allele)].append(allele)
        return it.product(*genes.values())
