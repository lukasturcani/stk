"""
Jumble
======

"""

import numpy as np

from stk.utilities import dedupe
import itertools as it
from .constructed_molecule import ConstructedMoleculeCrosser
from .utilities import get_constructed_molecule_key


class Jumble(ConstructedMoleculeCrosser):
    """
    Distributes all building blocks among offspring.

    Puts all the building blocks from each parent into one big pot
    and building blocks are drawn from the pot to generate the
    offspring. The offspring inherit the topology  graph of one of the
    parents.

    Examples
    --------
    *Producing Offspring With an Arbitrary Number of Parents*

    Note that any number of parents can be used for the crossover

    .. code-block:: python

        import stk

        # Create the molecules which will be crossed.
        bb1 = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])
        bb2 = stk.BuildingBlock('BrCNCBr', [stk.BromoFactory()])
        polymer1  = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear((bb1, bb2), 'AB', 2)
        )

        bb3 = stk.BuildingBlock('BrCCCBr', [stk.BromoFactory()])
        bb4 = stk.BuildingBlock('BrC[Si]CCCBr', [stk.BromoFactory()])
        polymer2  = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear((bb3, bb4), 'AB', 2)
        )

        bb5 = stk.BuildingBlock('BrC[Si]CBr', [stk.BromoFactory()])
        bb6 = stk.BuildingBlock('BrCCNNCCCBr', [stk.BromoFactory()])
        polymer3  = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear((bb5, bb6), 'AB', 2)
        )

        # Create the crosser.
        jumble = stk.Jumble(num_offspring_building_blocks=2)

        # Get the offspring molecules.
        cohort1 = tuple(jumble.cross(
            molecules=(polymer1, polymer2, polymer3),
        )

        # Get a second set of offspring molecules.
        cohort2 = tuple(jumble.cross(
            molecules=(polymer1, polymer2, polymer3),
        )

        # Make a third set of offspring molecules by crossing two of
        # the offspring molecules.
        offspring1, offspring2, *rest = cohort1
        cohort3 = list(jumble.cross(
            molecules=(offspring1, offspring2),
        )

    """

    def __init__(
        self,
        num_offspring_building_blocks,
        duplicate_building_blocks=False,
        random_yield_order=True,
        random_seed=None,
        input_database=None,
        output_database=None,
        get_key=get_constructed_molecule_key,
    ):
        """
        Initialize a :class:`.Jumble` instance.

        Parameters
        ----------
        num_offspring_building_blocks : :class:`int`
            The number of building blocks each offspring is made from.

        duplicate_building_blocks : :class:`bool`, optional
            Indicates whether the building blocks used to construct the
            offspring must all be unique.

        random_yield_order : :class:`bool`, optional
            Toggles if the offspring produced by the crosser get
            yielded in a random order.

        random_seed : :class:`int`, optional
            The random seed to use.

        input_database : :class:`.ConstructedMoleculeDatabase`, \
                optional
            Before constructing molecules, this database is checked
            if they exist already, and returns them, if so.

        output_database : :class:`.ConstructedMoleculeDatabase`, \
                optional
            All molecules which get constructed, get deposited into
            this database.

        get_key : :class:`callable`, optional
            Takes a single parameter, a :class:`.TopologyGraph`
            and returns a key which is used lookup the corresponding
            :class:`.ConstructedMolecule` in the `input_database`.
            By default, :func:`.get_constructed_molecule_key` will
            be used.

        """

        self._num_offspring_building_blocks = (
            num_offspring_building_blocks
        )
        self._duplicate_building_blocks = duplicate_building_blocks
        self._generator = np.random.RandomState(random_seed)
        self._input_database = input_database
        self._output_database = output_database
        self._get_key = get_key

    def _cross(self, records):
        building_blocks = self._get_building_blocks(records)

        if self._duplicate_building_blocks:
            combinations = it.combinations_with_replacement
        else:
            combinations = it.combinations
        building_block_groups = combinations(
            iterable=building_blocks,
            r=self._num_offspring_building_blocks
        )
        topology_graphs = (
            record.get_topology_graph() for record in records
        )
        product = list(it.product(building_block_groups, topology_graphs))
        self._generator.shuffle(product)

        parents = {
            (
                *tuple(sorted(
                    bb.get_identity_key()
                    for bb in mol.building_block_vertices
                )),
                repr(mol.topology_graph)
            )
            for mol in molecules
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

    def _get_building_blocks(self, records):
        building_blocks = defaultdict(list)
        for record in records:

