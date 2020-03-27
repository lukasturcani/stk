"""
Jumble
======

"""

import numpy as np

from stk.utilities import dedupe
from .constructed_molecule import ConstructedMoleculeCrosser


class Jumble(ConstructedMoleculeCrosser):
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

    def _cross(self, *mols):
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
