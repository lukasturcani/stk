"""
Collapser
=========

"""

from .optimizer import Optimizer
from .utilities import (
    get_mch_bond_topology,
    get_subunits,
    OptimizationIncompleteError,
)

import mchammer as mch


class Collapser(Optimizer):
    """
    Perform rigid-body collapse of molecule.

    Implementation: [1]_

    .. [1] https://github.com/andrewtarzia/MCHammer.

    """

    def __init__(
        self,
        step_size=0.1,
        distance_threshold=1.5,
        scale_steps=True,
        save_trajectory=False,
    ):
        """
        Initialize an instance of :class:`.Collapser`.

        Parameters
        ----------
        step_size : :class:`float`, optional
            The relative size of the step to take during collapse.
            Defaults to 0.1 Angstrom.

        distance_threshold : :class:`float`, optional
            Distance between distinct subunits to use as
            threshold for halting collapse in Angstrom.
            Defaults to 1.5 Angstrom.

        scale_steps : :class:`bool`, optional
            Whether to scale the step of each distict building block
            by their relative distance from the molecules centroid.
            Defaults to ``True``

        save_trajectory : :class:`bool`, optional
            `True` to save Collapser trajectory and information in
            :attr:`_trajectory_data`.
            Defaults to `False`.

        """
        self._optimizer = mch.Collapser(
            step_size=step_size,
            distance_threshold=distance_threshold,
            scale_steps=scale_steps,
        )

    def optimize(self, state):
        """
        Optimize the structure of a molecule under construction.

        Parameters
        ----------
        state : :class:`.ConstructionState`
            The molecule being constructed.

        Returns
        -------
        :class:`.ConstructionState`
            The optimized construction state.

        """

        # Define MCHammer molecule to optimize.
        long_bond_ids, mch_bonds = get_mch_bond_topology(state)

        mch_mol = mch.Molecule(
            atoms=(
                mch.Atom(
                    id=atom.get_id(),
                    element_string=atom.__class__.__name__,
                ) for atom in state.get_atoms()
            ),
            bonds=mch_bonds,
            position_matrix=state.get_position_matrix(),
        )

        # Run optimization.
        return state.with_position_matrix(
            self._optimizer.get_result(
                mol=mch_mol,
                bond_pair_ids=long_bond_ids,
                subunits=get_subunits(state),
            ).get_final_position_matrix()
        )
