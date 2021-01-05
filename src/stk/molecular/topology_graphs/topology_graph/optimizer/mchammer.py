"""
MCHammer
========

"""

from .optimizer import Optimizer
from .utilities import (
    get_mch_bond_topology,
    get_subunits,
    OptimizationIncompleteError,
)

import mchammer as mch


class MCHammer(Optimizer):
    """
    Perform Monte Carlo based optimisation of long-bonds in molecule.

    Implementation: [1]_

    .. [1] https://github.com/andrewtarzia/MCHammer.

    """

    def __init__(
        self,
        step_size=0.25,
        target_bond_length=1.2,
        num_steps=500,
        bond_epsilon=50,
        nonbond_epsilon=20,
        nonbond_sigma=1.2,
        nonbond_mu=3,
        beta=2,
        random_seed=None,
        save_trajectory=False,
    ):
        """
        Initialize an instance of :class:`.MCHammer`.

        Parameters
        ----------
        step_size : :class:`float`, optional
            The relative size of the step to take during step.
            Defaults to 0.25.

        target_bond_length : :class:`float`, optional
            Target equilibrium bond length for long bonds to minimize
            to.
            Defaults to 1.2 Angstrom.

        num_steps : :class:`int`, optional
            Number of MC moves to perform.
            Defaults to 500 steps.

        bond_epsilon : :class:`float`, optional
            Value of epsilon used in the bond potential in MC moves.
            Determines strength of the bond potential.
            Defaults to 50.

        nonbond_epsilon : :class:`float`, optional
            Value of epsilon used in the nonbond potential in MC moves.
            Determines strength of the nonbond potential.
            Defaults to 20.

        nonbond_sigma : :class:`float`, optional
            Value of sigma used in the nonbond potential in MC moves.
            Defaults to 1.2.

        nonbond_mu : :class:`float`, optional
            Value of mu used in the nonbond potential in MC moves.
            Determines the steepness of the nonbond potential.
            Defaults to 3.

        beta : :class:`float`, optional
            Value of beta used in the in MC moves. Beta takes the
            place of the inverse boltzmann temperature.
            Defaults to 2.

        random_seed : :class:`int`, optional
            Random seed to use for MC algorithm. Should only be set
            if exactly reproducible results are required, otherwise
            a system-based random seed should be used for proper
            sampling.

        save_trajectory : :class:`bool`, optional
            `True` to save MCHammer trajectory and information in
            :attr:`_trajectory_data`.
            Defaults to `False`.

        """
        self._optimizer = mch.Optimizer(
            step_size=step_size,
            target_bond_length=target_bond_length,
            num_steps=num_steps,
            bond_epsilon=bond_epsilon,
            nonbond_epsilon=nonbond_epsilon,
            nonbond_sigma=nonbond_sigma,
            nonbond_mu=nonbond_mu,
            beta=beta,
            random_seed=random_seed,
        )
        self._save_trajectory = save_trajectory
        self._trajectory_data = None

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
