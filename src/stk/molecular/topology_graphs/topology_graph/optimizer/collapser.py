"""
Collapser
=========

"""

from .optimizer import Optimizer
from .utilities import (
    get_mch_bond_topology,
    merge_subunits,
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
        subunits = mch_mol.get_subunits(bond_pair_ids=long_bond_ids)
        mch_result = self._optimizer.get_trajectory(
            mol=mch_mol,
            bond_pair_ids=long_bond_ids,
            # Merge subunits to match distinct BuildingBlocks in
            # stk ConstructedMolecule.
            subunits=merge_subunits(state, subunits),
        )

        if self._save_trajectory:
            self._trajectory_data = mch_result
        else:
            del mch_result
            del data

        return state.with_position_matrix(
            mch_result.get_final_position_matrix()
        )

    def get_trajectory_results(self):
        """
        Extract trajectory results after optimisation.

        Examples
        --------
        This optimisation's trajectory can be output in the
        :class:`mchammer.Result` data structure if `save_trajectory` is
        `True`. This data contains the structure at each step and
        properties of the structure that are being optimised. Here is
        an example of using this data

        .. code-block:: python

            coll_opt = stk.Collapser(
                scale_steps=False,
                save_trajectory=True,
            )
            polymer = stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(bb1, bb2),
                    repeating_unit='AB',
                    num_repeating_units=6,
                    optimizer=coll_opt,
                ),
            )
            polymer.write(f'polymer_opt.mol')

            mch_result = coll_opt.get_trajectory_results()

            # Output optimization log.
            with open(f'TESTING/optimization.out', 'w') as f:
                f.write(mch_result.get_log())

            # Output trajectory as separate xyz files for
            # visualisation. Note the use of a temporary stk.Molecule.
            for step, new_pos_mat in mch_result.get_trajectory():
                temp_ = polymer.with_position_matrix(new_pos_mat)
                temp_.write(
                    f'TESTING/traj_{step}.xyz'
                )

        Returns
        -------
        :class:`mchammer.Result`
            The trajectory data for the optimisation.

        Raises
        ------
        :class:`.OptimizationIncompleteError`
            Raises if optimization has not been completed.

        """

        if self._trajectory_data is None:
            raise OptimizationIncompleteError(
                'Optimization has not been run, so trajectory data is '
                'not available.'
            )

        return self._trajectory_data
