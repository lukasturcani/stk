"""
Collapser
=========

"""

from .optimizer import Optimizer
from .utilities import get_mch_bond_topology, get_subunits

import mchammer as mch


class Collapser(Optimizer):
    """
    Perform rigid-body collapse of molecule [1]_.

    Examples
    --------
    Optimisation with stk simply collects the final position matrix.
    The optimisation's trajectory can be output using the MCHammer
    implementation if required by the user [1]_.

    .. code-block:: python

        polymer = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(bb1, bb2),
                repeating_unit='AB',
                optimizer=stk.Collapser(scale_steps=False),
            ),
        )
        polymer.write(f'polymer_opt.mol')

    References
    ----------
    .. [1] https://github.com/andrewtarzia/MCHammer

    """

    def __init__(
        self,
        step_size=0.1,
        distance_threshold=1.5,
        scale_steps=True,
    ):
        """
        Initialize an instance of :class:`.Collapser`.

        Parameters
        ----------
        step_size : :class:`float`, optional
            The relative size of the step to take during collapse in
            Angstrom.

        distance_threshold : :class:`float`, optional
            Distance between distinct subunits to use as
            threshold for halting collapse in Angstrom.

        scale_steps : :class:`bool`, optional
            Whether to scale the step of each distict building block
            by their relative distance from the molecules centroid.

        """
        self._optimizer = mch.Collapser(
            step_size=step_size,
            distance_threshold=distance_threshold,
            scale_steps=scale_steps,
        )

    def optimize(self, state):
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
        mch_mol, result = self._optimizer.get_result(
            mol=mch_mol,
            bond_pair_ids=long_bond_ids,
            subunits=get_subunits(state),
        )
        return state.with_position_matrix(
            mch_mol.get_position_matrix()
        )
