"""
MCHammer
========

"""

import mchammer as mch

from .optimizer import Optimizer
from .utilities import get_long_bond_ids, get_mch_bonds, get_subunits


class MCHammer(Optimizer):
    """
    Performs Monte Carlo optimisation of long-bonds in molecules [1]_.

    Examples
    --------
    *Structure Optimization*

    Using :class:`.MCHammer` will lead to :class:`.ConstructedMolecule`
    structures without long bonds.

    .. testcode:: structure-optimization

        import stk

        bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock('O=CCC=O', [stk.AldehydeFactory()])

        polymer = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(bb1, bb2),
                repeating_unit='AB',
                num_repeating_units=6,
                optimizer=stk.MCHammer(),
            ),
        )

    Optimisation with :mod:`stk` simply collects the final position
    matrix. The optimisation's trajectory can be output using the
    :mod:`MCHammer` implementation if required by the user [1]_.

    The open-source optimization code :mod:`MCHammer` specializes in
    the `collapsing` of molecules with long bonds like those
    constructed by :mod:`stk`. This code is entirely nonphysical and
    is, therefore, completely general to any chemistry.

    References
    ----------
    .. [1] https://github.com/andrewtarzia/MCHammer

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
        random_seed=1000,
    ):
        """
        Initialize an instance of :class:`.MCHammer`.

        Parameters
        ----------
        step_size : :class:`float`, optional
            The relative size of the step to take during step.

        target_bond_length : :class:`float`, optional
            Target equilibrium bond length for long bonds to minimize
            to in Angstrom.

        num_steps : :class:`int`, optional
            Number of MC moves to perform.

        bond_epsilon : :class:`float`, optional
            Value of epsilon used in the bond potential in MC moves.
            Determines strength of the bond potential.

        nonbond_epsilon : :class:`float`, optional
            Value of epsilon used in the nonbond potential in MC moves.
            Determines strength of the nonbond potential.
            Larger values lead to a larger building block repulsion.

        nonbond_sigma : :class:`float`, optional
            Value of sigma used in the nonbond potential in MC moves.
            Larger values lead to building block repulsion at larger
            distances.

        nonbond_mu : :class:`float`, optional
            Value of mu used in the nonbond potential in MC moves.
            Determines the steepness of the nonbond potential.

        beta : :class:`float`, optional
            Value of beta used in the in MC moves. Beta takes the
            place of the inverse Boltzmann temperature.

        random_seed : :class:`int` or :class:`NoneType`, optional
            Random seed to use for MC algorithm. If
            ``None`` a system-based random seed will be used
            and results will not be reproducible between
            invocations.

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

    def optimize(self, state):
        # Define MCHammer molecule to optimize.
        mch_mol = mch.Molecule(
            atoms=(
                mch.Atom(
                    id=atom.get_id(),
                    element_string=atom.__class__.__name__,
                ) for atom in state.get_atoms()
            ),
            bonds=tuple(get_mch_bonds(state)),
            position_matrix=state.get_position_matrix(),
        )

        # Run optimization.
        mch_mol, results = self._optimizer.get_result(
            mol=mch_mol,
            bond_pair_ids=tuple(get_long_bond_ids(state)),
            subunits=get_subunits(state),
        )
        return state.with_position_matrix(
            position_matrix=mch_mol.get_position_matrix(),
        )
