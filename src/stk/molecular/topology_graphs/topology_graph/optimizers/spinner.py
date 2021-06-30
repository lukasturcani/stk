"""
Spinner
=======

"""

from .optimizer import Optimizer

import spindry as spd


class Spinner(Optimizer):
    """
    Performs Monte Carlo optimisation of host-guest complexes [1]_.

    Examples
    --------
    *Structure Optimization*

    Using :class:`.Spinner` will lead to :class:`.ConstructedMolecule`
    structures with better host-guest structures. Especially useful
    for multiple-guest systems and removing overlap.

    .. testcode:: structure-optimization

        import stk

        bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock(
            smiles='O=CC(C=O)C=O',
            functional_groups=[stk.AldehydeFactory()],
        )
        guest1 = stk.host_guest.Guest(stk.BuildingBlock('c1ccccc1'))
        guest2 = stk.host_guest.Guest(stk.BuildingBlock('C1CCCCC1'))
        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusSix(
                building_blocks=(bb1, bb2),
                optimizer=stk.MCHammer(),
            ),
        )
        complex = stk.ConstructedMolecule(
            topology_graph=stk.host_guest.Complex(
                host=stk.BuildingBlock.init_from_molecule(host),
                guests=(guest1, guest2),
                optimizer=stk.Spinner(),
            ),
        )

    Optimisation with :mod:`stk` simply collects the final position
    matrix. The optimisation's trajectory can be output using the
    :mod:`SpinDry` implementation if required by the user [1]_.
    This code is entirely nonphysical and is, therefore, completely
    general to any chemistry.

    References
    ----------
    .. [1] https://github.com/andrewtarzia/SpinDry

    """

    def __init__(
        self,
        step_size=1.5,
        rotation_step_size=5,
        num_conformers=50,
        max_attempts=1000,
        nonbond_epsilon=5,
        nonbond_sigma=1.2,
        beta=2,
        random_seed=1000,
    ):
        """
        Initialize an instance of :class:`.Spinner`.

        Parameters
        ----------
        step_size : :class:`float`
            The relative size of the step to take during step.

        rotation_step_size : :class:`float`
            The relative size of the rotation to take during step.

        num_conformers : :class:`int`
            Number of conformers to extract.

        max_attempts : :class:`int`
            Maximum number of MC moves to try to generate conformers.

        nonbond_epsilon : :class:`float`, optional
            Value of epsilon used in the nonbond potential in MC moves.
            Determines strength of the nonbond potential.
            Defaults to 20.

        nonbond_sigma : :class:`float`, optional
            Value of sigma used in the nonbond potential in MC moves.
            Defaults to 1.2.

        beta : :class:`float`, optional
            Value of beta used in the in MC moves. Beta takes the
            place of the inverse boltzmann temperature.
            Defaults to 2.

        random_seed : :class:`int` or :class:`NoneType`, optional
            Random seed to use for MC algorithm. Should only be set to
            ``None`` if system-based random seed is desired. Defaults
            to 1000.

        """

        self._optimizer = spd.Spinner(
            step_size=step_size,
            rotation_step_size=rotation_step_size,
            num_conformers=num_conformers,
            max_attempts=max_attempts,
            nonbond_epsilon=nonbond_epsilon,
            nonbond_sigma=nonbond_sigma,
            beta=beta,
            random_seed=random_seed,
        )

    def optimize(self, state):
        # Define SpinDry supramolecule to optimize.
        supramolecule = spd.SupraMolecule(
            atoms=(
                spd.Atom(
                    id=atom.get_id(),
                    element_string=atom.__class__.__name__,
                ) for atom in state.get_atoms()
            ),
            bonds=(
                spd.Bond(
                    id=i,
                    atom_ids=(
                        bond.get_atom1().get_id(),
                        bond.get_atom2().get_id(),
                    )
                ) for i, bond in enumerate(state.get_bonds())
            ),
            position_matrix=state.get_position_matrix(),
        )

        # Run optimization.
        conformer = self._optimizer.get_final_conformer(supramolecule)
        return state.with_position_matrix(
            position_matrix=conformer.get_position_matrix(),
        )
