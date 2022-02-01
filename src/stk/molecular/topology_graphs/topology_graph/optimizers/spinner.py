"""
Spinner
=======

"""

import spindry as spd

from ..construction_state import ConstructionState
from .optimizer import Optimizer


class Spinner(Optimizer):
    """
    Performs Monte Carlo optimisation of host-guest complexes [1]_.

    Examples:
        *Structure Optimization*

        Using :class:`.Spinner` will lead to
        :class:`.ConstructedMolecule` structures with better host-guest
        structures. Especially useful for multiple-guest systems and
        removing overlap.

        .. testcode:: structure-optimization

            import stk

            bb1 = stk.BuildingBlock(
                smiles='NCCN',
                functional_groups=[stk.PrimaryAminoFactory()],
            )
            bb2 = stk.BuildingBlock(
                smiles='O=CC(C=O)C=O',
                functional_groups=[stk.AldehydeFactory()],
            )
            guest1 = stk.host_guest.Guest(
                building_block=stk.BuildingBlock('c1ccccc1'),
            )
            guest2 = stk.host_guest.Guest(
                building_block=stk.BuildingBlock('C1CCCCC1'),
            )
            cage = stk.ConstructedMolecule(
                topology_graph=stk.cage.FourPlusSix(
                    building_blocks=(bb1, bb2),
                    optimizer=stk.MCHammer(),
                ),
            )

            complex = stk.ConstructedMolecule(
                topology_graph=stk.host_guest.Complex(
                    host=stk.BuildingBlock.init_from_molecule(cage),
                    guests=(guest1, guest2),
                    optimizer=stk.Spinner(),
                ),
            )

        Optimisation with :mod:`stk` simply collects the final position
        matrix. The optimisation's trajectory can be output using the
        :mod:`SpinDry` implementation if required by the user [1]_.
        This code is entirely nonphysical and is, therefore, completely
        general to any chemistry.

    References:
        .. [1] https://github.com/andrewtarzia/SpinDry

    """

    def __init__(
        self,
        step_size: float = 1.5,
        rotation_step_size: float = 5.,
        num_conformers: int = 50,
        max_attempts: int = 1000,
        nonbond_epsilon: float = 5.,
        beta: float = 2.,
        random_seed: int = 1000,
    ) -> None:
        """
        Initialize an instance of :class:`.Spinner`.

        Parameters:
            step_size: The relative size of the step to take during
                step.

            rotation_step_size: The relative size of the rotation to
                take during step.

            num_conformers: Number of conformers to extract.

            max_attempts: Maximum number of MC moves to try to generate
                conformers.

            nonbond_epsilon: Value of epsilon used in the nonbonded
                potential in MC moves. Determines strength of the
                nonbonded potential.

            beta: Value of beta used in the in MC moves. Beta takes the
                place of the inverse boltzmann temperature.

            random_seed: Random seed to use for MC algorithm.

        """

        self._optimizer = spd.Spinner(
            step_size=step_size,
            rotation_step_size=rotation_step_size,
            num_conformers=num_conformers,
            max_attempts=max_attempts,
            potential_function=spd.SpdPotential(
                nonbond_epsilon=nonbond_epsilon,
            ),
            beta=beta,
            random_seed=random_seed,
        )

    def optimize(self, state: ConstructionState) -> ConstructionState:
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

        conformer = self._optimizer.get_final_conformer(supramolecule)
        return state.with_position_matrix(
            position_matrix=conformer.get_position_matrix(),
        )
