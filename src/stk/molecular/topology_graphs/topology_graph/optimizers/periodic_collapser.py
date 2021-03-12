"""
PeriodicCollapser
=========

"""

import numpy as np
from .collapser import Collapser
from .utilities import get_mch_bonds, get_long_bond_ids, get_subunits

import mchammer as mch


class PeriodicCollapser(Collapser):
    """
    Performs rigid-body collapse of molecule [1]_.

    Examples
    --------
    *Structure Optimization*

    Using :class:`.PeriodicCollapser` will lead to
    :class:`.ConstructedMolecule` structures without long bonds and
    match the unit-cell to the new structure.

    .. code-block:: python

        FILL IN WITH TESTING.py

    Optimisation with :mod:`stk` simply collects the final position
    matrix and periodic info. The optimisation's trajectory can be
    output using the :mod:`MCHammer` implementation if required by the
    user [1]_.

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
        step_size=0.1,
        distance_threshold=1.5,
        scale_steps=True,
    ):
        """
        Initialize an instance of :class:`.PeriodicCollapser`.

        Parameters
        ----------
        step_size : :class:`float`, optional
            The relative size of the step to take during collapse in
            Angstrom.

        distance_threshold : :class:`float`, optional
            Distance between distinct building blocks to use as
            threshold for halting collapse in Angstrom.

        scale_steps : :class:`bool`, optional
            Whether to scale the step of each distinct building block
            by its relative distance from the molecules centroid.

        """

        self._optimizer = mch.Collapser(
            step_size=step_size,
            distance_threshold=distance_threshold,
            scale_steps=scale_steps,
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
            bonds=get_mch_bonds(state),
            position_matrix=state.get_position_matrix(),
        )

        # Run optimization.
        mch_mol, result = self._optimizer.get_result(
            mol=mch_mol,
            bond_pair_ids=tuple(get_long_bond_ids(state)),
            subunits=get_subunits(state),
        )

        print(result)

        old_pos_mat = state.get_position_matrix()
        new_pos_mat = mch_mol.get_position_matrix()
        print(old_pos_mat)
        print(new_pos_mat)
        old_x_extent = abs(
            max(old_pos_mat[:, 0])-min(old_pos_mat[:, 0])
        )
        old_y_extent = abs(
            max(old_pos_mat[:, 1])-min(old_pos_mat[:, 1])
        )
        old_z_extent = abs(
            max(old_pos_mat[:, 2])-min(old_pos_mat[:, 2])
        )
        print(old_x_extent, old_y_extent, old_z_extent)
        new_x_extent = abs(
            max(new_pos_mat[:, 0])-min(new_pos_mat[:, 0])
        )
        new_y_extent = abs(
            max(new_pos_mat[:, 1])-min(new_pos_mat[:, 1])
        )
        new_z_extent = abs(
            max(new_pos_mat[:, 2])-min(new_pos_mat[:, 2])
        )
        print(new_x_extent, new_y_extent, new_z_extent)
        x_ratio = new_x_extent/old_x_extent
        y_ratio = new_y_extent/old_y_extent
        z_ratio = new_z_extent/old_z_extent
        print(x_ratio, y_ratio, z_ratio)
        old_lattice = state.get_lattice_constants()
        print(old_lattice)
        new_lattice = np.array([
            old_lattice[0]*x_ratio,
            old_lattice[1]*y_ratio,
            old_lattice[2]*z_ratio,
        ])
        print(new_lattice)
        state = state.with_lattice_constants(new_lattice)

        return state.with_position_matrix(
            position_matrix=mch_mol.get_position_matrix()
        )
