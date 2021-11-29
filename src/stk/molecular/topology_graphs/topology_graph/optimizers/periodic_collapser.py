"""
Periodic Collapser
==================

"""

import mchammer as mch

from .optimizer import Optimizer
from .utilities import get_long_bond_ids, get_mch_bonds, get_subunits


class PeriodicCollapser(Optimizer):
    """
    Performs rigid-body collapse of molecules [1]_.

    This :class:`.Optimizer` will also update the `.PeriodicInfo`.

    Examples
    --------
    *Structure Optimization*

    Using :class:`.PeriodicCollapser` will lead to
    :class:`.ConstructedMolecule` structures without long bonds and
    match the unit-cell to the new structure.

    .. testcode:: structure-optimization

        import stk

        bb1 = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])
        bb2 = stk.BuildingBlock('BrCC(CBr)CBr', [stk.BromoFactory()])

        topology_graph = stk.cof.PeriodicHoneycomb(
            building_blocks=(bb1, bb2),
            lattice_size=(1, 2, 3),
            optimizer=stk.PeriodicCollapser(),
        )
        cof = stk.ConstructedMolecule(topology_graph)

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
        scale_steps=False,
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

        mch_mol, result = self._optimizer.get_result(
            mol=mch_mol,
            bond_pair_ids=tuple(get_long_bond_ids(state)),
            subunits=get_subunits(state),
        )

        old_pos_mat = state.get_position_matrix()
        new_pos_mat = mch_mol.get_position_matrix()
        old_extents = (
            abs(max(old_pos_mat[:, i])-min(old_pos_mat[:, i]))
            for i in range(3)
        )
        new_extents = (
            abs(max(new_pos_mat[:, i])-min(new_pos_mat[:, i]))
            for i in range(3)
        )
        ratios = (n/o for n, o in zip(new_extents, old_extents))
        old_lattice = state.get_lattice_constants()
        new_lattice = tuple(
            old_lattice[i]*ratio for i, ratio in enumerate(ratios)
        )
        state = state.with_lattice_constants(new_lattice)
        return state.with_position_matrix(
            position_matrix=mch_mol.get_position_matrix()
        )
