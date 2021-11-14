"""
Periodic Construction Result
============================

"""

from ....periodic_info import PeriodicInfo
from .construction_result import ConstructionResult


class PeriodicConstructionResult(ConstructionResult):
    """
    The result of :meth:`.TopologyGraph.construct` with periodic info.

    """

    __slots__ = [
        '_periodic_info',
    ]

    def __init__(
        self,
        construction_state,
        lattice_size,
    ):
        """
        Initialize a :class:`.ConstructionResult`.

        Parameters
        ----------
        construction_state : :class:`.ConstructionState`
            The state from which the result is initialized.

        lattice_size : :class:`tuple` of :class:`int`
            The size of the lattice in the x, y and z directions.

        """

        super().__init__(construction_state)
        self._periodic_info = PeriodicInfo(
            *(
                lattice_constant*dim
                for lattice_constant, dim in zip(
                    construction_state.get_lattice_constants(),
                    lattice_size,
                )
            )
        )

    def get_periodic_info(self):
        """
        Get the periodic cell information of the constructed molecule.

        Returns
        -------
        :class:`.PeriodicInfo`
            The periodic cell information of the constructed molecule.

        """

        return self._periodic_info
