"""
Periodic Construction Result
============================

"""

from .construction_result import ConstructionResult
from ...periodic_info import PeriodicInfo


class PeriodicConstructionResult(ConstructionResult):
    """
    The result of :meth:`.TopologyGraph.construct` with periodic info.

    """

    __slots__ = [
        '_periodic_info',
    ]

    def __init__(self, construction_state):
        """
        Initialize a :class:`.ConstructionResult`.

        Parameters
        ----------
        construction_state : :class:`.ConstructionState`
            The state from which the result is initialized.

        """

        super().__init__(construction_state)
        self._periodic_info = PeriodicInfo(
            *construction_state.get_lattice_constants()
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
