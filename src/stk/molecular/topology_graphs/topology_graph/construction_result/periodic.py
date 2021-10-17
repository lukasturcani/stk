"""
Periodic Construction Result
============================

"""

from ..construction_state import ConstructionState
from .construction_result import ConstructionResult
from ....periodic_info import PeriodicInfo


__all__ = (
    'PeriodicConstructionResult',
)


class PeriodicConstructionResult(ConstructionResult):
    """
    The result of :meth:`.TopologyGraph.construct` with periodic info.

    """

    __slots__ = [
        '_periodic_info',
    ]

    def __init__(
        self,
        construction_state: ConstructionState,
        lattice_size: tuple[int, int, int],
    ) -> None:
        """
        Initialize a :class:`.ConstructionResult`.

        Parameters:

            construction_state:
                The state from which the result is initialized.

            lattice_size:
                The size of the lattice in the x, y and z directions.

        """

        super().__init__(construction_state)
        lattice_constants = construction_state.get_lattice_constants()
        assert lattice_constants is not None
        self._periodic_info = PeriodicInfo(
            *(
                lattice_constant*dim
                for lattice_constant, dim in zip(
                    lattice_constants,
                    lattice_size,
                )
            )
        )

    def get_periodic_info(self) -> PeriodicInfo:
        """
        Get the periodic cell information of the constructed molecule.

        Returns:

            The periodic cell information of the constructed molecule.

        """

        return self._periodic_info
