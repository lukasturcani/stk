from dataclasses import dataclass

import numpy as np
import numpy.typing as npt


@dataclass(frozen=True, slots=True)
class PlacementVertex:
    position: npt.NDArray[np.float32]

    def place(
        self,
        matrix: npt.NDArray[np.float32],
        position_anchor: npt.NDArray[np.float32],
    ) -> npt.NDArray[np.float32]:
        pass


@dataclass(frozen=True, slots=True)
class OrientationVertex:
    position: npt.NDArray[np.float32]
    rotation_axis: npt.NDArray[np.float32]

    def place(
        self,
        matrix: npt.NDArray[np.float32],
        position_anchor: npt.NDArray[np.float32],
        rotation_anchor_axis: npt.NDArray[np.float32],
    ) -> npt.NDArray[np.float32]:
        pass


@dataclass(frozen=True, slots=True)
class RotationVertex:
    position: npt.NDArray[np.float32]
    rotation_axis: npt.NDArray[np.float32]
    rotation_target: npt.NDArray[np.float32]

    def place(
        self,
        matrix: npt.NDArray[np.float32],
        position_anchor: npt.NDArray[np.float32],
        rotation_anchor_axis: npt.NDArray[np.float32],
        rotation_anchor_target: npt.NDArray[np.float32],
    ) -> npt.NDArray[np.float32]:
        pass
