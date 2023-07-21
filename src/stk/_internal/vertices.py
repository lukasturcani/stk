from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from stk._internal.math import apply_rotation_between_vectors


@dataclass(frozen=True, slots=True)
class PlacementVertex:
    position: npt.NDArray[np.float32]

    def place(
        self,
        matrix: npt.NDArray[np.float32],
        position_anchor: npt.NDArray[np.float32],
    ) -> npt.NDArray[np.float32]:
        matrix = np.array(matrix, dtype=np.float32)
        matrix += self.position - position_anchor
        return matrix


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
        matrix = np.array(matrix, dtype=np.float32)
        matrix += self.position - position_anchor
        apply_rotation_between_vectors(
            matrix=matrix,
            start=rotation_anchor_axis,
            target=self.rotation_axis,
            origin=self.position,
        )
        return matrix


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
        matrix = np.array(matrix, dtype=np.float32)
        matrix += self.position - position_anchor
        return matrix
