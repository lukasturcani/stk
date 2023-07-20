from collections.abc import Sequence
from typing import Any, TypeAlias

import numpy as np
import numpy.typing as npt

Array: TypeAlias = npt.NDArray[Any] | Sequence[npt.NDArray[Any]]


def get_centroid(matrix: Array) -> npt.NDArray[np.float32]:
    return np.sum(matrix, axis=0, dtype=np.float32) / len(matrix)


def get_plane_normal(points: Array) -> npt.NDArray[np.float32]:
    return np.array(
        np.around(np.linalg.svd(points - get_centroid(points))[-1][2, :], 14),
        dtype=np.float32,
    )


def get_acute_vector(
    reference: npt.NDArray[np.float32],
    vector: npt.NDArray[np.float32],
) -> npt.NDArray[np.float32]:
    if (
        # vector_angle is NaN if reference is [0, 0, 0].
        not np.allclose(reference, [0, 0, 0], atol=1e-5)
        and get_vector_angle(vector, reference) > np.pi / 2
    ):
        return vector * -1
    return vector


def get_vector_angle(
    vector1: npt.NDArray[np.float32],
    vector2: npt.NDArray[np.float32],
) -> float:
    if np.all(np.equal(vector1, vector2)):
        return 0.0

    numerator = np.dot(vector1, vector2)
    denominator = np.linalg.norm(vector1) * np.linalg.norm(vector2)
    # This if statement prevents returns of NaN due to floating point
    # inaccuracy.
    term = numerator / denominator
    if term >= 1.0:
        return 0.0
    if term <= -1.0:
        return np.pi
    return np.arccos(term)


def get_orthogonal_vector(vector):
    ortho = np.array([0.0, 0.0, 0.0])
    for m, val in enumerate(vector):
        if not np.allclose(val, 0, atol=1e-8):
            n = (m + 1) % 3
            break
    ortho[n] = vector[m]
    ortho[m] = -vector[n]
    return ortho
