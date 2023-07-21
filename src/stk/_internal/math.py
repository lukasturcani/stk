from collections.abc import Sequence
from typing import Any, TypeAlias

import numpy as np
import numpy.typing as npt
from scipy.spatial.transform import Rotation

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


def get_orthogonal_vector(
    vector: npt.NDArray[np.float32],
) -> npt.NDArray[np.float32]:
    ortho = np.array([0.0, 0.0, 0.0])
    for m, val in enumerate(vector):
        if not np.allclose(val, 0, atol=1e-8):
            n = (m + 1) % 3
            break
    ortho[n] = vector[m]
    ortho[m] = -vector[n]
    return ortho


def apply_rotation_between_vectors(
    matrix: npt.NDArray[np.float32],
    start: npt.NDArray[np.float32],
    target: npt.NDArray[np.float32],
    origin: npt.NDArray[np.float32],
) -> None:
    matrix += -origin
    rotation_matrix = get_rotation_matrix(start, target)
    np.matmul(rotation_matrix, matrix.T, out=matrix.T)
    matrix += origin


def get_rotation_matrix(
    vector1: npt.NDArray[np.float32],
    vector2: npt.NDArray[np.float32],
) -> npt.NDArray[np.float32]:
    """
    Returns a rotation matrix which transforms `vector1` to `vector2`.

    Multiplying `vector1` by the rotation matrix returned by this
    function yields `vector2`.

    Parameters:
        vector1:
            The vector which needs to be transformed to `vector2`.
        vector2:
            The vector onto which `vector1` needs to be transformed.
    Returns:
        A rotation matrix which transforms `vector1` to `vector2`.

    References:
        http://tinyurl.com/kybj9ox
        http://tinyurl.com/gn6e8mz

    """

    # Make sure both inputs are unit vectors.
    vector1 = normalize_vector(vector1)
    vector2 = normalize_vector(vector2)

    # Hande the case where the input and output vectors are equal.
    if np.allclose(vector1, vector2, atol=1e-8):
        return np.identity(3, dtype=np.float32)

    # Handle the case where the rotation is 180 degrees.
    if np.allclose(vector1, np.multiply(vector2, -1), atol=1e-8):
        return get_rotation_matrix_arbitrary_axis(
            angle=np.pi,
            axis=get_orthogonal_vector(vector1),
        )

    v = np.cross(vector1, vector2)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    s = np.linalg.norm(v)
    c = np.dot(vector1, vector2)
    i = np.identity(3)
    mult_factor = (1 - c) / np.square(s)

    # Initialize as a scipy Rotation object, which normalizes the
    # matrix and allows for returns as quaternion or alternative
    # type in the future.
    return Rotation.from_matrix(
        i + vx + np.multiply(np.dot(vx, vx), mult_factor)
    ).as_matrix()


def normalize_vector(
    vector: npt.NDArray[np.float32],
) -> npt.NDArray[np.float32]:
    """
    Normalizes the given vector.

    A new vector is returned, the original vector is not modified.

    Parameters:
        vector:
            The vector to be normalized.
    Returns:
        The normalized vector.

    """

    return np.divide(vector, np.linalg.norm(vector))


def get_rotation_matrix_arbitrary_axis(
    angle: float,
    axis: npt.NDArray[np.float32],
) -> npt.NDArray[np.float32]:
    """
    Returns a rotation matrix of `angle` radians about `axis`.

    Parameters:
        angle:
            The size of the rotation in radians.
        axis:
            A 3 element aray which represents a vector. The vector is the
            axis about which the rotation is carried out. Must be of
            unit magnitude.
    Returns:
        A ``3x3`` array representing a rotation matrix.
    """

    a = np.cos(angle / 2)
    b, c, d = axis * np.sin(angle / 2)

    e11 = np.square(a) + np.square(b) - np.square(c) - np.square(d)
    e12 = 2 * (b * c - a * d)
    e13 = 2 * (b * d + a * c)

    e21 = 2 * (b * c + a * d)
    e22 = np.square(a) + np.square(c) - np.square(b) - np.square(d)
    e23 = 2 * (c * d - a * b)

    e31 = 2 * (b * d - a * c)
    e32 = 2 * (c * d + a * b)
    e33 = np.square(a) + np.square(d) - np.square(b) - np.square(c)

    # Initialize as a scipy Rotation object, which normalizes the
    # matrix and allows for returns as quaternion or alternative
    # type in the future.
    return Rotation.from_matrix(
        np.array([[e11, e12, e13], [e21, e22, e23], [e31, e32, e33]])
    ).as_matrix()
