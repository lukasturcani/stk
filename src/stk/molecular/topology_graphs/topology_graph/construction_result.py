from typing import NamedTuple


class ConstructionResult(NamedTuple):
    """
    The result of a :class:`.ConstructedMolecule` construction.

    Attributes
    ----------
    atoms : :class:`tuple` of :class:`.Atom`
        The atoms of the constructed molecule.

    bonds : :class:`tuple` of :class:`.Bond`
        The bonds of the constructed molecule.

    position_matrix : :class:`numpy.ndarray`
        The position matrix of the constructed molecule.

    """

    atoms: object
    bonds: object
    position_matrix: object
