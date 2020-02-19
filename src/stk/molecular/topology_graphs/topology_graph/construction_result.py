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

    atom_infos : :class:`tuple` of :class:`.AtomInfo`
        Holds additional data for each atom in the constructed
        molecule, for example, which building block it originated
        from.

    reaction_infos : :class:`tuple` of :class:`.ReactionInfo`
        Holds data on the reactions which were performed during
        construction.

    """

    atoms: object
    bonds: object
    position_matrix: object
    atom_infos: object
    reaction_infos: object
