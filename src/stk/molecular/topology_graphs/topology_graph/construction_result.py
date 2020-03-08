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

    bond_infos : :class:`tuple` of :class:`.BondInfo`
        Holds additional data for each bond in the constructed
        molecule, for example, which building block it originated from.

    """

    atoms: object
    bonds: object
    position_matrix: object
    atom_infos: object
    bond_infos: object

    @classmethod
    def init_from_construction_state(cls, construction_state):
        """
        Initialize a :class:`.ConstructionResult` instance.

        Parameters
        ----------
        construction_state : :class:`.ConstructionState`
            The construction state from which to create the result.

        """

        position_matrix = construction_state.get_position_matrix()
        position_matrix.setflags(write=False)
        return cls(
            atoms=tuple(construction_state.get_atoms()),
            bonds=tuple(construction_state.get_bonds()),
            position_matrix=position_matrix,
            atom_infos=tuple(construction_state.get_atom_infos()),
            bond_infos=tuple(
                construction_state.get_bond_infos()
            ),
        )
