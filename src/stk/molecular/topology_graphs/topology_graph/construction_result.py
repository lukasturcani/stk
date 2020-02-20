class ConstructionResult:
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

    building_block_counts : :class:`dict`
        Maps each :class:`.BuildingBlock` used during construction to
        the number of times it was used.

    """

    __slots__ = [
        'atoms',
        'bonds',
        'position_matrix',
        'atom_infos',
        'reaction_infos',
        'building_block_counts',
    ]

    def __init__(self, construction_state):
        """
        Initialize a :class:`.ConstructionResult` instance.

        Parameters
        ----------
        construction_state : :class:`.ConstructionState`
            The construction state from which to create the result.

        """

        self.atoms = construction_state.get_atoms()
        self.bonds = construction_state.get_bonds()
        self.position_matrix = construction_state.get_position_matrix()
        self.atom_infos = construction_state.get_atom_infos()
        self.reaction_infos = construction_state.get_reaction_infos()
        self.building_block_counts = (
            construction_state.get_building_block_counts()
        )
