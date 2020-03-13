from ......atoms import AtomInfo


class _AtomBatchData:
    """
    Holds data about a batch of atoms.

    """

    __slots__ = ['_atoms', '_atom_infos', '_atom_map', '_positions']

    def __init__(self, atoms, num_atoms):
        """
        Initialize a :class:`._AtomBatchData` instance.

        Parameters
        ----------
        atoms : :class:`iterable` of :class:`.NewAtom`
            The atoms, for which data should be created.

        num_atoms : :class:`int`
            The number of atoms in the molecule being constructed,
            before atoms in this batch are taken into account.

        """

        self._atoms = _atoms = []
        self._positions = positions = []
        self._atom_infos = atom_infos = []
        self._atom_map = atom_map = {}

        for id_, (atom, position) in enumerate(atoms, num_atoms):
            _atoms.append(atom.with_id(id_))
            atom_map[atom.get_id()] = _atoms[-1]
            atom_infos.append(AtomInfo(_atoms[-1], None, None))
            positions.append(position)

    def get_positions(self):
        """
        Yield the positions of atoms in the batch.

        Yields
        ------
        :class:`numpy.ndarray`
            The position of an atom in the batch.

        """

        yield from self._positions

    def get_atoms(self):
        """
        Yield the atoms in the batch.

        Yields
        ------
        :class:`.Atom`
            An atom.

        """

        yield from self._atoms

    def get_atom_infos(self):
        """
        Yield info about the atoms in the batch.

        Yields
        ------
        :class:`.AtomInfo`
            Info about an atom in the batch.

        """

        yield from self._atom_infos

    def get_atom_map(self):
        """
        Get a mapping from the old atom id to the new atom.

        Returns
        -------
        :class:`dict`
            Maps the id of an atom provided to the initializer, to
            the new atom held by the batch, which has an updated id.

        """

        return dict(self._atom_map)
