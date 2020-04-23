from ......atoms import AtomInfo


class _AtomBatch:
    """
    A batch of atoms.

    """

    __slots__ = ['_atoms', '_atom_infos', '_atom_map']

    def __init__(
        self,
        atoms,
        num_atoms,
        building_block,
        building_block_id,
    ):
        """
        Initialize an :class:`._AtomBatch` instance.

        Parameters
        ----------
        atoms : :class:`iterable` of :class:`.Atom`
            The atoms, which should be added to the batch.

        num_atoms : :class:`int`
            The number of atoms in the molecule being constructed,
            before atoms in this batch are taken into account.

        building_block : :class:`.BuildingBlock`
            The building block from which the atoms originate.

        building_block_id : :class:`.int`
            An id, unique to that building block and placement.

        """

        self._atoms = _atoms = []
        self._atom_infos = atom_infos = []
        self._atom_map = atom_map = {}

        for id_, atom in enumerate(atoms, num_atoms):
            _atoms.append(atom.with_id(id_))
            atom_map[atom.get_id()] = _atoms[-1]
            atom_infos.append(
                AtomInfo(
                    atom=_atoms[-1],
                    building_block=building_block,
                    building_block_id=building_block_id,
                )
            )

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
