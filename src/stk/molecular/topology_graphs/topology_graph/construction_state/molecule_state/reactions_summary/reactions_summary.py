from .atom_batch_data import _AtomBatchData
from .bond_batch_data import _BondBatchData


class _ReactionsSummary:
    """
    A summary of reaction results.

    """

    __slots__ = [
        '_num_atoms',
        '_atoms',
        '_atom_infos',
        '_positions',
        '_bonds',
        '_bond_infos',
        '_deleted_ids',
    ]

    def __init__(self, num_atoms, reaction_results):
        """
        Initialize a :class:`.ReactionsSummary` instance.

        Parameters
        ----------
        num_atoms : :class:`int`
            The number of atoms in molecule being constructed,
            before this summary is taken into account.

        reaction_results : :class:`iterable`
            Holds the :class:`.ReactionResult` instances to be
            summarized.

        """

        # This will get updated as reaction results are added to the
        # summary.
        self._num_atoms = num_atoms
        self._atoms = []
        self._atom_infos = []
        self._positions = []
        self._bonds = []
        self._bond_infos = []
        self._deleted_ids = set()

        for result in reaction_results:
            self._with_reaction_result(result)
            self._num_atoms += len(result.get_new_atoms())

    def _with_reaction_result(self, result):
        """
        Add the `result` to the summary.

        Parameters
        ----------
        result : :class:`.ReactionResult`
            The result to add to the summary.

        Returns
        -------
        None : :class:`NoneType`

        """

        atom_batch_data = _AtomBatchData(
            atoms=result.get_new_atoms(),
            num_atoms=self._num_atoms,
        )
        self._with_atom_batch_data(atom_batch_data)

        bond_batch_data = _BondBatchData(
            bonds=result.get_new_bonds(),
            atom_map=atom_batch_data.get_atom_map(),
        )
        self._with_bond_batch_data(bond_batch_data)

        self._deleted_ids.update(result.get_deleted_ids())

    def _with_atom_batch_data(self, batch_data):
        """
        Add data about a batch of atoms to the summary.

        Parameters
        ----------
        batch_data : :class:`._AtomBatchData`
            Data about a batch of atoms.

        Returns
        -------
        None : :class:`NoneType`

        """

        self._atoms.extend(batch_data.get_atoms())
        self._atom_infos.extend(batch_data.get_atom_infos())
        self._positions.extend(batch_data.get_positions())

    def _with_bond_batch_data(self, batch_data):
        """
        Add data about a batch of bonds to the summary.

        Parameters
        ----------
        batch_data : :class:`.BondBatchData`
            Data about a batch of bonds.

        Returns
        -------
        None : :class:`NoneType`

        """

        self._bonds.extend(batch_data.get_bonds())
        self._bond_infos.extend(batch_data.get_bond_infos())

    def get_atoms(self):
        """
        Yield the atoms in the summary.

        Yields
        ------
        :class:`.Atom`
            An atom.

        """

        yield from self._atoms

    def get_atom_infos(self):
        """
        Yield infos about atoms in the summary.

        Yields
        ------
        :class:`.AtomInfo`
            Info about an atom.

        """

        yield from self._atom_infos

    def get_bonds(self):
        """
        Yield the bonds in the summary.

        Yields
        ------
        :class:`.Bond`
            A bond.

        """

        yield from self._bonds

    def get_bond_infos(self):
        """
        Yield infos about the bonds in the summary.

        Yields
        ------
        :class:`.BondInfo`
            Info about a bond.

        """

        yield from self._bond_infos

    def get_deleted_ids(self):
        """
        Yield the ids of deletable atoms held by the summary.

        Yields
        ------
        :class:`int`
            The id of an atom which should be deleted.

        """

        yield from self._deleted_ids

    def get_positions(self):
        """
        Yield the positions of atoms held by the summary.

        Yields
        ------
        :class:`numpy.ndarray`
            The position of an atom held by the summary.

        """

        yield from self._positions
