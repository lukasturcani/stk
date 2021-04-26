"""
Reactions Summary
=================

"""

from typing import NamedTuple

from .atom_batch import _AtomBatch
from .bond_batch import _BondBatch


class _BondId(NamedTuple):
    """
    Identifies a bond in a molecule.

    Attributes
    ----------
    atom1_id
        The id of the first :class:`.Atom` in the bond.

    atom2_id
        The id of the sceond :class:`.Atom` in the bond.

    """

    atom1_id: int
    atom2_id: int


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
        '_deleted_atom_ids',
        '_deleted_bond_ids',
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
        self._deleted_atom_ids = set()
        self._deleted_bond_ids = set()

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

        atom_batch = _AtomBatch(
            atoms=result.get_new_atoms(),
            num_atoms=self._num_atoms,
        )
        self._with_atom_batch(atom_batch)

        bond_batch = _BondBatch(
            bonds=result.get_new_bonds(),
            atom_map=atom_batch.get_atom_map(),
        )
        self._with_bond_batch(bond_batch)

        self._deleted_atom_ids.update(
            atom.get_id() for atom in result.get_deleted_atoms()
        )

        self._deleted_bond_ids.update(
            _BondId(
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
            )
            for bond in result.get_deleted_bonds()
        )

    def _with_atom_batch(self, batch):
        """
        Add a batch of atoms to the summary.

        Parameters
        ----------
        batch : :class:`._AtomBatch`
            A batch of atoms.

        Returns
        -------
        None : :class:`NoneType`

        """

        self._atoms.extend(batch.get_atoms())
        self._atom_infos.extend(batch.get_atom_infos())
        self._positions.extend(batch.get_positions())

    def _with_bond_batch(self, batch):
        """
        Add a batch of bonds to the summary.

        Parameters
        ----------
        batch : :class:`.BondBatch`
            A batch of bonds.

        Returns
        -------
        None : :class:`NoneType`

        """

        self._bonds.extend(batch.get_bonds())
        self._bond_infos.extend(batch.get_bond_infos())

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

    def get_deleted_atom_ids(self):
        """
        Yield the ids of deletable atoms held by the summary.

        Yields
        ------
        :class:`int`
            The id of an atom which should be deleted.

        """

        yield from self._deleted_atom_ids

    def get_deleted_bond_ids(self):
        """
        Yield the atom ids of bonds to be deleted held by the summary.

        Yields
        ------
        :class:`tuple`
            A tuple of the atom ids of the bond which should be
            deleted.

        """

        yield from self._deleted_bond_ids

    def get_positions(self):
        """
        Yield the positions of atoms held by the summary.

        Yields
        ------
        :class:`numpy.ndarray`
            The position of an atom held by the summary.

        """

        yield from self._positions
