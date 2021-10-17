"""
Reactions Summary
=================

"""

from collections import abc
from typing import NamedTuple
import numpy as np

from .atom_batch import AtomBatch
from .bond_batch import BondBatch
from ......atom import Atom
from ......atom_info import AtomInfo
from ......bond import Bond
from ......bond_info import BondInfo
from ......reactions import ReactionResult


__all__ = (
    'ReactionsSummary',
    'BondId',
)


class BondId(NamedTuple):
    """
    Identifies a bond in a molecule.

    Attributes:

        atom1_id:
            The id of the first :class:`.Atom` in the bond.

        atom2_id:
            The id of the second :class:`.Atom` in the bond.

    """

    atom1_id: int
    atom2_id: int


class ReactionsSummary:
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

    _atoms: list[Atom]
    _atom_infos: list[AtomInfo]
    _positions: list[np.ndarray]
    _bonds: list[Bond]
    _bond_infos: list[BondInfo]
    _deleted_atom_ids: set[int]
    _deleted_bond_ids: set[BondId]

    def __init__(
        self,
        num_atoms: int,
        reaction_results: abc.Iterable[ReactionResult],
    ) -> None:
        """
        Initialize a :class:`.ReactionsSummary` instance.

        Parameters:

            num_atoms:
                The number of atoms in molecule being constructed,
                before this summary is taken into account.

            reaction_results:
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

    def _with_reaction_result(
        self,
        result: ReactionResult,
    ) -> None:
        """
        Add the `result` to the summary.

        Parameters:

            result:
                The result to add to the summary.

        """

        atom_batch = AtomBatch(
            atoms=result.get_new_atoms(),
            num_atoms=self._num_atoms,
        )
        self._with_atom_batch(atom_batch)

        bond_batch = BondBatch(
            bonds=result.get_new_bonds(),
            atom_map=atom_batch.get_atom_map(),
        )
        self._with_bond_batch(bond_batch)

        self._deleted_atom_ids.update(
            atom.get_id() for atom in result.get_deleted_atoms()
        )

        self._deleted_bond_ids.update(
            BondId(
                atom1_id=bond.get_atom1().get_id(),
                atom2_id=bond.get_atom2().get_id(),
            )
            for bond in result.get_deleted_bonds()
        )

    def _with_atom_batch(
        self,
        batch: AtomBatch,
    ) -> None:
        """
        Add a batch of atoms to the summary.

        Parameters:

            batch:
                A batch of atoms.

        """

        self._atoms.extend(batch.get_atoms())
        self._atom_infos.extend(batch.get_atom_infos())
        self._positions.extend(batch.get_positions())

    def _with_bond_batch(
        self,
        batch: BondBatch,
    ) -> None:
        """
        Add a batch of bonds to the summary.

        Parameters:

            batch:
                A batch of bonds.

        """

        self._bonds.extend(batch.get_bonds())
        self._bond_infos.extend(batch.get_bond_infos())

    def get_atoms(self) -> abc.Iterator[Atom]:
        """
        Yield the atoms in the summary.

        Yields:

            An atom.

        """

        yield from self._atoms

    def get_atom_infos(self) -> abc.Iterator[AtomInfo]:
        """
        Yield infos about atoms in the summary.

        Yields:

            Info about an atom.

        """

        yield from self._atom_infos

    def get_bonds(self) -> abc.Iterator[Bond]:
        """
        Yield the bonds in the summary.

        Yields:

            A bond.

        """

        yield from self._bonds

    def get_bond_infos(self) -> abc.Iterator[BondInfo]:
        """
        Yield infos about the bonds in the summary.

        Yields:

            Info about a bond.

        """

        yield from self._bond_infos

    def get_deleted_atom_ids(self) -> abc.Iterator[int]:
        """
        Yield the ids of deletable atoms held by the summary.

        Yields:

            The id of an atom which should be deleted.

        """

        yield from self._deleted_atom_ids

    def get_deleted_bond_ids(self) -> abc.Iterator[BondId]:
        """
        Yield the atom ids of bonds to be deleted held by the summary.

        Yields:

            A tuple of the atom ids of the bond which should be
            deleted.

        """

        yield from self._deleted_bond_ids

    def get_positions(self) -> abc.Iterator[np.ndarray]:
        """
        Yield the positions of atoms held by the summary.

        Yields:

            The position of an atom held by the summary.

        """

        yield from self._positions
