"""
Bond
====

"""

from __future__ import annotations

from typing import TypeVar

from ..atoms import Atom

_T = TypeVar('_T', bound='Bond')


class Bond:
    """
    Represents an atomic bond.

    Examples:

        *Changing the Atoms of a Bond*

        You want to substitute the atoms in the bond for other atoms.
        You can do this by using :meth:`.with_atoms` to create a clone
        of the bond, which holds the replacement atoms

        .. testcode:: changing-the-atoms-of-a-bond

            import stk
            bond = stk.Bond(stk.C(0), stk.C(12), 2)
            # Replace C(0) with H(13) but keep C(12).
            clone = bond.with_atoms({0: stk.H(13)})

        .. testcode:: changing-the-atoms-of-a-bond
            :hide:

            assert clone.get_atom1().get_id() == 13
            assert isinstance(clone.get_atom1(), stk.H)
            assert clone.get_atom2().get_id() == 12
            assert isinstance(clone.get_atom2(), stk.C)

    """

    def __init__(
        self,
        atom1: Atom,
        atom2: Atom,
        order: int,
        periodicity: tuple[int, int, int] = (0, 0, 0),
    ) -> None:
        """
        Initialize a :class:`Bond`.

        Parameters
        ----------
        atom1:
            The first atom in the bond.

        atom2:
            The second atom in the bond.

        order:
            The bond order.

        periodicity:

            The directions across which the bond is periodic. For
            example, ``(1, 0, -1)`` means that when going from
            `atom1` to `atom2` the bond is
            periodic across the x axis in the positive direction, is
            not periodic across the y axis and is periodic across the z
            axis in the negative direction.

        """

        self._atom1 = atom1
        self._atom2 = atom2
        self._order = order
        self._periodicity = periodicity

    def get_atom1(self) -> Atom:
        """
        Get the first atom of the bond.

        Returns:

            The first atom of the bond.

        """

        return self._atom1

    def get_atom2(self) -> Atom:
        """
        Get the second atom of the bond.

        Returns:

            The second atom of the bond.

        """

        return self._atom2

    def get_order(self) -> int:
        """
        Get the bond order of the bond.

        Returns:

            The bond order.

        """

        return self._order

    def get_periodicity(self) -> tuple[int, int, int]:
        """
        Get the periodicity of the bond.

        Returns:

            The directions across which the bond is periodic. For
            example, ``(1, 0, -1)`` means that when going from
            `atom1` to `atom2` the bond is
            periodic across the x axis in the positive direction, is
            not periodic across the y axis and is periodic across the z
            axis in the negative direction.

        """

        return self._periodicity

    def clone(self) -> Bond:
        """
        Return a clone.

        Returns:

            The clone.

        """

        clone = self.__class__.__new__(self.__class__)
        Bond.__init__(
            self=clone,
            atom1=self._atom1,
            atom2=self._atom2,
            order=self._order,
            periodicity=self._periodicity,
        )
        return clone

    def _with_atoms(self: _T, atom_map: dict[int, Atom]) -> _T:
        """
        Modify the bond.

        """

        self._atom1 = atom_map.get(self._atom1.get_id(), self._atom1)
        self._atom2 = atom_map.get(self._atom2.get_id(), self._atom2)
        return self

    def with_atoms(self, atom_map: dict[int, Atom]) -> Bond:
        """
        Return a clone holding different atoms.

        Parameters:

            atom_map:
                Maps the id of an atom in the bond to the new atom the
                clone should hold. If the id of an atom in the bond is
                not found in `atom_map`, the atom will not be replaced
                in the clone.

        Returns:

            The clone.

        """

        return self.clone()._with_atoms(atom_map)

    def with_ids(self, id_map: dict[int, int]) -> Bond:
        """
        Return a clone holding different atom ids.

        Parameters:

            id_map:
                Maps the id of an atom in the bond to the new id the
                clone should hold. If the id of an atom in the bond is
                not found in `id_map`, the atom id will not be replaced
                in the clone.

        Returns:

            The clone.

        """

        return self.clone()._with_ids(id_map)

    def _with_ids(self, id_map: dict[int, int]) -> Bond:

        id1 = self._atom1.get_id()
        if id1 in id_map:
            self._atom1 = self._atom1.with_id(id_map[id1])

        id2 = self._atom2.get_id()
        if id2 in id_map:
            self._atom2 = self._atom2.with_id(id_map[id2])

        return self

    def is_periodic(self) -> bool:
        """
        Return ``True`` if the bond is periodic.

        Returns:

            ``True`` if the bond is periodic.

        """

        return any(direction != 0 for direction in self._periodicity)

    def __repr__(self) -> str:
        periodicity = (
            f', periodicity={self._periodicity}'
            if self.is_periodic()
            else ''
        )

        cls_name = self.__class__.__name__
        return (
            f'{cls_name}({self._atom1!r}, {self._atom2!r}, '
            f'{self._order}{periodicity})'
        )

    def __str__(self) -> str:
        return repr(self)
