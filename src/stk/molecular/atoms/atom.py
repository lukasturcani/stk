"""
Atom
====

"""

from __future__ import annotations


class Atom:
    """
    An abstract base class for atoms.

    A subclass is made for each element. The name of each subclass is
    the periodic table symbol of that element.

    Atoms of a particular element can be made with this
    class or with the subclass representing that element.

    Examples:

        *Initialization.*

        Initialization of an :class:`.Atom` can happen in one of two
        ways. The atom can be initialized through the :class:`.Atom`
        class or through the class representing the element.

        .. testcode:: initialization

            import stk

            # h0 is an instance of the H class.
            h0 = stk.Atom(id=0, atomic_number=1)

            # h1 is also an instance of the H class.
            h1 = stk.H(id=1)

        When the class corresponding to the element is used directly,
        the ``atomic_number`` is not provided. Here are a few more
        examples.

        .. testcode:: initialization

            # Both he0 and he1 are instances of the He class.
            he0 = stk.Atom(id=2, atomic_number=2)
            he1 = stk.He(id=3)

            # Both c0 and c1 are instances of the
            # C class.
            c0 = stk.Atom(id=4, atomic_number=6)
            c1 = stk.C(id=5)

    """

    # Maps each atomic number (int) to the relevant Atom subclass.
    _elements: dict[int, type[Atom]] = {}

    def __init__(
        self,
        id: int,
        atomic_number: int,
        charge: int = 0,
    ) -> None:
        """
        Initialize an :class:`Atom`.

        Parameters:

            id: The id of the atom.

            atomic_number: The atomic number.

            charge: The formal charge.

        """

        self._elements[atomic_number].__init__(self, id, charge)
        self.__class__ = self._elements[atomic_number]

    def get_id(self) -> int:
        """
        Get the id of the atom.

        Returns:
            The id.

        """

        raise NotImplementedError()

    def with_id(self, id: int) -> Atom:
        """
        Get a clone but with a different id.

        Parameters:
            id: The id of the clone.

        Returns:
            A clone with a new id.

        """

        raise NotImplementedError()

    def get_atomic_number(self) -> int:
        """
        Get the atomic number of the atom.

        Returns:
            The atomic number.

        """

        raise NotImplementedError()

    def get_charge(self) -> int:
        """
        Get the charge of the atom.

        Returns:
            The charge.

        """

        raise NotImplementedError()

    def clone(self) -> Atom:
        """
        Return a clone.

        Returns:
            The clone.

        """

        raise NotImplementedError()
