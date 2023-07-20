from dataclasses import dataclass

import rdkit.Chem.AllChem as rdkit


@dataclass(frozen=True, slots=True)
class Atom:
    """
    Represents atoms.

    Examples:

        *Creating atoms*

        You can create atoms by providing an id and an atomic number
        or chemical symbol:

        .. testcode:: creating-atoms

            import stk

            helium = stk.Atom.new(0, 2)
            also_helium = stk.Atom.new(0, 'He')

        Charged atoms can be created by providing an optional charge:

        .. testcode:: creating-atoms

            charged_helium = stk.Atom.new(0, 'He', 1)

        *Getting the element symbol*

        .. doctest:: getting-element

            >>> import stk
            >>> atom = stk.Atom.new(0, 1)
            >>> atom.get_symbol()
            'H'

    Parameters:
        id: The id of the atom.
        atomic_number: The atomic number of the atom.
        charge: The charge of the atom.
    """

    id: int
    """The id of the atom."""
    atomic_number: int
    """The atomic number of the atom."""
    charge: int = 0
    """The charge of the atom."""

    @staticmethod
    def new(id: int, element: int | str, charge: int = 0) -> "Atom":
        """
        Create a new atom.

        Parameters:
            id: The id of the atom.
            atomic_number: The atomic number of the atom.
            charge: The charge of the atom.
        Returns:
            A new atom.
        """
        atomic_number = (
            element
            if isinstance(element, int)
            else rdkit.GetPeriodicTable().GetAtomicNumber(element)
        )
        return Atom(id, atomic_number, charge)

    def get_symbol(self) -> str:
        """
        Get the element symbol of the atom.

        Returns:
            The element symbol.
        """
        return rdkit.GetPeriodicTable().GetElementSymbol(self.atomic_number)
