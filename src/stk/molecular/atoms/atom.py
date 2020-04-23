"""
Atom
====

"""


class Atom:
    """
    An abstract base class for atoms.

    A subclass is made for each element. The name of each subclass is
    the periodic table symbol of that element.

    Atoms of a particular element can be made with this
    class or with the subclass representing that element.

    Examples
    --------
    *Initialization.*

    Initialization of an :class:`.Atom` can happen in one of two ways.
    The atom can be initialized through the :class:`.Atom` class or
    through the class representing the element.

    .. code-block:: python

        import stk

        # h0 is an instance of the H class.
        h0 = stk.Atom(id=0, atomic_number=1)

        # h1 is also an instance of the H class.
        h1 = stk.H(id=1)

    When the class correspnding to the element is used directly, the
    ``atomic_number`` is not provided. Here are a few more examples.

    .. code-block:: python

        # Both he0 and he1 are instances of the He class.
        he0 = stk.Atom(id=2, atomic_number=2)
        he1 = stk.He(id=3)

        # Both c0 and c1 are instances of the
        # C class.
        c0 = stk.Atom(id=4, atomic_number=6)
        c1 = stk.C(id=5)

    """

    # Maps each atomic number (int) to the relevant Atom subclass.
    _elements = {}

    def __init_subclass__(cls, **kwargs):
        # Replace the default __init__() method of the subclass with
        # _subclass_init(). This is because the default __init__()
        # method takes an atomic_number parameter, but
        # _subclass_init() does not.
        cls.__init__ = cls._subclass_init
        cls._elements[cls._atomic_number] = cls

    @staticmethod
    def _subclass_init(self, id, charge=0):
        """
        Initialize an atom of the element.

        Parameters
        ----------
        id : :class:`int`
            The id of the atom.

        charge : :class:`int`
            The formal charge.

        """

        Atom.__init__(self, id, self._atomic_number, charge)

    def __init__(self, id, atomic_number, charge=0):
        """
        Initialize an :class:`Atom`.

        Parameters
        ----------
        id : :class:`int`
            The id of the atom.

        atomic_number : :class:`int`
            The atomic number.

        charge : :class:`int`
            The formal charge.

        """

        self.__class__ = self._elements[atomic_number]
        self._id = id
        self._charge = charge

    def get_id(self):
        """
        Get the id of the atom.

        Returns
        -------
        :class:`int`
            The id.

        """

        return self._id

    def _with_id(self, id):
        """
        Modify the atom.

        """

        self._id = id
        return self

    def with_id(self, id):
        """
        Get a clone but with a different id.

        Returns
        -------
        :class:`.Atom`
            A clone with a new id. Has the same type as the original
            atom.

        """

        return self.clone()._with_id(id)

    def get_atomic_number(self):
        """
        Get the atomic number of the atom.

        Returns
        -------
        :class:`int`
            The atomic number.

        """

        return self._atomic_number

    def get_charge(self):
        """
        Get the charge of the atom.

        Returns
        -------
        :class:`int`
            The charge.

        """

        return self._charge

    def __repr__(self):
        charge = (
            f', charge={self._charge}' if self._charge != 0 else ''
        )
        return f'{self.__class__.__name__}({self._id}{charge})'

    def __str__(self):
        return repr(self)

    def clone(self):
        """
        Return a clone.

        Returns
        -------
        :class:`.Atom`
            The clone. It has the same type as the original atom.

        """

        clone = self.__class__.__new__(self.__class__)
        clone._id = self._id
        clone._charge = self._charge
        return clone
