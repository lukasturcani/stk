"""
Defines classes which represent bonds.

"""


class Bond:
    """
    Represents an atomic bond.

    Attributes
    ----------
    atom1 : :class:`.Atom`
        The first atom in the bond.

    atom2 : :class:`.Atom`
        The second atom in the bond.

    order : :class:`int`
        The bond order.

    Examples
    --------
    *Adding additional attributes.*

    Each bond can be given additional attributes. For example

    .. code-block:: python

        import stk

        bond1 = stk.Bond(stk.H(0), stk.H(1), 1)
        bond1.attr1 = 12.2
        bond1.attr2 = 'something'

    A :class:`.Bond` can also initialized with additional attributes
    directly

    .. code-block:: python

        bond2 = stk.Bond(
            atom1=stk.C(2),
            atom2=stk.C(3),
            order=1,
            attr1=123,
            attr2='hi'
        )

        bond2.attr1  # Holds 123.
        bond2.attr2  # Holds 'hi'.

    *Printing*

    To print a brief summary of the bond you can run

    .. code-block:: python

        # Prints Bond(C(2), C(3), 1).
        print(bond2)

    To print a complete description of the bond, including additional
    attributes, you can run

    .. code-block:: python

        # Prints Bond(C(2), C(3), 1, attr1=123, attr2='hi')
        print(repr(bond2))

    If the atoms have additional attributes, they will be printed too

    .. code-block:: python

        # Add a custom attribute to an atom.
        bond2.atom1.alpha = 1

        # Prints Bond(C(2, alpha=1), C(3), 1, attr1=123, attr2='hi')
        print(repr(bond2))

    If private attributes are added to the bond, they will not be
    printed

    .. code-block:: python

        bond2._attr3 = 'private'
        # Prints Bond(C(2, alpha=1), C(3), 1, attr1=123, attr2='hi')
        print(repr(bond2))

    """

    def __init__(self, atom1, atom2, order, **kwargs):
        """
        Initialize a :class:`Bond`.

        Attributes
        ----------
        atom1 : :class:`.Atom`
            The first atom in the bond.

        atom2 : :class:`.Atom`
            The second atom in the bond.

        order : :class:`int`
            The bond order.

        **kwargs : :class:`object`
            Additional attributes to be added to the bond.

        """

        self.atom1 = atom1
        self.atom2 = atom2
        self.order = order
        for attr, val in kwargs.items():
            setattr(self, attr, val)

    def clone(self, atom_map=None):
        """
        Return a clone.

        Parameters
        ----------
        atom_map : :class:`dict`, optional
            If the clone should hold different :class:`.Atom`
            instances, then a :class:`dict` should be provided which
            maps atoms in the current :class:`.Bond` to the
            atoms which should be used in the clone. Only atoms which
            need to be remapped need to be present in the `atom_map`.

        Returns
        -------
        :class:`.Bond`
            A clone of the bond.

        """

        if atom_map is None:
            atom_map = {}

        obj = self.__class__.__new__(self.__class__)
        for attr, val in vars(self).items():
            if not attr.startswith('_'):
                setattr(obj, attr, val)
        obj.atom1 = atom_map.get(obj.atom1, obj.atom1)
        obj.atom2 = atom_map.get(obj.atom2, obj.atom2)
        return obj

    def __repr__(self):
        cls_name = self.__class__.__name__
        if isinstance(self.order, float) and self.order.is_integer():
            self.order = int(self.order)

        return (
            f'{cls_name}({self.atom1!r}, {self.atom2!r}, {self.order})'
        )

    def __str__(self):
        cls_name = self.__class__.__name__
        if isinstance(self.order, float) and self.order.is_integer():
            self.order = int(self.order)

        return (
            f'{cls_name}({self.atom1}, {self.atom2}, {self.order})'
        )


class PeriodicBond(Bond):
    """
    Represents a periodic bond.

    Attributes
    ----------
    direction : :class:`list` of :class:`int`
        The directions across which the bond is periodic. For example,
        ``[1, 0, -1]`` means that when going from
        :attr:`~.Bond.atom1` to :attr:`~.Bond.atom2` the bond is
        periodic across the x axis in the positive direction, is not
        periodic across the y axis and is periodic across the z axis
        in the negative direction.

    """

    def __init__(self, atom1, atom2, order, direction):
        """
        Initialize a :class:`PeriodicBond`.

        Parameters
        ----------
        atom1 : :class:`Atom`
            The first atom in the bond.

        atom2 : :class:`Atom`
            The second atom in the bond.

        order : :class:`int`
            The bond order.

        """

        super().__init__(atom1, atom2, order)
        self.direction = direction
