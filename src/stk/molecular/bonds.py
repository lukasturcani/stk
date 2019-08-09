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

    periodicity : :class:`tuple` of :class:`int`
        The directions across which the bond is periodic. For example,
        ``(1, 0, -1)`` means that when going from
        :attr:`~.Bond.atom1` to :attr:`~.Bond.atom2` the bond is
        periodic across the x axis in the positive direction, is not
        periodic across the y axis and is periodic across the z axis
        in the negative direction.

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

    To print a brief summary of a bond you can run

    .. code-block:: python

        # Prints Bond(C(2), C(3), 1).
        print(bond2)

    To print a complete description of a bond, including additional
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

    If private attributes were added to a bond, they will not be
    printed

    .. code-block:: python

        bond2._attr3 = 'private'
        # Prints Bond(C(2, alpha=1), C(3), 1, attr1=123, attr2='hi')
        print(repr(bond2))

    """

    def __init__(
        self,
        atom1,
        atom2,
        order,
        periodicity=(0, 0, 0),
        **kwargs
    ):
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

        periodicity : :class:`tuple` of :class:`int`, optional
            The directions across which the bond is periodic. For
            example, ``(1, 0, -1)`` means that when going from
            :attr:`~.Bond.atom1` to :attr:`~.Bond.atom2` the bond is
            periodic across the x axis in the positive direction, is
            not periodic across the y axis and is periodic across the z
            axis in the negative direction.

        **kwargs : :class:`object`, optional
            Additional attributes to be added to the bond.

        """

        self.atom1 = atom1
        self.atom2 = atom2
        self.order = order
        self.periodicity = periodicity
        for attr, val in kwargs.items():
            setattr(self, attr, val)

    def clone(self, atom_map=None):
        """
        Return a clone.

        Private attributes are not passed to the clone.

        Parameters
        ----------
        atom_map : :class:`dict`, optional
            If the clone should hold different :class:`.Atom`
            instances, then a :class:`dict` should be provided, which
            maps atoms in the current :class:`.Bond` to the
            atoms which should be used in the clone. Only atoms which
            need to be remapped need to be present in the `atom_map`.

        Returns
        -------
        :class:`.Bond`
            The clone.

        Examples
        --------
        .. code-block:: python

            import stk

            c0 = stk.C(0)
            c1 = stk.C(1)
            bond = stk.Bond(c0, c1, 1, custom_attr=12, _private_attr=1)

            # bond_clone holds c0 and c1 in its atom1 and atom2
            # attributes, respectively. It also has a custom_attr
            # with a value of 12 but it does not have a _private_attr
            # attribute.
            bond_clone = bond.clone()

        It is possible to make sure that the clone holds different
        atoms

        .. code-block:: python

            li2 = stk.Li(2)
            n3 = stk.N(3)

            # clone2 is also a clone, except that it holds
            # li2 in the atom2 attribute. Its atom1 attribute still
            # holds c0.
            clone2 = bond.clone(atom_map={
                c1: li2
            })

            # clone3 is also a clone, except that it holds n3 and
            # li2 in its atom1 and atom2 attributes, respectively.
            clone3 = bond.clone(atom_map={
                c0: n3,
                c1: li2
            })

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

    def is_periodic(self):
        """
        Return ``True`` if the bond is periodic.

        Returns
        -------
        :class:`bool`
            ``True`` if the bond is periodic.

        """

        return any(direction != 0 for direction in self.periodicity)

    def __repr__(self):
        if isinstance(self.order, float) and self.order.is_integer():
            self.order = int(self.order)

        periodicity = (
            '' if not self.is_periodic()
            else f', periodicity={self.periodicity}'
        )

        mandatory = {'atom1', 'atom2', 'order', 'periodicity'}
        attrs = ', '.join(
            f'{attr}={val!r}' for attr, val in vars(self).items()
            if attr not in mandatory and not attr.startswith('_')
        )
        cls_name = self.__class__.__name__
        return (
            f'{cls_name}({self.atom1!r}, {self.atom2!r}, '
            f'{self.order}{periodicity}{", " if attrs else ""}{attrs})'
        )

    def __str__(self):
        cls_name = self.__class__.__name__
        if isinstance(self.order, float) and self.order.is_integer():
            self.order = int(self.order)

        periodicity = (
            '' if not self.is_periodic()
            else f', {self.periodicity}'
        )
        return (
            f'{cls_name}({self.atom1}, {self.atom2}, '
            f'{self.order}{periodicity})'
        )
