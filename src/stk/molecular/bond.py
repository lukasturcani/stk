"""
Defines :class:`Bond`.

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

    """

    def __init__(self, atom1, atom2, order, periodicity=(0, 0, 0)):
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

        """

        self.atom1 = atom1
        self.atom2 = atom2
        self.order = order
        self.periodicity = periodicity

    def clone(self, atom_map=None):
        """
        Return a clone.

        Private attributes are not passed to the clone.

        Parameters
        ----------
        atom_map : :class:`dict`, optional
            If the clone should hold specific :class:`.Atom`
            instances, then a :class:`dict` should be provided, which
            maps atom ids of atoms in the current :class:`.Bond` to the
            :class:`.Atom` which the clone should hold.

        Returns
        -------
        :class:`.Bond`
            The clone.

        Examples
        --------
        .. code-block:: python

            import stk

            c0 = stk.C(0)
            c5 = stk.C(5)
            bond = stk.Bond(c0, c5, 1)

            # bond_clone holds clones of c0 and c5 in its atom1 and
            # atom2 attributes, respectively.
            bond_clone = bond.clone()

        It is possible to make sure that the clone holds different
        atoms

        .. code-block:: python

            li2 = stk.Li(2)
            n3 = stk.N(3)

            # clone2 is also a clone, except that it holds
            # li2 in the atom2 attribute. Its atom1 attribute holds a
            # clone of c0.
            clone2 = bond.clone(atom_map={
                c5.id: li2,
            })

            # clone3 is also a clone, except that it holds n3 and
            # li2 in its atom1 and atom2 attributes, respectively.
            clone3 = bond.clone(atom_map={
                c0.id: n3,
                c5.id: li2,
            })

        """

        if atom_map is None:
            atom_map = {}

        clone = self.__class__.__new__(self.__class__)
        for attr, val in vars(self).items():
            if not attr.startswith('_'):
                setattr(clone, attr, val)
        clone.atom1 = atom_map.get(self.atom1.id, self.atom1.clone())
        clone.atom2 = atom_map.get(self.atom2.id, self.atom2.clone())
        return clone

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
            f', periodicity={self.periodicity}'
            if self.is_periodic()
            else ''
        )

        cls_name = self.__class__.__name__
        return (
            f'{cls_name}({self.atom1!r}, {self.atom2!r}, '
            f'{self.order}{periodicity})'
        )

    def __str__(self):
        return repr(self)
