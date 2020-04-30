"""
Bond
====

"""


class Bond:
    """
    Represents an atomic bond.

    """

    def __init__(self, atom1, atom2, order, periodicity=(0, 0, 0)):
        """
        Initialize a :class:`Bond`.

        Parameters
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
            `atom1` to `atom2` the bond is
            periodic across the x axis in the positive direction, is
            not periodic across the y axis and is periodic across the z
            axis in the negative direction.

        """

        self._atom1 = atom1
        self._atom2 = atom2
        self._order = order
        self._periodicity = periodicity

    def get_atom1(self):
        """
        Get the first atom of the bond.

        Returns
        -------
        :class:`.Atom`
            The first atom of the bond.

        """

        return self._atom1

    def get_atom2(self):
        """
        Get the second atom of the bond.

        Returns
        -------
        :class:`.Atom`
            The second atom of the bond.

        """

        return self._atom2

    def get_order(self):
        """
        Get the bond order of the bond.

        Returns
        -------
        :class:`int`
            The bond order.

        """

        return self._order

    def get_periodicity(self):
        """
        Get the periodicity of the bond.

        Returns
        -------
        :class:`tuple` of :class:`int`
            The directions across which the bond is periodic. For
            example, ``(1, 0, -1)`` means that when going from
            `atom1` to `atom2` the bond is
            periodic across the x axis in the positive direction, is
            not periodic across the y axis and is periodic across the z
            axis in the negative direction.

        """

        return self._periodicity

    def clone(self):
        """
        Return a clone.

        Returns
        -------
        :class:`.Bond`
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

    def _with_atoms(self, atom_map):
        """
        Modify the bond.

        """

        self._atom1 = atom_map.get(self._atom1.get_id(), self._atom1)
        self._atom2 = atom_map.get(self._atom2.get_id(), self._atom2)
        return self

    def with_atoms(self, atom_map):
        """
        Return a clone holding different atoms.

        Parameters
        ----------
        atom_map : :class:`dict`
            Maps the id of an atom in the bond to the new atom the
            clone should hold. If the id of an atom in the bond is not
            found in `atom_map`, the atom will not be replaced in the
            clone.

        Returns
        -------
        :class:`.Bond`
            The clone.

        Examples
        --------
        .. code-block:: python

            import stk
            bond = stk.Bond(stk.C(0), stk.C(12), 2)
            # Replace C(0) with H(13) but keep C(12).
            clone = bond.with_atoms({0: stk.H(13)})

        """

        return self.clone()._with_atoms(atom_map)

    def is_periodic(self):
        """
        Return ``True`` if the bond is periodic.

        Returns
        -------
        :class:`bool`
            ``True`` if the bond is periodic.

        """

        return any(direction != 0 for direction in self._periodicity)

    def __repr__(self):
        if isinstance(self._order, float) and self._order.is_integer():
            self._order = int(self._order)

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

    def __str__(self):
        return repr(self)
