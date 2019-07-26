class Bond:
    """
    Represents an atomic bond.

    Attributes
    ----------
    atom1 : :class:`Atom`
        The first atom in the bond.

    atom2 : :class:`Atom`
        The second atom in the bond.

    order : :class:`int`
        The bond order.

    Methods
    -------
    :meth:`__init__`
    :meth:`clone`

    """

    def __init__(self, atom1, atom2, order):
        """
        Intialize a :class:`Bond`.

        Attributes
        ----------
        atom1 : :class:`Atom`
            The first atom in the bond.

        atom2 : :class:`Atom`
            The second atom in the bond.

        order : :class:`int`
            The bond order.

        """

        self.atom1 = atom1
        self.atom2 = atom2
        self.order = order

    def clone(self, atom_map=None):
        """
        Clone the bond.

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
        return f'{cls_name}({self.atom1}, {self.atom2}, {self.order})'

    def __str__(self):
        return repr(self)


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
