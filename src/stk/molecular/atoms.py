"""
Defines a class for each type of element.

Read the docstring of :class:`.Atom`.

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
        # The init method for subclasses takes a slightly different
        # form.
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
            A clone with a new id.

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

    def get_mass(self):
        """
        Get the mass of the atom.

        Returns
        -------
        :class:`float`
            The mass.

        """

        return self._mass

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

        Private attributes are not passed to the clone.

        Returns
        -------
        :class:`.Atom`
            The clone.

        """

        clone = self.__class__.__new__(self.__class__)
        clone._id = self._id
        clone._charge = self._charge
        for attr, val in vars(self).items():
            if not attr.startswith('_'):
                setattr(clone, attr, val)
        return clone


class H(Atom):
    _atomic_number = 1
    _mass = 1.008


class He(Atom):
    _atomic_number = 2
    _mass = 4.003


class Li(Atom):
    _atomic_number = 3
    _mass = 6.941


class Be(Atom):
    _atomic_number = 4
    _mass = 9.012


class B(Atom):
    _atomic_number = 5
    _mass = 10.811


class C(Atom):
    _atomic_number = 6
    _mass = 12.011


class N(Atom):
    _atomic_number = 7
    _mass = 14.007


class O(Atom):
    _atomic_number = 8
    _mass = 15.999


class F(Atom):
    _atomic_number = 9
    _mass = 18.998


class Ne(Atom):
    _atomic_number = 10
    _mass = 20.180


class Na(Atom):
    _atomic_number = 11
    _mass = 22.990


class Mg(Atom):
    _atomic_number = 12
    _mass = 24.305


class Al(Atom):
    _atomic_number = 13
    _mass = 26.982


class Si(Atom):
    _atomic_number = 14
    _mass = 28.086


class P(Atom):
    _atomic_number = 15
    _mass = 30.974


class S(Atom):
    _atomic_number = 16
    _mass = 32.066


class Cl(Atom):
    _atomic_number = 17
    _mass = 35.453


class Ar(Atom):
    _atomic_number = 18
    _mass = 39.948


class K(Atom):
    _atomic_number = 19
    _mass = 39.098


class Ca(Atom):
    _atomic_number = 20
    _mass = 40.078


class Sc(Atom):
    _atomic_number = 21
    _mass = 44.956


class Ti(Atom):
    _atomic_number = 22
    _mass = 47.867


class V(Atom):
    _atomic_number = 23
    _mass = 50.942


class Cr(Atom):
    _atomic_number = 24
    _mass = 51.996


class Mn(Atom):
    _atomic_number = 25
    _mass = 54.938


class Fe(Atom):
    _atomic_number = 26
    _mass = 55.845


class Co(Atom):
    _atomic_number = 27
    _mass = 58.933


class Ni(Atom):
    _atomic_number = 28
    _mass = 58.693


class Cu(Atom):
    _atomic_number = 29
    _mass = 63.546


class Zn(Atom):
    _atomic_number = 30
    _mass = 65.38


class Ga(Atom):
    _atomic_number = 31
    _mass = 69.723


class Ge(Atom):
    _atomic_number = 32
    _mass = 72.631


class As(Atom):
    _atomic_number = 33
    _mass = 74.922


class Se(Atom):
    _atomic_number = 34
    _mass = 78.971


class Br(Atom):
    _atomic_number = 35
    _mass = 79.904


class Kr(Atom):
    _atomic_number = 36
    _mass = 83.798


class Rb(Atom):
    _atomic_number = 37
    _mass = 85.468


class Sr(Atom):
    _atomic_number = 38
    _mass = 87.62


class Y(Atom):
    _atomic_number = 39
    _mass = 88.906


class Zr(Atom):
    _atomic_number = 40
    _mass = 91.224


class Nb(Atom):
    _atomic_number = 41
    _mass = 92.906


class Mo(Atom):
    _atomic_number = 42
    _mass = 95.95


class Tc(Atom):
    _atomic_number = 43
    _mass = 98.907


class Ru(Atom):
    _atomic_number = 44
    _mass = 101.07


class Rh(Atom):
    _atomic_number = 45
    _mass = 102.906


class Pd(Atom):
    _atomic_number = 46
    _mass = 106.42


class Ag(Atom):
    _atomic_number = 47
    _mass = 107.868


class Cd(Atom):
    _atomic_number = 48
    _mass = 112.414


class In(Atom):
    _atomic_number = 49
    _mass = 114.818


class Sn(Atom):
    _atomic_number = 50
    _mass = 118.711


class Sb(Atom):
    _atomic_number = 51
    _mass = 121.760


class Te(Atom):
    _atomic_number = 52
    _mass = 127.6


class I(Atom):
    _atomic_number = 53
    _mass = 126.904


class Xe(Atom):
    _atomic_number = 54
    _mass = 131.293


class Cs(Atom):
    _atomic_number = 55
    _mass = 132.905


class Ba(Atom):
    _atomic_number = 56
    _mass = 137.328


class La(Atom):
    _atomic_number = 57
    _mass = 138.905


class Ce(Atom):
    _atomic_number = 58
    _mass = 140.116


class Pr(Atom):
    _atomic_number = 59
    _mass = 140.908


class Nd(Atom):
    _atomic_number = 60
    _mass = 144.243


class Pm(Atom):
    _atomic_number = 61
    _mass = 144.913


class Sm(Atom):
    _atomic_number = 62
    _mass = 150.36


class Eu(Atom):
    _atomic_number = 63
    _mass = 151.964


class Gd(Atom):
    _atomic_number = 64
    _mass = 157.25


class Tb(Atom):
    _atomic_number = 65
    _mass = 158.925


class Dy(Atom):
    _atomic_number = 66
    _mass = 162.500


class Ho(Atom):
    _atomic_number = 67
    _mass = 164.930


class Er(Atom):
    _atomic_number = 68
    _mass = 167.259


class Tm(Atom):
    _atomic_number = 69
    _mass = 168.934


class Yb(Atom):
    _atomic_number = 70
    _mass = 173.055


class Lu(Atom):
    _atomic_number = 71
    _mass = 174.967


class Hf(Atom):
    _atomic_number = 72
    _mass = 178.49


class Ta(Atom):
    _atomic_number = 73
    _mass = 180.948


class W(Atom):
    _atomic_number = 74
    _mass = 183.84


class Re(Atom):
    _atomic_number = 75
    _mass = 186.207


class Os(Atom):
    _atomic_number = 76
    _mass = 190.23


class Ir(Atom):
    _atomic_number = 77
    _mass = 192.217


class Pt(Atom):
    _atomic_number = 78
    _mass = 195.085


class Au(Atom):
    _atomic_number = 79
    _mass = 196.967


class Hg(Atom):
    _atomic_number = 80
    _mass = 200.592


class Tl(Atom):
    _atomic_number = 81
    _mass = 204.383


class Pb(Atom):
    _atomic_number = 82
    _mass = 207.2


class Bi(Atom):
    _atomic_number = 83
    _mass = 208.980


class Po(Atom):
    _atomic_number = 84
    _mass = 208.982


class At(Atom):
    _atomic_number = 85
    _mass = 209.987


class Rn(Atom):
    _atomic_number = 86
    _mass = 222.018


class Fr(Atom):
    _atomic_number = 87
    _mass = 223.020


class Ra(Atom):
    _atomic_number = 88
    _mass = 226.025


class Ac(Atom):
    _atomic_number = 89
    _mass = 227.028


class Th(Atom):
    _atomic_number = 90
    _mass = 232.038


class Pa(Atom):
    _atomic_number = 91
    _mass = 231.036


class U(Atom):
    _atomic_number = 92
    _mass = 238.029


class Np(Atom):
    _atomic_number = 93
    _mass = 237.048


class Pu(Atom):
    _atomic_number = 94
    _mass = 244.064


class Am(Atom):
    _atomic_number = 95
    _mass = 243.061


class Cm(Atom):
    _atomic_number = 96
    _mass = 247.070


class Bk(Atom):
    _atomic_number = 97
    _mass = 247.070


class Cf(Atom):
    _atomic_number = 98
    _mass = 251.080


class Es(Atom):
    _atomic_number = 99
    _mass = 254


class Fm(Atom):
    _atomic_number = 100
    _mass = 257.095


class Md(Atom):
    _atomic_number = 101
    _mass = 258.1


class No(Atom):
    _atomic_number = 102
    _mass = 259.101


class Lr(Atom):
    _atomic_number = 103
    _mass = 262


class Rf(Atom):
    _atomic_number = 104
    _mass = 261


class Db(Atom):
    _atomic_number = 105
    _mass = 262


class Sg(Atom):
    _atomic_number = 106
    _mass = 266


class Bh(Atom):
    _atomic_number = 107
    _mass = 264


class Hs(Atom):
    _atomic_number = 108
    _mass = 269


class Mt(Atom):
    _atomic_number = 109
    _mass = 278


class Ds(Atom):
    _atomic_number = 110
    _mass = 281


class Rg(Atom):
    _atomic_number = 111
    _mass = 280


class Cn(Atom):
    _atomic_number = 112
    _mass = 285


class Nh(Atom):
    _atomic_number = 113
    _mass = 286


class Fl(Atom):
    _atomic_number = 114
    _mass = 289


class Mc(Atom):
    _atomic_number = 115
    _mass = 289


class Lv(Atom):
    _atomic_number = 116
    _mass = 293


class Ts(Atom):
    _atomic_number = 117
    _mass = 294


class Og(Atom):
    _atomic_number = 118
    _mass = 294
