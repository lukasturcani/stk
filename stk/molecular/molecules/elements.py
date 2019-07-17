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

    """

    def __init__(self, atom1, atom2, order):
        self.atom1 = atom1
        self.atom2 = atom2
        self.order = order

    def __repr__(self):
        return f'Bond({self.atom1}, {self.atom2}, {self.order})'

    def __str__(self):
        return repr(self)


class Atom:
    """
    Base class for atoms.

    A subclass is made for each element. The name of each elemental
    class is the periodic table symbol.

    Attributes
    ----------
    atomic_number : :class:`int`
        A class attribute. Specifies the atomic number.

    mass : :class:`int`
        A class attribute. Specifies the standard atmoic weight.

    id : :class:`int`
        The id of the atom. Should be equal to its index in
        :attr:`.Molecule.atoms`.

    charge : :class:`int`
        The formal charge of the atom.

    _elements : :class:`dict`
        Maps an atomic number to the class for that element.

    """

    _elements = {}

    def __init_subclass__(cls, **kwargs):
        cls._elements[cls.atomic_number] = cls

    def __init__(self, id, atomic_number, charge=0, **kwargs):
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

        **kwargs : :class:`object`
            Additional attributes to be added to the atom.

        """

        self.id = id
        self.charge = charge
        self.__class__ = self._elements[atomic_number]
        for attr, val in kwargs.items():
            setattr(self, attr, val)

    def __repr__(self):
        charge = f', charge={self.charge}' if self.charge != 0 else ''
        mandatory = {'charge', 'id'}
        attrs = ', '.join(
            f'{attr}={val}' for attr, val in vars(self).items()
            if attr not in mandatory and not attr.startswith('_')
        )
        return (
            f'{self.__class__.__name__}'
            f'({self.id}{charge}{attrs})'
        )

    def __str__(self):
        return repr(self)


class _Atom(Atom):
    atomic_number = float('nan')

    def __init__(self, id, charge=0, **kwargs):
        super().__init__(id, self.atomic_number, charge, **kwargs)


class H(_Atom):
    atomic_number = 1
    mass = 1.008


class He(_Atom):
    atomic_number = 2
    mass = 4.003


class Li(_Atom):
    atomic_number = 3
    mass = 6.941


class Be(_Atom):
    atomic_number = 4
    mass = 9.012


class B(_Atom):
    atomic_number = 5
    mass = 10.811


class C(_Atom):
    atomic_number = 6
    mass = 12.011


class N(_Atom):
    atomic_number = 7
    mass = 14.007


class O(_Atom):
    atomic_number = 8
    mass = 15.999


class F(_Atom):
    atomic_number = 9
    mass = 18.998


class Ne(_Atom):
    atomic_number = 10
    mass = 20.180


class Na(_Atom):
    atomic_nubmer = 11
    mass = 22.990


class Mg(_Atom):
    atomic_number = 12
    mass = 24.305


class Al(_Atom):
    atomic_number = 13
    mass = 26.982


class Si(_Atom):
    atomic_number = 14
    mass = 28.086


class P(_Atom):
    atomic_number = 15
    mass = 30.974


class S(_Atom):
    atomic_number = 16
    mass = 32.066


class Cl(_Atom):
    atomic_nubmer = 17
    mass = 35.453


class Ar(_Atom):
    atomic_number = 18
    mass = 39.948


class K(_Atom):
    atomic_nubmer = 19
    mass = 39.098


class Ca(_Atom):
    atomic_number = 20
    mass = 40.078


class Sc(_Atom):
    atomic_number = 21
    mass = 44.956


class Ti(_Atom):
    atomic_number = 22
    mass = 47.867


class V(_Atom):
    atomic_number = 23
    mass = 50.942


class Cr(_Atom):
    atomic_number = 24
    mass = 51.996


class Mn(_Atom):
    atomic_number = 25
    mass = 54.938


class Fe(_Atom):
    atomic_number = 26
    mass = 55.845


class Co(_Atom):
    atomic_number = 27
    mass = 58.933


class Ni(_Atom):
    atomic_number = 28
    mass = 58.693


class Cu(_Atom):
    atomic_number = 29
    mass = 63.546


class Zn(_Atom):
    atomic_number = 30
    mass = 65.38


class Ga(_Atom):
    atomic_number = 31
    mass = 69.723


class Ge(_Atom):
    atomic_number = 32
    mass = 72.631


class As(_Atom):
    atomic_number = 33
    mass = 74.922


class Se(_Atom):
    atomic_number = 34
    mass = 78.971


class Br(_Atom):
    atomic_number = 35
    mass = 79.904


class Kr(_Atom):
    atomic_number = 36
    mass = 83.798


class Rb(_Atom):
    atomic_number = 37
    mass = 85.468


class Sr(_Atom):
    atomic_number = 38
    mass = 87.62


class Y(_Atom):
    atomic_number = 39
    mass = 88.906


class Zr(_Atom):
    atomic_number = 40
    mass = 91.224


class Nb(_Atom):
    atomic_number = 41
    mass = 92.906


class Mo(_Atom):
    atomic_number = 42
    mass = 95.95


class Tc(_Atom):
    atomic_number = 43
    mass = 98.907


class Ru(_Atom):
    atomic_number = 44
    mass = 101.07


class Rh(_Atom):
    atomic_number = 45
    mass = 102.906


class Pd(_Atom):
    atomic_number = 46
    mass = 106.42


class Ag(_Atom):
    atomic_number = 47
    mass = 107.868


class Cd(_Atom):
    atomic_number = 48
    mass = 112.414


class In(_Atom):
    atomic_number = 49
    mass = 114.818


class Sn(_Atom):
    atomic_number = 50
    mass = 118.711


class Sb(_Atom):
    atomic_number = 51
    mass = 121.760


class Te(_Atom):
    atomic_number = 52
    mass = 127.6


class I(_Atom):
    atomic_number = 53
    mass = 126.904


class Xe(_Atom):
    atomic_number = 54
    mass = 131.293


class Cs(_Atom):
    atomic_number = 55
    mass = 132.905


class Ba(_Atom):
    atomic_number = 56
    mass = 137.328


class La(_Atom):
    atomic_number = 57
    mass = 138.905


class Ce(_Atom):
    atomic_number = 58
    mass = 140.116


class Pr(_Atom):
    atomic_number = 59
    mass = 140.908


class Nd(_Atom):
    atomic_number = 60
    mass = 144.243


class Pm(_Atom):
    atomic_number = 61
    mass = 144.913


class Sm(_Atom):
    atomic_number = 62
    mass = 150.36


class Eu(_Atom):
    atomic_number = 63
    mass = 151.964


class Gd(_Atom):
    atomic_number = 64
    mass = 157.25


class Tb(_Atom):
    atomic_number = 65
    mass = 158.925


class Dy(_Atom):
    atomic_number = 66
    mass = 162.500


class Ho(_Atom):
    atomic_number = 67
    mass = 164.930


class Er(_Atom):
    atomic_number = 68
    mass = 167.259


class Tm(_Atom):
    atomic_number = 69
    mass = 168.934


class Yb(_Atom):
    atomic_nubmer = 70
    mass = 173.055


class Lu(_Atom):
    atomic_number = 71
    mass = 174.967


class Hf(_Atom):
    atomic_number = 72
    mass = 178.49


class Ta(_Atom):
    atomic_number = 73
    mass = 180.948


class W(_Atom):
    atomic_number = 74
    mass = 183.84


class Re(_Atom):
    atomic_number = 75
    mass = 186.207


class Os(_Atom):
    atomic_number = 76
    mass = 190.23


class Ir(_Atom):
    atomic_number = 77
    mass = 192.217


class Pt(_Atom):
    atomic_number = 78
    mass = 195.085


class Au(_Atom):
    atomic_number = 79
    mass = 196.967


class Hg(_Atom):
    atomic_number = 80
    mass = 200.592


class Tl(_Atom):
    atomic_number = 81
    mass = 204.383


class Pb(_Atom):
    atomic_number = 82
    mass = 207.2


class Bi(_Atom):
    atomic_number = 83
    mass = 208.980


class Po(_Atom):
    atomic_number = 84
    mass = 208.982


class At(_Atom):
    atomic_number = 85
    mass = 209.987


class Rn(_Atom):
    atomic_number = 86
    mass = 222.018


class Fr(_Atom):
    atomic_number = 87
    mass = 223.020


class Ra(_Atom):
    atomic_number = 88
    mass = 226.025


class Ac(_Atom):
    atomic_number = 89
    mass = 227.028


class Th(_Atom):
    atomic_number = 90
    mass = 232.038


class Pa(_Atom):
    atomic_number = 91
    mass = 231.036


class U(_Atom):
    atomic_number = 92
    mass = 238.029


class Np(_Atom):
    atomic_number = 93
    mass = 237.048


class Pu(_Atom):
    atomic_number = 94
    mass = 244.064


class Am(_Atom):
    atomic_number = 95
    mass = 243.061


class Cm(_Atom):
    atomic_number = 96
    mass = 247.070


class Bk(_Atom):
    atomic_number = 97
    mass = 247.070


class Cf(_Atom):
    atomic_number = 98
    mass = 251.080


class Es(_Atom):
    atomic_number = 99
    mass = 254


class Fm(_Atom):
    atomic_number = 100
    mass = 257.095


class Md(_Atom):
    atomic_number = 101
    mass = 258.1


class No(_Atom):
    atomic_number = 102
    mass = 259.101


class Lr(_Atom):
    atomic_number = 103
    mass = 262


class Rf(_Atom):
    atomic_number = 104
    mass = 261


class Db(_Atom):
    atomic_number = 105
    mass = 262


class Sg(_Atom):
    atomic_number = 106
    mass = 266


class Bh(_Atom):
    atomic_number = 107
    mass = 264


class Hs(_Atom):
    atomic_number = 108
    mass = 269


class Mt(_Atom):
    atomic_number = 109
    mass = 278


class Ds(_Atom):
    atomic_number = 110
    mass = 281


class Rg(_Atom):
    atomic_number = 111
    mass = 280


class Cn(_Atom):
    atomic_number = 112
    mass = 285


class Nh(_Atom):
    atomic_number = 113
    mass = 286


class Fl(_Atom):
    atomic_number = 114
    mass = 289


class Mc(_Atom):
    atomic_number = 115
    mass = 289


class Lv(_Atom):
    atomic_number = 116
    mass = 293


class Ts(_Atom):
    atomic_number = 117
    mass = 294


class Og(_Atom):
    atomic_number = 118
    mass = 294
