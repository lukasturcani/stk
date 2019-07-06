class Bond:
    """
    Represents an atomic bond.

    Attributes
    ----------
    order : :class:`int`
        The bond order.

    atom1 : :class:`Atom`
        The first atom in the bond.

    atom2 : :class:`Atom`
        The second atom in the bond.

    """

    def __init__(self, order, atom1, atom2):
        self.order = order
        self.atom1 = atom1
        self.atom2 = atom2


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

    _elements : :class:`dict`
        Maps an atomic number to the class for that element.

    """

    _elements = {}

    def __init_subclass__(cls, **kwargs):
        cls._elements[cls.atomic_number] = cls

    def __init__(self, atomic_number):
        self.__class__ = self._elements[atomic_number]


class H(Atom):
    atomic_number = 1
    mass = 1.008


class He(Atom):
    atomic_number = 2
    mass = 4.003


class Li(Atom):
    atomic_number = 3
    mass = 6.941


class Be(Atom):
    atomic_number = 4
    mass = 9.012


class B(Atom):
    atomic_number = 5
    mass = 10.811


class C(Atom):
    atomic_number = 6
    mass = 12.011


class N(Atom):
    atomic_number = 7
    mass = 14.007


class O(Atom):
    atomic_number = 8
    mass = 15.999


class F(Atom):
    atomic_number = 9
    mass = 18.998


class Ne(Atom):
    atomic_number = 10
    mass = 20.180


class Na(Atom):
    atomic_nubmer = 11
    mass = 22.990


class Mg(Atom):
    atomic_number = 12
    mass = 24.305


class Al(Atom):
    atomic_number = 13
    mass = 26.982


class Si(Atom):
    atomic_number = 14
    mass = 28.086


class P(Atom):
    atomic_number = 15
    mass = 30.974


class S(Atom):
    atomic_number = 16
    mass = 32.066


class Cl(Atom):
    atomic_nubmer = 17
    mass = 35.453


class Ar(Atom):
    atomic_number = 18
    mass = 39.948


class K(Atom):
    atomic_nubmer = 19
    mass = 39.098


class Ca(Atom):
    atomic_number = 20
    mass = 40.078


class Sc(Atom):
    atomic_number = 21
    mass = 44.956


class Ti(Atom):
    atomic_number = 22
    mass = 47.867


class V(Atom):
    atomic_number = 23
    mass = 50.942


class Cr(Atom):
    atomic_number = 24
    mass = 51.996


class Mn(Atom):
    atomic_number = 25
    mass = 54.938


class Fe(Atom):
    atomic_number = 26
    mass = 55.845


class Co(Atom):
    atomic_number = 27
    mass = 58.933


class Ni(Atom):
    atomic_number = 28
    mass = 58.693


class Cu(Atom):
    atomic_number = 29
    mass = 63.546


class Zn(Atom):
    atomic_number = 30
    mass = 65.38


class Ga(Atom):
    atomic_number = 31
    mass = 69.723


class Ge(Atom):
    atomic_number = 32
    mass = 72.631


class As(Atom):
    atomic_number = 33
    mass = 74.922


class Se(Atom):
    atomic_number = 34
    mass = 78.971


class Br(Atom):
    atomic_number = 35
    mass = 79.904


class Kr(Atom):
    atomic_number = 36
    mass = 83.798


class Rb(Atom):
    atomic_number = 37
    mass = 85.468


class Sr(Atom):
    atomic_number = 38
    mass = 87.62


class Y(Atom):
    atomic_number = 39
    mass = 88.906


class Zr(Atom):
    atomic_number = 40
    mass = 91.224


class Nb(Atom):
    atomic_number = 41
    mass = 92.906


class Mo(Atom):
    atomic_number = 42
    mass = 95.95


class Tc(Atom):
    atomic_number = 43
    mass = 98.907


class Ru(Atom):
    atomic_number = 44
    mass = 101.07


class Rh(Atom):
    atomic_number = 45
    mass = 102.906


class Pd(Atom):
    atomic_number = 46
    mass = 106.42


class Ag(Atom):
    atomic_number = 47
    mass = 107.868


class Cd(Atom):
    atomic_number = 48
    mass = 112.414


class In(Atom):
    atomic_number = 49
    mass = 114.818


class Sn(Atom):
    atomic_number = 50
    mass = 118.711


class Sb(Atom):
    atomic_number = 51
    mass = 121.760


class Te(Atom):
    atomic_number = 52
    mass = 127.6


class I(Atom):
    atomic_number = 53
    mass = 126.904


class Xe(Atom):
    atomic_number = 54
    mass = 131.293


class Cs(Atom):
    atomic_number = 55
    mass = 132.905


class Ba(Atom):
    atomic_number = 56
    mass = 137.328


class La(Atom):
    atomic_number = 57
    mass = 138.905


class Ce(Atom):
    atomic_number = 58
    mass = 140.116


class Pr(Atom):
    atomic_number = 59
    mass = 140.908


class Nd(Atom):
    atomic_number = 60
    mass = 144.243


class Pm(Atom):
    atomic_number = 61
    mass = 144.913


class Sm(Atom):
    atomic_number = 62
    mass = 150.36


class Eu(Atom):
    atomic_number = 63
    mass = 151.964


class Gd(Atom):
    atomic_number = 64
    mass = 157.25


class Tb(Atom):
    atomic_number = 65
    mass = 158.925


class Dy(Atom):
    atomic_number = 66
    mass = 162.500


class Ho(Atom):
    atomic_number = 67
    mass = 164.930


class Er(Atom):
    atomic_number = 68
    mass = 167.259


class Tm(Atom):
    atomic_number = 69
    mass = 168.934


class Yb(Atom):
    atomic_nubmer = 70
    mass = 173.055


class Lu(Atom):
    atomic_number = 71
    mass = 174.967


class Hf(Atom):
    atomic_number = 72
    mass = 178.49


class Ta(Atom):
    atomic_number = 73
    mass = 180.948


class W(Atom):
    atomic_number = 74
    mass = 183.84


class Re(Atom):
    atomic_number = 75
    mass = 186.207


class Os(Atom):
    atomic_number = 76
    mass = 190.23


class Ir(Atom):
    atomic_number = 77
    mass = 192.217


class Pt(Atom):
    atomic_number = 78
    mass = 195.085


class Au(Atom):
    atomic_number = 79
    mass = 196.967


class Hg(Atom):
    atomic_number = 80
    mass = 200.592


class Tl(Atom):
    atomic_number = 81
    mass = 204.383


class Pb(Atom):
    atomic_number = 82
    mass = 207.2


class Bi(Atom):
    atomic_number = 83
    mass = 208.980


class Po(Atom):
    atomic_number = 84
    mass = 208.982


class At(Atom):
    atomic_number = 85
    mass = 209.987


class Rn(Atom):
    atomic_number = 86
    mass = 222.018


class Fr(Atom):
    atomic_number = 87
    mass = 223.020


class Ra(Atom):
    atomic_number = 88
    mass = 226.025


class Ac(Atom):
    atomic_number = 89
    mass = 227.028


class Th(Atom):
    atomic_number = 90
    mass = 232.038


class Pa(Atom):
    atomic_number = 91
    mass = 231.036


class U(Atom):
    atomic_number = 92
    mass = 238.029


class Np(Atom):
    atomic_number = 93
    mass = 237.048


class Pu(Atom):
    atomic_number = 94
    mass = 244.064


class Am(Atom):
    atomic_number = 95
    mass = 243.061


class Cm(Atom):
    atomic_number = 96
    mass = 247.070


class Bk(Atom):
    atomic_number = 97
    mass = 247.070


class Cf(Atom):
    atomic_number = 98
    mass = 251.080


class Es(Atom):
    atomic_number = 99
    mass = 254


class Fm(Atom):
    atomic_number = 100
    mass = 257.095


class Md(Atom):
    atomic_number = 101
    mass = 258.1


class No(Atom):
    atomic_number = 102
    mass = 259.101


class Lr(Atom):
    atomic_number = 103
    mass = 262


class Rf(Atom):
    atomic_number = 104
    mass = 261


class Db(Atom):
    atomic_number = 105
    mass = 262


class Sg(Atom):
    atomic_number = 106
    mass = 266


class Bh(Atom):
    atomic_number = 107
    mass = 264


class Hs(Atom):
    atomic_number = 108
    mass = 269


class Mt(Atom):
    atomic_number = 109
    mass = 278


class Ds(Atom):
    atomic_number = 110
    mass = 281


class Rg(Atom):
    atomic_number = 111
    mass = 280


class Cn(Atom):
    atomic_number = 112
    mass = 285


class Nh(Atom):
    atomic_number = 113
    mass = 286


class Fl(Atom):
    atomic_number = 114
    mass = 289


class Mc(Atom):
    atomic_number = 115
    mass = 289


class Lv(Atom):
    atomic_number = 116
    mass = 293


class Ts(Atom):
    atomic_number = 117
    mass = 294


class Og(Atom):
    atomic_number = 118
    mass = 294
