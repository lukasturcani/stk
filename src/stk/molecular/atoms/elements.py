"""
Elements
========

Defines an :class:`.Atom` class for each element.

"""

from .atom import Atom


class AtomImpl(Atom):
    """
    A partial implementation of the :class:`.Atom` interface.

    """

    def __init__(
        self,
        id: int,
        charge: int = 0,
    ) -> None:

        Atom.__init__(self, id, self._atomic_number, charge)


class H(AtomImpl):
    _atomic_number = 1


class He(AtomImpl):
    _atomic_number = 2


class Li(AtomImpl):
    _atomic_number = 3


class Be(AtomImpl):
    _atomic_number = 4


class B(AtomImpl):
    _atomic_number = 5


class C(AtomImpl):
    _atomic_number = 6


class N(AtomImpl):
    _atomic_number = 7


# "O" is a valid elemental symbol.
class O(AtomImpl):  # noqa
    _atomic_number = 8


class F(AtomImpl):
    _atomic_number = 9


class Ne(AtomImpl):
    _atomic_number = 10


class Na(AtomImpl):
    _atomic_number = 11


class Mg(AtomImpl):
    _atomic_number = 12


class Al(AtomImpl):
    _atomic_number = 13


class Si(AtomImpl):
    _atomic_number = 14


class P(AtomImpl):
    _atomic_number = 15


class S(AtomImpl):
    _atomic_number = 16


class Cl(AtomImpl):
    _atomic_number = 17


class Ar(AtomImpl):
    _atomic_number = 18


class K(AtomImpl):
    _atomic_number = 19


class Ca(AtomImpl):
    _atomic_number = 20


class Sc(AtomImpl):
    _atomic_number = 21


class Ti(AtomImpl):
    _atomic_number = 22


class V(AtomImpl):
    _atomic_number = 23


class Cr(AtomImpl):
    _atomic_number = 24


class Mn(AtomImpl):
    _atomic_number = 25


class Fe(AtomImpl):
    _atomic_number = 26


class Co(AtomImpl):
    _atomic_number = 27


class Ni(AtomImpl):
    _atomic_number = 28


class Cu(AtomImpl):
    _atomic_number = 29


class Zn(AtomImpl):
    _atomic_number = 30


class Ga(AtomImpl):
    _atomic_number = 31


class Ge(AtomImpl):
    _atomic_number = 32


class As(AtomImpl):
    _atomic_number = 33


class Se(AtomImpl):
    _atomic_number = 34


class Br(AtomImpl):
    _atomic_number = 35


class Kr(AtomImpl):
    _atomic_number = 36


class Rb(AtomImpl):
    _atomic_number = 37


class Sr(AtomImpl):
    _atomic_number = 38


class Y(AtomImpl):
    _atomic_number = 39


class Zr(AtomImpl):
    _atomic_number = 40


class Nb(AtomImpl):
    _atomic_number = 41


class Mo(AtomImpl):
    _atomic_number = 42


class Tc(AtomImpl):
    _atomic_number = 43


class Ru(AtomImpl):
    _atomic_number = 44


class Rh(AtomImpl):
    _atomic_number = 45


class Pd(AtomImpl):
    _atomic_number = 46


class Ag(AtomImpl):
    _atomic_number = 47


class Cd(AtomImpl):
    _atomic_number = 48


class In(AtomImpl):
    _atomic_number = 49


class Sn(AtomImpl):
    _atomic_number = 50


class Sb(AtomImpl):
    _atomic_number = 51


class Te(AtomImpl):
    _atomic_number = 52


# "I" is a valid elemental symbol.
class I(AtomImpl):  # noqa
    _atomic_number = 53


class Xe(AtomImpl):
    _atomic_number = 54


class Cs(AtomImpl):
    _atomic_number = 55


class Ba(AtomImpl):
    _atomic_number = 56


class La(AtomImpl):
    _atomic_number = 57


class Ce(AtomImpl):
    _atomic_number = 58


class Pr(AtomImpl):
    _atomic_number = 59


class Nd(AtomImpl):
    _atomic_number = 60


class Pm(AtomImpl):
    _atomic_number = 61


class Sm(AtomImpl):
    _atomic_number = 62


class Eu(AtomImpl):
    _atomic_number = 63


class Gd(AtomImpl):
    _atomic_number = 64


class Tb(AtomImpl):
    _atomic_number = 65


class Dy(AtomImpl):
    _atomic_number = 66


class Ho(AtomImpl):
    _atomic_number = 67


class Er(AtomImpl):
    _atomic_number = 68


class Tm(AtomImpl):
    _atomic_number = 69


class Yb(AtomImpl):
    _atomic_number = 70


class Lu(AtomImpl):
    _atomic_number = 71


class Hf(AtomImpl):
    _atomic_number = 72


class Ta(AtomImpl):
    _atomic_number = 73


class W(AtomImpl):
    _atomic_number = 74


class Re(AtomImpl):
    _atomic_number = 75


class Os(AtomImpl):
    _atomic_number = 76


class Ir(AtomImpl):
    _atomic_number = 77


class Pt(AtomImpl):
    _atomic_number = 78


class Au(AtomImpl):
    _atomic_number = 79


class Hg(AtomImpl):
    _atomic_number = 80


class Tl(AtomImpl):
    _atomic_number = 81


class Pb(AtomImpl):
    _atomic_number = 82


class Bi(AtomImpl):
    _atomic_number = 83


class Po(AtomImpl):
    _atomic_number = 84


class At(AtomImpl):
    _atomic_number = 85


class Rn(AtomImpl):
    _atomic_number = 86


class Fr(AtomImpl):
    _atomic_number = 87


class Ra(AtomImpl):
    _atomic_number = 88


class Ac(AtomImpl):
    _atomic_number = 89


class Th(AtomImpl):
    _atomic_number = 90


class Pa(AtomImpl):
    _atomic_number = 91


class U(AtomImpl):
    _atomic_number = 92


class Np(AtomImpl):
    _atomic_number = 93


class Pu(AtomImpl):
    _atomic_number = 94


class Am(AtomImpl):
    _atomic_number = 95


class Cm(AtomImpl):
    _atomic_number = 96


class Bk(AtomImpl):
    _atomic_number = 97


class Cf(AtomImpl):
    _atomic_number = 98


class Es(AtomImpl):
    _atomic_number = 99


class Fm(AtomImpl):
    _atomic_number = 100


class Md(AtomImpl):
    _atomic_number = 101


class No(AtomImpl):
    _atomic_number = 102


class Lr(AtomImpl):
    _atomic_number = 103


class Rf(AtomImpl):
    _atomic_number = 104


class Db(AtomImpl):
    _atomic_number = 105


class Sg(AtomImpl):
    _atomic_number = 106


class Bh(AtomImpl):
    _atomic_number = 107


class Hs(AtomImpl):
    _atomic_number = 108


class Mt(AtomImpl):
    _atomic_number = 109


class Ds(AtomImpl):
    _atomic_number = 110


class Rg(AtomImpl):
    _atomic_number = 111


class Cn(AtomImpl):
    _atomic_number = 112


class Nh(AtomImpl):
    _atomic_number = 113


class Fl(AtomImpl):
    _atomic_number = 114


class Mc(AtomImpl):
    _atomic_number = 115


class Lv(AtomImpl):
    _atomic_number = 116


class Ts(AtomImpl):
    _atomic_number = 117


class Og(AtomImpl):
    _atomic_number = 118
