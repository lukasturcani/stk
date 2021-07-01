"""
Elements
========

Defines an :class:`.Atom` class for each element.

"""

from __future__ import annotations
from typing import ClassVar

from .atom import Atom


class AtomImpl(Atom):
    """
    An implementation of the :class:`.Atom` interface.

    """

    _atomic_number: ClassVar[int]

    def __init_subclass__(cls: type[AtomImpl], **kwargs) -> None:
        cls._elements[cls._atomic_number] = cls

    def __init__(
        self,
        id: int,
        charge: int = 0,
    ) -> None:

        self._id = id
        self._charge = charge

    def get_id(self) -> int:
        return self._id

    def _with_id(self, id: int) -> Atom:
        """
        Modify the atom.

        """

        self._id = id
        return self

    def with_id(self, id: int) -> Atom:
        return self.clone()._with_id(id)

    def get_atomic_number(self) -> int:
        return self._atomic_number

    def get_charge(self) -> int:
        return self._charge

    def clone(self) -> AtomImpl:
        return type(self)(self._id, self._charge)

    def __repr__(self) -> str:
        charge = (
            f', charge={self._charge}' if self._charge != 0 else ''
        )
        return f'{self.__class__.__name__}({self._id}{charge})'

    def __str__(self) -> str:
        return repr(self)


class H(AtomImpl):
    _atomic_number = 1

    def clone(self) -> H:
        return type(self)(self._id, self._charge)


class He(AtomImpl):
    _atomic_number = 2

    def clone(self) -> He:
        return type(self)(self._id, self._charge)


class Li(AtomImpl):
    _atomic_number = 3

    def clone(self) -> Li:
        return type(self)(self._id, self._charge)


class Be(AtomImpl):
    _atomic_number = 4

    def clone(self) -> Be:
        return type(self)(self._id, self._charge)


class B(AtomImpl):
    _atomic_number = 5

    def clone(self) -> B:
        return type(self)(self._id, self._charge)


class C(AtomImpl):
    _atomic_number = 6

    def clone(self) -> C:
        return type(self)(self._id, self._charge)


class N(AtomImpl):
    _atomic_number = 7

    def clone(self) -> N:
        return type(self)(self._id, self._charge)


# "O" is a valid elemental symbol.
class O(AtomImpl):  # noqa
    _atomic_number = 8

    def clone(self) -> O:  # noqa
        return type(self)(self._id, self._charge)


class F(AtomImpl):
    _atomic_number = 9

    def clone(self) -> F:
        return type(self)(self._id, self._charge)


class Ne(AtomImpl):
    _atomic_number = 10

    def clone(self) -> Ne:
        return type(self)(self._id, self._charge)


class Na(AtomImpl):
    _atomic_number = 11

    def clone(self) -> Na:
        return type(self)(self._id, self._charge)


class Mg(AtomImpl):
    _atomic_number = 12

    def clone(self) -> Mg:
        return type(self)(self._id, self._charge)


class Al(AtomImpl):
    _atomic_number = 13

    def clone(self) -> Al:
        return type(self)(self._id, self._charge)


class Si(AtomImpl):
    _atomic_number = 14

    def clone(self) -> Si:
        return type(self)(self._id, self._charge)


class P(AtomImpl):
    _atomic_number = 15

    def clone(self) -> P:
        return type(self)(self._id, self._charge)


class S(AtomImpl):
    _atomic_number = 16

    def clone(self) -> S:
        return type(self)(self._id, self._charge)


class Cl(AtomImpl):
    _atomic_number = 17

    def clone(self) -> Cl:
        return type(self)(self._id, self._charge)


class Ar(AtomImpl):
    _atomic_number = 18

    def clone(self) -> Ar:
        return type(self)(self._id, self._charge)


class K(AtomImpl):
    _atomic_number = 19

    def clone(self) -> K:
        return type(self)(self._id, self._charge)


class Ca(AtomImpl):
    _atomic_number = 20

    def clone(self) -> Ca:
        return type(self)(self._id, self._charge)


class Sc(AtomImpl):
    _atomic_number = 21

    def clone(self) -> Sc:
        return type(self)(self._id, self._charge)


class Ti(AtomImpl):
    _atomic_number = 22

    def clone(self) -> Ti:
        return type(self)(self._id, self._charge)


class V(AtomImpl):
    _atomic_number = 23

    def clone(self) -> V:
        return type(self)(self._id, self._charge)


class Cr(AtomImpl):
    _atomic_number = 24

    def clone(self) -> Cr:
        return type(self)(self._id, self._charge)


class Mn(AtomImpl):
    _atomic_number = 25

    def clone(self) -> Mn:
        return type(self)(self._id, self._charge)


class Fe(AtomImpl):
    _atomic_number = 26

    def clone(self) -> Fe:
        return type(self)(self._id, self._charge)


class Co(AtomImpl):
    _atomic_number = 27

    def clone(self) -> Co:
        return type(self)(self._id, self._charge)


class Ni(AtomImpl):
    _atomic_number = 28

    def clone(self) -> Ni:
        return type(self)(self._id, self._charge)


class Cu(AtomImpl):
    _atomic_number = 29

    def clone(self) -> Cu:
        return type(self)(self._id, self._charge)


class Zn(AtomImpl):
    _atomic_number = 30

    def clone(self) -> Zn:
        return type(self)(self._id, self._charge)


class Ga(AtomImpl):
    _atomic_number = 31

    def clone(self) -> Ga:
        return type(self)(self._id, self._charge)


class Ge(AtomImpl):
    _atomic_number = 32

    def clone(self) -> Ge:
        return type(self)(self._id, self._charge)


class As(AtomImpl):
    _atomic_number = 33

    def clone(self) -> As:
        return type(self)(self._id, self._charge)


class Se(AtomImpl):
    _atomic_number = 34

    def clone(self) -> Se:
        return type(self)(self._id, self._charge)


class Br(AtomImpl):
    _atomic_number = 35

    def clone(self) -> Br:
        return type(self)(self._id, self._charge)


class Kr(AtomImpl):
    _atomic_number = 36

    def clone(self) -> Kr:
        return type(self)(self._id, self._charge)


class Rb(AtomImpl):
    _atomic_number = 37

    def clone(self) -> Rb:
        return type(self)(self._id, self._charge)


class Sr(AtomImpl):
    _atomic_number = 38

    def clone(self) -> Sr:
        return type(self)(self._id, self._charge)


class Y(AtomImpl):
    _atomic_number = 39

    def clone(self) -> Y:
        return type(self)(self._id, self._charge)


class Zr(AtomImpl):
    _atomic_number = 40

    def clone(self) -> Zr:
        return type(self)(self._id, self._charge)


class Nb(AtomImpl):
    _atomic_number = 41

    def clone(self) -> Nb:
        return type(self)(self._id, self._charge)


class Mo(AtomImpl):
    _atomic_number = 42

    def clone(self) -> Mo:
        return type(self)(self._id, self._charge)


class Tc(AtomImpl):
    _atomic_number = 43

    def clone(self) -> Tc:
        return type(self)(self._id, self._charge)


class Ru(AtomImpl):
    _atomic_number = 44

    def clone(self) -> Ru:
        return type(self)(self._id, self._charge)


class Rh(AtomImpl):
    _atomic_number = 45

    def clone(self) -> Rh:
        return type(self)(self._id, self._charge)


class Pd(AtomImpl):
    _atomic_number = 46

    def clone(self) -> Pd:
        return type(self)(self._id, self._charge)


class Ag(AtomImpl):
    _atomic_number = 47

    def clone(self) -> Ag:
        return type(self)(self._id, self._charge)


class Cd(AtomImpl):
    _atomic_number = 48

    def clone(self) -> Cd:
        return type(self)(self._id, self._charge)


class In(AtomImpl):
    _atomic_number = 49

    def clone(self) -> In:
        return type(self)(self._id, self._charge)


class Sn(AtomImpl):
    _atomic_number = 50

    def clone(self) -> Sn:
        return type(self)(self._id, self._charge)


class Sb(AtomImpl):
    _atomic_number = 51

    def clone(self) -> Sb:
        return type(self)(self._id, self._charge)


class Te(AtomImpl):
    _atomic_number = 52

    def clone(self) -> Te:
        return type(self)(self._id, self._charge)


# "I" is a valid elemental symbol.
class I(AtomImpl):  # noqa
    _atomic_number = 53

    def clone(self) -> I:  # noqa
        return type(self)(self._id, self._charge)


class Xe(AtomImpl):
    _atomic_number = 54

    def clone(self) -> Xe:
        return type(self)(self._id, self._charge)


class Cs(AtomImpl):
    _atomic_number = 55

    def clone(self) -> Cs:
        return type(self)(self._id, self._charge)


class Ba(AtomImpl):
    _atomic_number = 56

    def clone(self) -> Ba:
        return type(self)(self._id, self._charge)


class La(AtomImpl):
    _atomic_number = 57

    def clone(self) -> La:
        return type(self)(self._id, self._charge)


class Ce(AtomImpl):
    _atomic_number = 58

    def clone(self) -> Ce:
        return type(self)(self._id, self._charge)


class Pr(AtomImpl):
    _atomic_number = 59

    def clone(self) -> Pr:
        return type(self)(self._id, self._charge)


class Nd(AtomImpl):
    _atomic_number = 60

    def clone(self) -> Nd:
        return type(self)(self._id, self._charge)


class Pm(AtomImpl):
    _atomic_number = 61

    def clone(self) -> Pm:
        return type(self)(self._id, self._charge)


class Sm(AtomImpl):
    _atomic_number = 62

    def clone(self) -> Sm:
        return type(self)(self._id, self._charge)


class Eu(AtomImpl):
    _atomic_number = 63

    def clone(self) -> Eu:
        return type(self)(self._id, self._charge)


class Gd(AtomImpl):
    _atomic_number = 64

    def clone(self) -> Gd:
        return type(self)(self._id, self._charge)


class Tb(AtomImpl):
    _atomic_number = 65

    def clone(self) -> Tb:
        return type(self)(self._id, self._charge)


class Dy(AtomImpl):
    _atomic_number = 66

    def clone(self) -> Dy:
        return type(self)(self._id, self._charge)


class Ho(AtomImpl):
    _atomic_number = 67

    def clone(self) -> Ho:
        return type(self)(self._id, self._charge)


class Er(AtomImpl):
    _atomic_number = 68

    def clone(self) -> Er:
        return type(self)(self._id, self._charge)


class Tm(AtomImpl):
    _atomic_number = 69

    def clone(self) -> Tm:
        return type(self)(self._id, self._charge)


class Yb(AtomImpl):
    _atomic_number = 70

    def clone(self) -> Yb:
        return type(self)(self._id, self._charge)


class Lu(AtomImpl):
    _atomic_number = 71

    def clone(self) -> Lu:
        return type(self)(self._id, self._charge)


class Hf(AtomImpl):
    _atomic_number = 72

    def clone(self) -> Hf:
        return type(self)(self._id, self._charge)


class Ta(AtomImpl):
    _atomic_number = 73

    def clone(self) -> Ta:
        return type(self)(self._id, self._charge)


class W(AtomImpl):
    _atomic_number = 74

    def clone(self) -> W:
        return type(self)(self._id, self._charge)


class Re(AtomImpl):
    _atomic_number = 75

    def clone(self) -> Re:
        return type(self)(self._id, self._charge)


class Os(AtomImpl):
    _atomic_number = 76

    def clone(self) -> Os:
        return type(self)(self._id, self._charge)


class Ir(AtomImpl):
    _atomic_number = 77

    def clone(self) -> Ir:
        return type(self)(self._id, self._charge)


class Pt(AtomImpl):
    _atomic_number = 78

    def clone(self) -> Pt:
        return type(self)(self._id, self._charge)


class Au(AtomImpl):
    _atomic_number = 79

    def clone(self) -> Au:
        return type(self)(self._id, self._charge)


class Hg(AtomImpl):
    _atomic_number = 80

    def clone(self) -> Hg:
        return type(self)(self._id, self._charge)


class Tl(AtomImpl):
    _atomic_number = 81

    def clone(self) -> Tl:
        return type(self)(self._id, self._charge)


class Pb(AtomImpl):
    _atomic_number = 82

    def clone(self) -> Pb:
        return type(self)(self._id, self._charge)


class Bi(AtomImpl):
    _atomic_number = 83

    def clone(self) -> Bi:
        return type(self)(self._id, self._charge)


class Po(AtomImpl):
    _atomic_number = 84

    def clone(self) -> Po:
        return type(self)(self._id, self._charge)


class At(AtomImpl):
    _atomic_number = 85

    def clone(self) -> At:
        return type(self)(self._id, self._charge)


class Rn(AtomImpl):
    _atomic_number = 86

    def clone(self) -> Rn:
        return type(self)(self._id, self._charge)


class Fr(AtomImpl):
    _atomic_number = 87

    def clone(self) -> Fr:
        return type(self)(self._id, self._charge)


class Ra(AtomImpl):
    _atomic_number = 88

    def clone(self) -> Ra:
        return type(self)(self._id, self._charge)


class Ac(AtomImpl):
    _atomic_number = 89

    def clone(self) -> Ac:
        return type(self)(self._id, self._charge)


class Th(AtomImpl):
    _atomic_number = 90

    def clone(self) -> Th:
        return type(self)(self._id, self._charge)


class Pa(AtomImpl):
    _atomic_number = 91

    def clone(self) -> Pa:
        return type(self)(self._id, self._charge)


class U(AtomImpl):
    _atomic_number = 92

    def clone(self) -> U:
        return type(self)(self._id, self._charge)


class Np(AtomImpl):
    _atomic_number = 93

    def clone(self) -> Np:
        return type(self)(self._id, self._charge)


class Pu(AtomImpl):
    _atomic_number = 94

    def clone(self) -> Pu:
        return type(self)(self._id, self._charge)


class Am(AtomImpl):
    _atomic_number = 95

    def clone(self) -> Am:
        return type(self)(self._id, self._charge)


class Cm(AtomImpl):
    _atomic_number = 96

    def clone(self) -> Cm:
        return type(self)(self._id, self._charge)


class Bk(AtomImpl):
    _atomic_number = 97

    def clone(self) -> Bk:
        return type(self)(self._id, self._charge)


class Cf(AtomImpl):
    _atomic_number = 98

    def clone(self) -> Cf:
        return type(self)(self._id, self._charge)


class Es(AtomImpl):
    _atomic_number = 99

    def clone(self) -> Es:
        return type(self)(self._id, self._charge)


class Fm(AtomImpl):
    _atomic_number = 100

    def clone(self) -> Fm:
        return type(self)(self._id, self._charge)


class Md(AtomImpl):
    _atomic_number = 101

    def clone(self) -> Md:
        return type(self)(self._id, self._charge)


class No(AtomImpl):
    _atomic_number = 102

    def clone(self) -> No:
        return type(self)(self._id, self._charge)


class Lr(AtomImpl):
    _atomic_number = 103

    def clone(self) -> Lr:
        return type(self)(self._id, self._charge)


class Rf(AtomImpl):
    _atomic_number = 104

    def clone(self) -> Rf:
        return type(self)(self._id, self._charge)


class Db(AtomImpl):
    _atomic_number = 105

    def clone(self) -> Db:
        return type(self)(self._id, self._charge)


class Sg(AtomImpl):
    _atomic_number = 106

    def clone(self) -> Sg:
        return type(self)(self._id, self._charge)


class Bh(AtomImpl):
    _atomic_number = 107

    def clone(self) -> Bh:
        return type(self)(self._id, self._charge)


class Hs(AtomImpl):
    _atomic_number = 108

    def clone(self) -> Hs:
        return type(self)(self._id, self._charge)


class Mt(AtomImpl):
    _atomic_number = 109

    def clone(self) -> Mt:
        return type(self)(self._id, self._charge)


class Ds(AtomImpl):
    _atomic_number = 110

    def clone(self) -> Ds:
        return type(self)(self._id, self._charge)


class Rg(AtomImpl):
    _atomic_number = 111

    def clone(self) -> Rg:
        return type(self)(self._id, self._charge)


class Cn(AtomImpl):
    _atomic_number = 112

    def clone(self) -> Cn:
        return type(self)(self._id, self._charge)


class Nh(AtomImpl):
    _atomic_number = 113

    def clone(self) -> Nh:
        return type(self)(self._id, self._charge)


class Fl(AtomImpl):
    _atomic_number = 114

    def clone(self) -> Fl:
        return type(self)(self._id, self._charge)


class Mc(AtomImpl):
    _atomic_number = 115

    def clone(self) -> Mc:
        return type(self)(self._id, self._charge)


class Lv(AtomImpl):
    _atomic_number = 116

    def clone(self) -> Lv:
        return type(self)(self._id, self._charge)


class Ts(AtomImpl):
    _atomic_number = 117

    def clone(self) -> Ts:
        return type(self)(self._id, self._charge)


class Og(AtomImpl):
    _atomic_number = 118

    def clone(self) -> Og:
        return type(self)(self._id, self._charge)
