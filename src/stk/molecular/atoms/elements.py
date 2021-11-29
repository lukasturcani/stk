"""
Elements
========

Defines an :class:`.Atom` class for each element.

"""

from __future__ import annotations

from typing import ClassVar, TypeVar

from .atom import Atom

_T = TypeVar('_T', bound='AtomImpl')


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

    def _with_id(self: _T, id: int) -> _T:
        """
        Modify the atom.

        """

        self._id = id
        return self

    def with_id(self, id: int) -> AtomImpl:
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

    def with_id(self, id: int) -> H:
        return self.clone()._with_id(id)


class He(AtomImpl):
    _atomic_number = 2

    def clone(self) -> He:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> He:
        return self.clone()._with_id(id)


class Li(AtomImpl):
    _atomic_number = 3

    def clone(self) -> Li:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Li:
        return self.clone()._with_id(id)


class Be(AtomImpl):
    _atomic_number = 4

    def clone(self) -> Be:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Be:
        return self.clone()._with_id(id)


class B(AtomImpl):
    _atomic_number = 5

    def clone(self) -> B:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> B:
        return self.clone()._with_id(id)


class C(AtomImpl):
    _atomic_number = 6

    def clone(self) -> C:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> C:
        return self.clone()._with_id(id)


class N(AtomImpl):
    _atomic_number = 7

    def clone(self) -> N:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> N:
        return self.clone()._with_id(id)


# "O" is a valid elemental symbol.
class O(AtomImpl):  # noqa
    _atomic_number = 8

    def clone(self) -> O:  # noqa
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> O:  # noqa
        return self.clone()._with_id(id)


class F(AtomImpl):
    _atomic_number = 9

    def clone(self) -> F:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> F:
        return self.clone()._with_id(id)


class Ne(AtomImpl):
    _atomic_number = 10

    def clone(self) -> Ne:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Ne:
        return self.clone()._with_id(id)


class Na(AtomImpl):
    _atomic_number = 11

    def clone(self) -> Na:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Na:
        return self.clone()._with_id(id)


class Mg(AtomImpl):
    _atomic_number = 12

    def clone(self) -> Mg:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Mg:
        return self.clone()._with_id(id)


class Al(AtomImpl):
    _atomic_number = 13

    def clone(self) -> Al:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Al:
        return self.clone()._with_id(id)


class Si(AtomImpl):
    _atomic_number = 14

    def clone(self) -> Si:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Si:
        return self.clone()._with_id(id)


class P(AtomImpl):
    _atomic_number = 15

    def clone(self) -> P:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> P:
        return self.clone()._with_id(id)


class S(AtomImpl):
    _atomic_number = 16

    def clone(self) -> S:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> S:
        return self.clone()._with_id(id)


class Cl(AtomImpl):
    _atomic_number = 17

    def clone(self) -> Cl:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Cl:
        return self.clone()._with_id(id)


class Ar(AtomImpl):
    _atomic_number = 18

    def clone(self) -> Ar:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Ar:
        return self.clone()._with_id(id)


class K(AtomImpl):
    _atomic_number = 19

    def clone(self) -> K:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> K:
        return self.clone()._with_id(id)


class Ca(AtomImpl):
    _atomic_number = 20

    def clone(self) -> Ca:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Ca:
        return self.clone()._with_id(id)


class Sc(AtomImpl):
    _atomic_number = 21

    def clone(self) -> Sc:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Sc:
        return self.clone()._with_id(id)


class Ti(AtomImpl):
    _atomic_number = 22

    def clone(self) -> Ti:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Ti:
        return self.clone()._with_id(id)


class V(AtomImpl):
    _atomic_number = 23

    def clone(self) -> V:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> V:
        return self.clone()._with_id(id)


class Cr(AtomImpl):
    _atomic_number = 24

    def clone(self) -> Cr:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Cr:
        return self.clone()._with_id(id)


class Mn(AtomImpl):
    _atomic_number = 25

    def clone(self) -> Mn:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Mn:
        return self.clone()._with_id(id)


class Fe(AtomImpl):
    _atomic_number = 26

    def clone(self) -> Fe:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Fe:
        return self.clone()._with_id(id)


class Co(AtomImpl):
    _atomic_number = 27

    def clone(self) -> Co:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Co:
        return self.clone()._with_id(id)


class Ni(AtomImpl):
    _atomic_number = 28

    def clone(self) -> Ni:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Ni:
        return self.clone()._with_id(id)


class Cu(AtomImpl):
    _atomic_number = 29

    def clone(self) -> Cu:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Cu:
        return self.clone()._with_id(id)


class Zn(AtomImpl):
    _atomic_number = 30

    def clone(self) -> Zn:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Zn:
        return self.clone()._with_id(id)


class Ga(AtomImpl):
    _atomic_number = 31

    def clone(self) -> Ga:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Ga:
        return self.clone()._with_id(id)


class Ge(AtomImpl):
    _atomic_number = 32

    def clone(self) -> Ge:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Ge:
        return self.clone()._with_id(id)


class As(AtomImpl):
    _atomic_number = 33

    def clone(self) -> As:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> As:
        return self.clone()._with_id(id)


class Se(AtomImpl):
    _atomic_number = 34

    def clone(self) -> Se:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Se:
        return self.clone()._with_id(id)


class Br(AtomImpl):
    _atomic_number = 35

    def clone(self) -> Br:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Br:
        return self.clone()._with_id(id)


class Kr(AtomImpl):
    _atomic_number = 36

    def clone(self) -> Kr:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Kr:
        return self.clone()._with_id(id)


class Rb(AtomImpl):
    _atomic_number = 37

    def clone(self) -> Rb:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Rb:
        return self.clone()._with_id(id)


class Sr(AtomImpl):
    _atomic_number = 38

    def clone(self) -> Sr:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Sr:
        return self.clone()._with_id(id)


class Y(AtomImpl):
    _atomic_number = 39

    def clone(self) -> Y:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Y:
        return self.clone()._with_id(id)


class Zr(AtomImpl):
    _atomic_number = 40

    def clone(self) -> Zr:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Zr:
        return self.clone()._with_id(id)


class Nb(AtomImpl):
    _atomic_number = 41

    def clone(self) -> Nb:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Nb:
        return self.clone()._with_id(id)


class Mo(AtomImpl):
    _atomic_number = 42

    def clone(self) -> Mo:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Mo:
        return self.clone()._with_id(id)


class Tc(AtomImpl):
    _atomic_number = 43

    def clone(self) -> Tc:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Tc:
        return self.clone()._with_id(id)


class Ru(AtomImpl):
    _atomic_number = 44

    def clone(self) -> Ru:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Ru:
        return self.clone()._with_id(id)


class Rh(AtomImpl):
    _atomic_number = 45

    def clone(self) -> Rh:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Rh:
        return self.clone()._with_id(id)


class Pd(AtomImpl):
    _atomic_number = 46

    def clone(self) -> Pd:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Pd:
        return self.clone()._with_id(id)


class Ag(AtomImpl):
    _atomic_number = 47

    def clone(self) -> Ag:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Ag:
        return self.clone()._with_id(id)


class Cd(AtomImpl):
    _atomic_number = 48

    def clone(self) -> Cd:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Cd:
        return self.clone()._with_id(id)


class In(AtomImpl):
    _atomic_number = 49

    def clone(self) -> In:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> In:
        return self.clone()._with_id(id)


class Sn(AtomImpl):
    _atomic_number = 50

    def clone(self) -> Sn:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Sn:
        return self.clone()._with_id(id)


class Sb(AtomImpl):
    _atomic_number = 51

    def clone(self) -> Sb:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Sb:
        return self.clone()._with_id(id)


class Te(AtomImpl):
    _atomic_number = 52

    def clone(self) -> Te:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Te:
        return self.clone()._with_id(id)


# "I" is a valid elemental symbol.
class I(AtomImpl):  # noqa
    _atomic_number = 53

    def clone(self) -> I:  # noqa
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> I:  # noqa
        return self.clone()._with_id(id)


class Xe(AtomImpl):
    _atomic_number = 54

    def clone(self) -> Xe:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Xe:
        return self.clone()._with_id(id)


class Cs(AtomImpl):
    _atomic_number = 55

    def clone(self) -> Cs:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Cs:
        return self.clone()._with_id(id)


class Ba(AtomImpl):
    _atomic_number = 56

    def clone(self) -> Ba:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Ba:
        return self.clone()._with_id(id)


class La(AtomImpl):
    _atomic_number = 57

    def clone(self) -> La:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> La:
        return self.clone()._with_id(id)


class Ce(AtomImpl):
    _atomic_number = 58

    def clone(self) -> Ce:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Ce:
        return self.clone()._with_id(id)


class Pr(AtomImpl):
    _atomic_number = 59

    def clone(self) -> Pr:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Pr:
        return self.clone()._with_id(id)


class Nd(AtomImpl):
    _atomic_number = 60

    def clone(self) -> Nd:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Nd:
        return self.clone()._with_id(id)


class Pm(AtomImpl):
    _atomic_number = 61

    def clone(self) -> Pm:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Pm:
        return self.clone()._with_id(id)


class Sm(AtomImpl):
    _atomic_number = 62

    def clone(self) -> Sm:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Sm:
        return self.clone()._with_id(id)


class Eu(AtomImpl):
    _atomic_number = 63

    def clone(self) -> Eu:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Eu:
        return self.clone()._with_id(id)


class Gd(AtomImpl):
    _atomic_number = 64

    def clone(self) -> Gd:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Gd:
        return self.clone()._with_id(id)


class Tb(AtomImpl):
    _atomic_number = 65

    def clone(self) -> Tb:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Tb:
        return self.clone()._with_id(id)


class Dy(AtomImpl):
    _atomic_number = 66

    def clone(self) -> Dy:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Dy:
        return self.clone()._with_id(id)


class Ho(AtomImpl):
    _atomic_number = 67

    def clone(self) -> Ho:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Ho:
        return self.clone()._with_id(id)


class Er(AtomImpl):
    _atomic_number = 68

    def clone(self) -> Er:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Er:
        return self.clone()._with_id(id)


class Tm(AtomImpl):
    _atomic_number = 69

    def clone(self) -> Tm:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Tm:
        return self.clone()._with_id(id)


class Yb(AtomImpl):
    _atomic_number = 70

    def clone(self) -> Yb:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Yb:
        return self.clone()._with_id(id)


class Lu(AtomImpl):
    _atomic_number = 71

    def clone(self) -> Lu:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Lu:
        return self.clone()._with_id(id)


class Hf(AtomImpl):
    _atomic_number = 72

    def clone(self) -> Hf:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Hf:
        return self.clone()._with_id(id)


class Ta(AtomImpl):
    _atomic_number = 73

    def clone(self) -> Ta:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Ta:
        return self.clone()._with_id(id)


class W(AtomImpl):
    _atomic_number = 74

    def clone(self) -> W:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> W:
        return self.clone()._with_id(id)


class Re(AtomImpl):
    _atomic_number = 75

    def clone(self) -> Re:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Re:
        return self.clone()._with_id(id)


class Os(AtomImpl):
    _atomic_number = 76

    def clone(self) -> Os:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Os:
        return self.clone()._with_id(id)


class Ir(AtomImpl):
    _atomic_number = 77

    def clone(self) -> Ir:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Ir:
        return self.clone()._with_id(id)


class Pt(AtomImpl):
    _atomic_number = 78

    def clone(self) -> Pt:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Pt:
        return self.clone()._with_id(id)


class Au(AtomImpl):
    _atomic_number = 79

    def clone(self) -> Au:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Au:
        return self.clone()._with_id(id)


class Hg(AtomImpl):
    _atomic_number = 80

    def clone(self) -> Hg:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Hg:
        return self.clone()._with_id(id)


class Tl(AtomImpl):
    _atomic_number = 81

    def clone(self) -> Tl:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Tl:
        return self.clone()._with_id(id)


class Pb(AtomImpl):
    _atomic_number = 82

    def clone(self) -> Pb:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Pb:
        return self.clone()._with_id(id)


class Bi(AtomImpl):
    _atomic_number = 83

    def clone(self) -> Bi:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Bi:
        return self.clone()._with_id(id)


class Po(AtomImpl):
    _atomic_number = 84

    def clone(self) -> Po:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Po:
        return self.clone()._with_id(id)


class At(AtomImpl):
    _atomic_number = 85

    def clone(self) -> At:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> At:
        return self.clone()._with_id(id)


class Rn(AtomImpl):
    _atomic_number = 86

    def clone(self) -> Rn:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Rn:
        return self.clone()._with_id(id)


class Fr(AtomImpl):
    _atomic_number = 87

    def clone(self) -> Fr:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Fr:
        return self.clone()._with_id(id)


class Ra(AtomImpl):
    _atomic_number = 88

    def clone(self) -> Ra:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Ra:
        return self.clone()._with_id(id)


class Ac(AtomImpl):
    _atomic_number = 89

    def clone(self) -> Ac:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Ac:
        return self.clone()._with_id(id)


class Th(AtomImpl):
    _atomic_number = 90

    def clone(self) -> Th:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Th:
        return self.clone()._with_id(id)


class Pa(AtomImpl):
    _atomic_number = 91

    def clone(self) -> Pa:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Pa:
        return self.clone()._with_id(id)


class U(AtomImpl):
    _atomic_number = 92

    def clone(self) -> U:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> U:
        return self.clone()._with_id(id)


class Np(AtomImpl):
    _atomic_number = 93

    def clone(self) -> Np:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Np:
        return self.clone()._with_id(id)


class Pu(AtomImpl):
    _atomic_number = 94

    def clone(self) -> Pu:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Pu:
        return self.clone()._with_id(id)


class Am(AtomImpl):
    _atomic_number = 95

    def clone(self) -> Am:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Am:
        return self.clone()._with_id(id)


class Cm(AtomImpl):
    _atomic_number = 96

    def clone(self) -> Cm:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Cm:
        return self.clone()._with_id(id)


class Bk(AtomImpl):
    _atomic_number = 97

    def clone(self) -> Bk:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Bk:
        return self.clone()._with_id(id)


class Cf(AtomImpl):
    _atomic_number = 98

    def clone(self) -> Cf:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Cf:
        return self.clone()._with_id(id)


class Es(AtomImpl):
    _atomic_number = 99

    def clone(self) -> Es:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Es:
        return self.clone()._with_id(id)


class Fm(AtomImpl):
    _atomic_number = 100

    def clone(self) -> Fm:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Fm:
        return self.clone()._with_id(id)


class Md(AtomImpl):
    _atomic_number = 101

    def clone(self) -> Md:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Md:
        return self.clone()._with_id(id)


class No(AtomImpl):
    _atomic_number = 102

    def clone(self) -> No:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> No:
        return self.clone()._with_id(id)


class Lr(AtomImpl):
    _atomic_number = 103

    def clone(self) -> Lr:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Lr:
        return self.clone()._with_id(id)


class Rf(AtomImpl):
    _atomic_number = 104

    def clone(self) -> Rf:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Rf:
        return self.clone()._with_id(id)


class Db(AtomImpl):
    _atomic_number = 105

    def clone(self) -> Db:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Db:
        return self.clone()._with_id(id)


class Sg(AtomImpl):
    _atomic_number = 106

    def clone(self) -> Sg:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Sg:
        return self.clone()._with_id(id)


class Bh(AtomImpl):
    _atomic_number = 107

    def clone(self) -> Bh:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Bh:
        return self.clone()._with_id(id)


class Hs(AtomImpl):
    _atomic_number = 108

    def clone(self) -> Hs:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Hs:
        return self.clone()._with_id(id)


class Mt(AtomImpl):
    _atomic_number = 109

    def clone(self) -> Mt:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Mt:
        return self.clone()._with_id(id)


class Ds(AtomImpl):
    _atomic_number = 110

    def clone(self) -> Ds:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Ds:
        return self.clone()._with_id(id)


class Rg(AtomImpl):
    _atomic_number = 111

    def clone(self) -> Rg:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Rg:
        return self.clone()._with_id(id)


class Cn(AtomImpl):
    _atomic_number = 112

    def clone(self) -> Cn:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Cn:
        return self.clone()._with_id(id)


class Nh(AtomImpl):
    _atomic_number = 113

    def clone(self) -> Nh:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Nh:
        return self.clone()._with_id(id)


class Fl(AtomImpl):
    _atomic_number = 114

    def clone(self) -> Fl:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Fl:
        return self.clone()._with_id(id)


class Mc(AtomImpl):
    _atomic_number = 115

    def clone(self) -> Mc:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Mc:
        return self.clone()._with_id(id)


class Lv(AtomImpl):
    _atomic_number = 116

    def clone(self) -> Lv:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Lv:
        return self.clone()._with_id(id)


class Ts(AtomImpl):
    _atomic_number = 117

    def clone(self) -> Ts:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Ts:
        return self.clone()._with_id(id)


class Og(AtomImpl):
    _atomic_number = 118

    def clone(self) -> Og:
        return type(self)(self._id, self._charge)

    def with_id(self, id: int) -> Og:
        return self.clone()._with_id(id)
