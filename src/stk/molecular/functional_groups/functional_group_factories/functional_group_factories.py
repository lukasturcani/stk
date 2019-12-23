import itertools as it

from .. import FunctionalGroup, _FunctionalGroup
from . import SmartsFunctionalGroupFactory


class Amine(FunctionalGroup, _FunctionalGroup):
    """
    Represents an amine functional group.

    """

    pass


class AmineFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Amine` instances.

    """

    functional_group_smarts = '[N]([H])[H]'
    bonder_smarts = '[$([N]([H])[H])]'
    deleter_smarts = '[$([H][N][H])]'

    def __init__(self, num_deleters=2):
        """
        Initialize an :class:`.AmineFactory`.

        Parameters
        ----------
        num_deleters : :class:`int`, optional
            The number of deleter atoms the created :class:`.Amine`
            instances will have, maximum is ``2``. If ``0`` or ``1``
            is used, the :class:`.Amine` instances created by the
            factory will lose ``0`` or ``1`` hydrogen atoms during
            construction, respectively.

        """

        self._num_deleters = num_deleters

    def get_functional_groups(self, molecule):
        for ids in self._get_ids(molecule):
            atoms = tuple(molecule.get_atoms(
                atom_ids=it.chain(ids.bonder_ids, ids.deleter_ids)
            ))
            yield Amine(
                atoms=atoms,
                bonders=(atoms[0], ),
                deleters=atoms[1:1+self._num_deleters],
            )

    def __repr__(self):
        return (
            f'{self.__class__.__name__}'
            f'(num_deleters={self._num_deleters})'
        )


class PrimaryAmine(FunctionalGroup, _FunctionalGroup):
    pass


class PrimaryAmineFactory(SmartsFunctionalGroupFactory):
    pass


class Aldehyde(FunctionalGroup, _FunctionalGroup):
    def __init__(self):
        pass


class AldehydeFactory(SmartsFunctionalGroupFactory):
    pass


class CarboxylicAcid(FunctionalGroup, _FunctionalGroup):
    pass


class CarboxylicAcidFactory(SmartsFunctionalGroupFactory):
    pass


class Amide(FunctionalGroup, _FunctionalGroup):
    pass


class AmideFactory(SmartsFunctionalGroupFactory):
    pass


class Thioacid(FunctionalGroup, _FunctionalGroup):
    pass


class ThioacidFactory(SmartsFunctionalGroupFactory):
    pass


class Bromine(FunctionalGroup, _FunctionalGroup):
    pass


class BromineFactory(SmartsFunctionalGroupFactory):
    pass


class Iodine(FunctionalGroup, _FunctionalGroup):
    pass


class IodineFactory(SmartsFunctionalGroupFactory):
    pass


class Alkyne(FunctionalGroup, _FunctionalGroup):
    pass


class AlkyneFactory(SmartsFunctionalGroupFactory):
    pass


class TerminalAlkene(FunctionalGroup, _FunctionalGroup):
    pass


class TerminalAlkeneFactory(SmartsFunctionalGroupFactory):
    pass


class BoronicAcid(FunctionalGroup, _FunctionalGroup):
    pass


class BoronicAcidFactory(SmartsFunctionalGroupFactory):
    pass


class SecondaryAmine(FunctionalGroup, _FunctionalGroup):
    pass


class SecondaryAmineFactory(SmartsFunctionalGroupFactory):
    pass


class Diol(FunctionalGroup, _FunctionalGroup):
    pass


class DiolFactory(SmartsFunctionalGroupFactory):
    pass


class Difluorene(FunctionalGroup, _FunctionalGroup):
    pass


class DifluoreneFactory(SmartsFunctionalGroupFactory):
    pass


class Dibromine(FunctionalGroup, _FunctionalGroup):
    pass


class DibromineFactory(SmartsFunctionalGroupFactory):
    pass


class RingAmine(FunctionalGroup, _FunctionalGroup):
    pass


class RingAmineFactory(SmartsFunctionalGroupFactory):
    pass
