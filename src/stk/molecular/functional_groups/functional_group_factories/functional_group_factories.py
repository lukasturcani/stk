import itertools as it

from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from .. import FunctionalGroup, _FunctionalGroup


class Amine(_FunctionalGroup, FunctionalGroup):
    """
    Represents an amine functional group.

    """

    pass


class AmineFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Amine` instances.

    """

    _functional_group_smarts = '[N]([H])[H]'
    _bonder_smarts = ['[$([N]([H])[H])]']
    _deleter_smarts = ['[$([H][N][H])]']*2

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


class SecondaryAmineFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Amine` instances.

    """

    _functional_group = Amine
    _functional_group_smarts = '[H][N]([#6])[#6]'
    _bonder_smarts = ['[$([N]([H])([#6])[#6])]']
    _deleter_smarts = ['[$([H][N]([#6])[#6])]']


class Aldehyde(_FunctionalGroup, FunctionalGroup):
    """
    Represents an aldehyde functional group.

    """

    pass


class AldehydeFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Aldehyde` instances.

    """

    _functional_group = Aldehyde
    _functional_group_smarts = '[C](=[O])[H]'
    _bonder_smarts = ['[$([C](=[O])[H])]']
    _deleter_smarts = ['[$([O]=[C][H])]']


class CarboxylicAcid(_FunctionalGroup, FunctionalGroup):
    """
    Represents a carboxylic acid functional group.

    """

    pass


class CarboxylicAcidFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.CarboxylicAcid` instances.

    """

    _functional_group = CarboxylicAcid
    _functional_group_smarts = '[C](=[O])[O][H]'
    _bonder_smarts = ['[$([C](=[O])[O][H])]']
    _deleter_smarts = [
        '[$([H][O][C](=[O]))]',
        '[$([O]([H])[C](=[O]))]',
    ]


class Amide(_FunctionalGroup, FunctionalGroup):
    """
    Represents an amide functional group.

    """

    pass


class AmideFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Amide` instances.

    """

    _functional_group = Amide
    _functional_group_smarts = '[C](=[O])[N]([H])[H]'
    _bonder_smarts = ['[$([C](=[O])[N]([H])[H])]']
    _deleter_smarts = (
        ['[$([N]([H])([H])[C](=[O]))]'] +
        ['[$([H][N]([H])[C](=[O]))]']*2
    )


class Thioacid(_FunctionalGroup, FunctionalGroup):
    """
    Represents a thioacid functional group.

    """

    pass


class ThioacidFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Thioacid` instances.

    """

    _functional_group = Thioacid
    _functional_group_smarts = '[C](=[O])[S][H]'
    _bonder_smarts = ['[$([C](=[O])[S][H])]']
    _deleter_smarts = [
        '[$([H][S][C](=[O]))]',
        '[$([S]([H])[C](=[O]))]',
    ]


class Alcohol(_FunctionalGroup, FunctionalGroup):
    """
    Represents an alcohol functional group.

    """

    pass


class AlcoholFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Alcohol` instances.

    """

    _functional_group = Alcohol
    _functional_group_smarts = '[O][H]'
    _bonder_smarts = ['[$([O][H])]']
    _deleter_smarts = ['[$([H][O])]']


class Thiol(_FunctionalGroup, FunctionalGroup):
    """
    Represents a thiol functional group.

    """

    pass


class ThiolFactory(SmartsFunctionalGroupFactory):
    _functional_group = Thiol
    _functional_group_smarts = '[S][H]'
    _bonder_smarts = ['[$([S][H])]']
    _deleter_smarts = ['[$([H][S])]']


class Fluoro(_FunctionalGroup, FunctionalGroup):
    """
    Represents a fluoro functional group.

    """

    pass


class FluoroFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Fluoro` instances.

    """

    _functional_group = Fluoro
    _functional_group_smarts = '*[F]'
    _bonder_smarts = ['[$(*[F])]']
    _deleter_smarts = ['[$([F]*)]']


class Bromo(_FunctionalGroup, FunctionalGroup):
    """
    Represents a bromo functional group.

    """


class BromoFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Bromo` instances.

    """

    _functional_group = Bromo
    _functional_group_smarts = '*[Br]'
    _bonder_smarts = ['[$(*[Br])]']
    _deleter_smarts = ['[$([Br]*)]']


class Iodo(_FunctionalGroup, FunctionalGroup):
    """
    Represents an iodo functional group.

    """

    pass


class IodoFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Iodo` instances.

    """

    _functional_group = Iodo
    _functional_group_smarts = '*[I]'
    _bonder_smarts = ['[$(*[I])]']
    _deleter_smarts = ['[$([I]*)]']


class TerminalAlkyne(_FunctionalGroup, FunctionalGroup):
    """
    Represents an alkyne functional group.

    """

    pass


class TerminalAlkyneFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Alkyne` instances.

    """

    _functional_group = TerminalAlkyne
    _functional_group_smarts = '[C]#[C][H]'
    _bonder_smarts = ['[$([C]([H])#[C])]']
    _deleter_smarts = ['[$([H][C]#[C])]', '[$([C](#[C])[H])]']

    def __init__(self, delete_carbon=False):
        """
        Initialize a :class:`.TerminalAlkyneFactory`.

        Parameters
        ----------
        delete_carbon : :class:`bool`, optional
            If ``True``, the terminal carbon is also added to the
            deleter atoms.

        """

        self._delete_carbon = delete_carbon

    def get_functional_groups(self, molecule):
        for ids in self._get_ids(molecule):
            atoms = tuple(molecule.get_atoms(ids.atom_ids))
            bonder_ids = set(ids.bonder_ids)
            deleter_ids = set(ids.deleter_ids)
            yield TerminalAlkyne(
                atoms=atoms,
                bonders=tuple(
                    a for a in atoms if a.id in bonder_ids
                ),
                deleters=tuple(
                    a for a in atoms
                    if a.id in deleter_ids and self._is_deleter(a)
                ),
            )

    def _is_deleter(self, atom):
        if atom.atomic_number == 1:
            return True
        elif self._delete_carbon:
            return True
        else:
            return False

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'delete_carbon={self._delete_carbon})'
        )


class TerminalAlkene(_FunctionalGroup, FunctionalGroup):
    """
    Represents a terminal alkene functional group.

    """

    pass


class TerminalAlkeneFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.TerminalAlkene` instances.

    """

    _functional_group = TerminalAlkene
    _functional_group_smarts = '[C]=[C]([H])[H]'
    _bonder_smarts = ['[$([C]=[C]([H])[H])]']
    _deleter_smarts = (
        ['[$([H][C]([H])=[C])]']*2 +
        ['[$([C](=[C])([H])[H])]']
    )


class BoronicAcid(_FunctionalGroup, FunctionalGroup):
    """
    Represents a boronic acid functional group.

    """

    pass


class BoronicAcidFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.BoronicAcid` instances.

    """

    _functional_group = BoronicAcid
    _functional_group_smarts = '[B]([O][H])[O][H]'
    _bonder_smarts = ['[$([B]([O][H])[O][H])]']
    _deleter_smarts = (
        ['[$([O]([H])[B][O][H])]']*2 +
        ['[$([H][O][B][O][H])]']*2
    )


class Diol(_FunctionalGroup, FunctionalGroup):
    """
    Represents a diol functional group.

    """

    pass


class DiolFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Diol` instances.

    """

    _functional_group = Diol
    _functional_group_smarts = '[H][O][#6]~[#6][O][H]'
    _bonder_smarts = ['[$([O]([H])[#6]~[#6][O][H])]']*2
    _deleter_smarts = ['[$([H][O][#6]~[#6][O][H])]']*2


class Difluoro(_FunctionalGroup, FunctionalGroup):
    """
    Represents a difluoro functional group.

    """

    pass


class DifluoroFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Difluoro` instances.

    """

    _functional_group = Difluoro
    _functional_group_smarts = '[F][#6]~[#6][F]'
    _bonder_smarts = ['[$([#6]([F])~[#6][F])]']*2
    _deleter_smarts = ['[$([F][#6]~[#6][F])]']*2


class Dibromo(_FunctionalGroup, FunctionalGroup):
    """
    Represents a dibromo functional group.

    """

    pass


class DibromoFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Dibromo` instances.

    """

    _functional_group = Dibromo
    _functional_group_smarts = '[Br][#6]~[#6][Br]'
    _bonder_smarts = ['[$([#6]([Br])~[#6][Br])]']*2
    _deleter_smarts = ['[$([Br][#6]~[#6][Br])]']*2


class RingAmine(_FunctionalGroup, FunctionalGroup):
    """
    Represents an amine bonded to a ring.

    """

    pass


class RingAmineFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.RingAmine` functional groups.

    """

    _functional_group = RingAmine
    _functional_group_smarts = '[N]([H])([H])[#6]~[#6]([H])~[#6R1]'
    _bonder_smarts = [
        '[$([N]([H])([H])[#6]~[#6]([H])~[#6R1])]',
        '[$([#6]([H])(~[#6R1])~[#6][N]([H])[H])]',
    ]
    _deleter_smarts = (
            ['[$([H][N]([H])[#6]~[#6]([H])~[#6R1])]']*2 +
            ['[$([H][#6](~[#6R1])~[#6][N]([H])[H])]']
    )
