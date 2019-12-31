from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from .. import FunctionalGroup_


class PrimaryAmine(FunctionalGroup_):
    """
    Represents a primary amine functional group.

    """

    def __init__(
        self,
        nitrogen,
        hydrogen1,
        hydrogen2,
        r,
        bonders,
        deleters,
    ):
        self._nitrogen = nitrogen
        self._hydrogen1 = hydrogen1
        self._hydrogen2 = hydrogen2
        self._r = r
        atoms = (nitrogen, hydrogen1, hydrogen2, r)
        super().__init__(atoms, bonders, deleters)

    def get_nitrogen(self):
        return self._nitrogen.clone()

    def get_hydrogen1(self):
        return self._hydrogen1.clone()

    def get_hydrogen2(self):
        return self._hydrogen2.clone()

    def get_r(self):
        return self._r.clone()

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._nitrogen}, {self._hydrogen1}, {self._hydrogen2}, '
            f'{self._r}, bonders={self._bonders}, '
            f'deleters={self._deleters}'
            ')'
        )


class PrimaryAmineFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.PrimaryAmine` instances.

    """

    _functional_group_smarts = '[*][N]([H])[H]'
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
            atoms = tuple(molecule.get_atoms(ids.atom_ids))
            yield PrimaryAmine(
                nitrogen=atoms[1],
                hydrogen1=atoms[2],
                hydrogen2=atoms[3],
                r=atoms[0],
                bonders=(atoms[1], ),
                deleters=atoms[2:2+self._num_deleters],
            )

    def __repr__(self):
        return (
            f'{self.__class__.__name__}'
            f'(num_deleters={self._num_deleters})'
        )


class SecondaryAmine(FunctionalGroup_):
    """
    Represents a secondary amine functional group.

    """

    def __init__(self, nitrogen, hydrogen, r1, r2, bonders, deleters):
        self._nitrogen = nitrogen
        self._hydrogen = hydrogen
        self._r1 = r1
        self._r2 = r2
        atoms = (nitrogen, hydrogen, r1, r2)
        super().__init__(atoms, bonders, deleters)

    def get_nitrogen(self):
        return self._nitrogen.clone()

    def get_hydrogen(self):
        return self._hydrogen.clone()

    def get_r1(self):
        return self._r1.clone()

    def get_r2(self):
        return self._r2.clone()


class SecondaryAmineFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.SecondaryAmine` instances.

    """

    _functional_group = SecondaryAmine
    _functional_group_smarts = '[H][N]([#6])[#6]'
    _bonder_smarts = ['[$([N]([H])([#6])[#6])]']
    _deleter_smarts = ['[$([H][N]([#6])[#6])]']

    def get_functional_groups(self, molecule):
        for ids in self._get_ids(molecule):
            atoms = tuple(molecule.get_atoms(ids.atom_ids))
            yield SecondaryAmine(
                nitrogen=atoms[1],
                hydrogen=atoms[0],
                r1=atoms[2],
                r2=atoms[3],
                bonders=tuple(molecule.get_atoms(ids.bonder_ids)),
                deleters=tuple(molecule.get_atoms(ids.deleter_ids)),
            )


class Aldehyde(FunctionalGroup_):
    """
    Represents an aldehyde functional group.

    """

    def __init__(self, carbon, oxygen, hydrogen, r, bonders, deleters):
        self._carbon = carbon
        self._oxygen = oxygen
        self._hydrogen = hydrogen
        self._r = r
        atoms = (carbon, oxygen, hydrogen, r)
        super().__init__(atoms, bonders, deleters)

    def get_carbon(self):
        return self._carbon.clone()

    def get_oxygen(self):
        return self._oxygen.clone()

    def get_hydrogen(self):
        return self._hydrogen.clone()

    def get_r(self):
        return self._r.clone()


class AldehydeFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Aldehyde` instances.

    """

    _functional_group = Aldehyde
    _functional_group_smarts = '[*][C](=[O])[H]'
    _bonder_smarts = ['[$([C](=[O])[H])]']
    _deleter_smarts = ['[$([O]=[C][H])]']

    def get_functional_groups(self, molecule):
        for ids in self._get_ids(molecule):
            atoms = tuple(molecule.get_atoms(ids.atom_ids))
            yield Aldehyde(
                carbon=atoms[1],
                oxygen=atoms[2],
                hydrogen=atoms[3],
                r=atoms[0],
                bonders=tuple(molecule.get_atoms(ids.bonder_ids)),
                deleters=tuple(molecule.get_atoms(ids.deleter_ids)),
            )


class CarboxylicAcid(FunctionalGroup_):
    """
    Represents a carboxylic acid functional group.

    """

    def __init__(
        self,
        carbon,
        oxygen1,
        oxygen2,
        hydrogen,
        r,
        bonders,
        deleters,
    ):
        self._carbon = carbon
        self._oxygen1 = oxygen1
        self._oxygen2 = oxygen2
        self._hydrogen = hydrogen
        self._r = r
        atoms = (carbon, oxygen1, oxygen2, hydrogen, r)
        super().__init__(atoms, bonders, deleters)

    def get_carbon(self):
        return self._carbon.clone()

    def get_oxygen1(self):
        return self._oxygen1.clone()

    def get_oxygen2(self):
        return self._oxygen2.clone()

    def get_hydrogen(self):
        return self._hydrogen.clone()

    def get_r(self):
        return self._r.clone()


class CarboxylicAcidFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.CarboxylicAcid` instances.

    """

    _functional_group = CarboxylicAcid
    _functional_group_smarts = '[*][C](=[O])[O][H]'
    _bonder_smarts = ['[$([C](=[O])[O][H])]']
    _deleter_smarts = [
        '[$([H][O][C](=[O]))]',
        '[$([O]([H])[C](=[O]))]',
    ]

    def get_functional_groups(self, molecule):
        for ids in self._get_ids(molecule):
            atoms = tuple(molecule.get_atoms(ids.atom_ids))
            yield CarboxylicAcid(
                carbon=atoms[1],
                oxygen1=atoms[2],
                oxygen2=atoms[3],
                hydrogen=atoms[4],
                r=atoms[0],
                bonders=tuple(molecule.get_atoms(ids.bonder_ids)),
                deleters=tuple(molecule.get_atoms(ids.deleter_ids)),
            )


class Amide(FunctionalGroup_):
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


class Thioacid(FunctionalGroup_):
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


class Alcohol(FunctionalGroup_):
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


class Thiol(FunctionalGroup_):
    """
    Represents a thiol functional group.

    """

    pass


class ThiolFactory(SmartsFunctionalGroupFactory):
    _functional_group = Thiol
    _functional_group_smarts = '[S][H]'
    _bonder_smarts = ['[$([S][H])]']
    _deleter_smarts = ['[$([H][S])]']


class Fluoro(FunctionalGroup_):
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


class Bromo(FunctionalGroup_):
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


class Iodo(FunctionalGroup_):
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


class TerminalAlkyne(FunctionalGroup_):
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
    _bonder_smarts = ()
    _deleter_smarts = ()

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
        if delete_carbon:
            self._bonder_smarts = ['[$([C]#[C][H])]']
            self._deleter_smarts = [
                '[$([H][C]#[C])]',
                '[$([C]([H])#[C])]',
            ]
        else:
            self._bonder_smarts = ['[$([C]([H])#[C])]']
            self._deleter_smarts = ['[$([H][C]#[C])]']

        self._set_queries(self)

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


class TerminalAlkene(FunctionalGroup_):
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


class BoronicAcid(FunctionalGroup_):
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


class Diol(FunctionalGroup_):
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


class Difluoro(FunctionalGroup_):
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


class Dibromo(FunctionalGroup_):
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


class RingAmine(FunctionalGroup_):
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
