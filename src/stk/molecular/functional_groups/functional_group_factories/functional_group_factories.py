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
        atom,
        bonders,
        deleters,
    ):
        self._nitrogen = nitrogen
        self._hydrogen1 = hydrogen1
        self._hydrogen2 = hydrogen2
        self._atom = atom
        atoms = (nitrogen, hydrogen1, hydrogen2, atom)
        super().__init__(atoms, bonders, deleters)

    def get_nitrogen(self):
        return self._nitrogen.clone()

    def get_hydrogen1(self):
        return self._hydrogen1.clone()

    def get_hydrogen2(self):
        return self._hydrogen2.clone()

    def get_atom(self):
        return self._atom.clone()

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._nitrogen}, {self._hydrogen1}, {self._hydrogen2}, '
            f'{self._atom}, bonders={self._bonders}, '
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
                atom=atoms[0],
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

    def __init__(
        self,
        nitrogen,
        hydrogen,
        atom1,
        atom2,
        bonders,
        deleters,
    ):
        self._nitrogen = nitrogen
        self._hydrogen = hydrogen
        self._atom1 = atom1
        self._atom2 = atom2
        atoms = (nitrogen, hydrogen, atom1, atom2)
        super().__init__(atoms, bonders, deleters)

    def get_nitrogen(self):
        return self._nitrogen.clone()

    def get_hydrogen(self):
        return self._hydrogen.clone()

    def get_atom1(self):
        return self._atom1.clone()

    def get_atom2(self):
        return self._atom2.clone()


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
                atom1=atoms[2],
                atom2=atoms[3],
                bonders=tuple(molecule.get_atoms(ids.bonder_ids)),
                deleters=tuple(molecule.get_atoms(ids.deleter_ids)),
            )


class Aldehyde(FunctionalGroup_):
    """
    Represents an aldehyde functional group.

    """

    def __init__(
        self,
        carbon,
        oxygen,
        hydrogen,
        atom,
        bonders,
        deleters,
    ):
        self._carbon = carbon
        self._oxygen = oxygen
        self._hydrogen = hydrogen
        self._atom = atom
        atoms = (carbon, oxygen, hydrogen, atom)
        super().__init__(atoms, bonders, deleters)

    def get_carbon(self):
        return self._carbon.clone()

    def get_oxygen(self):
        return self._oxygen.clone()

    def get_hydrogen(self):
        return self._hydrogen.clone()

    def get_atom(self):
        return self._atom.clone()


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
                atom=atoms[0],
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
        atom,
        bonders,
        deleters,
    ):
        self._carbon = carbon
        self._oxygen1 = oxygen1
        self._oxygen2 = oxygen2
        self._hydrogen = hydrogen
        self._atom = atom
        atoms = (carbon, oxygen1, oxygen2, hydrogen, atom)
        super().__init__(atoms, bonders, deleters)

    def get_carbon(self):
        return self._carbon.clone()

    def get_oxygen1(self):
        return self._oxygen1.clone()

    def get_oxygen2(self):
        return self._oxygen2.clone()

    def get_hydrogen(self):
        return self._hydrogen.clone()

    def get_atom(self):
        return self._atom.clone()


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
                atom=atoms[0],
                bonders=tuple(molecule.get_atoms(ids.bonder_ids)),
                deleters=tuple(molecule.get_atoms(ids.deleter_ids)),
            )


class Amide(FunctionalGroup_):
    """
    Represents an amide functional group.

    """

    def __init__(
        self,
        carbon,
        oxygen,
        nitrogen,
        hydrogen1,
        hydrogen2,
        atom,
        bonders,
        deleters,
    ):
        self._carbon = carbon
        self._oxygen = oxygen
        self._nitrogen = nitrogen
        self._hydrogen1 = hydrogen1
        self._hydrogen2 = hydrogen2
        self._atom = atom
        atoms = (carbon, oxygen, nitrogen, hydrogen1, hydrogen2, atom)
        super().__init__(atoms, bonders, deleters)

    def get_carbon(self):
        return self._carbon.clone()

    def get_oxygen(self):
        return self._oxygen.clone()

    def get_nitrogen(self):
        return self._nitrogen.clone()

    def get_hydrogen1(self):
        return self._hydrogen1.clone()

    def get_hydrogen2(self):
        return self._hydrogen2.clone()

    def get_atom(self):
        return self._atom.clone()


class AmideFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Amide` instances.

    """

    _functional_group = Amide
    _functional_group_smarts = '[*][C](=[O])[N]([H])[H]'
    _bonder_smarts = ['[$([C](=[O])[N]([H])[H])]']
    _deleter_smarts = (
        ['[$([N]([H])([H])[C](=[O]))]'] +
        ['[$([H][N]([H])[C](=[O]))]']*2
    )

    def get_functional_groups(self, molecule):
        for ids in self._get_ids(molecule):
            atoms = tuple(molecule.get_atoms(ids.atom_ids))
            yield Amide(
                carbon=atoms[1],
                oxygen=atoms[2],
                nitrogen=atoms[3],
                hydrogen1=atoms[4],
                hydrogen2=atoms[5],
                atom=atoms[0],
                bonders=tuple(molecule.get_atoms(ids.bonder_ids)),
                deleters=tuple(molecule.get_atoms(ids.deleter_ids)),
            )


class Thioacid(FunctionalGroup_):
    """
    Represents a thioacid functional group.

    """

    def __init__(
        self,
        carbon,
        oxygen,
        sulfur,
        hydrogen,
        atom,
        bonders,
        deleters
    ):
        self._carbon = carbon
        self._oxygen = oxygen
        self._sulfur = sulfur
        self._hydrogen = hydrogen
        self._atom = atom
        atoms = (carbon, oxygen, sulfur, hydrogen, atom)
        super().__init__(atoms, bonders, deleters)

    def get_carbon(self):
        return self._carbon.clone()

    def get_oxygen(self):
        return self._oxygen.clone()

    def get_sulfur(self):
        return self._sulfur.clone()

    def get_hydrogen(self):
        return self._hydrogen.clone()

    def get_atom(self):
        return self._atom.clone()


class ThioacidFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Thioacid` instances.

    """

    _functional_group = Thioacid
    _functional_group_smarts = '[*][C](=[O])[S][H]'
    _bonder_smarts = ['[$([C](=[O])[S][H])]']
    _deleter_smarts = [
        '[$([H][S][C](=[O]))]',
        '[$([S]([H])[C](=[O]))]',
    ]

    def get_functional_groups(self, molecule):
        for ids in self._get_ids(molecule):
            atoms = tuple(molecule.get_atoms(ids.atom_ids))
            yield Thioacid(
                carbon=atoms[1],
                oxygen=atoms[2],
                sulfur=atoms[3],
                hydrogen=atoms[4],
                atom=atoms[0],
                bonders=tuple(molecule.get_atoms(ids.bonder_ids)),
                deleters=tuple(molecule.get_atoms(ids.deleter_ids)),
            )


class Alcohol(FunctionalGroup_):
    """
    Represents an alcohol functional group.

    """

    def __init__(self, oxygen, hydrogen, atom, bonders, deleters):
        self._oxygen = oxygen
        self._hydrogen = hydrogen
        self._atom = atom
        atoms = (oxygen, hydrogen, atom)
        super().__init__(atoms, bonders, deleters)

    def get_oxygen(self):
        return self._oxygen.clone()

    def get_hydrogen(self):
        return self._hydrogen.clone()

    def get_atom(self):
        return self._atom.clone()


class AlcoholFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Alcohol` instances.

    """

    _functional_group = Alcohol
    _functional_group_smarts = '[*][O][H]'
    _bonder_smarts = ['[$([O][H])]']
    _deleter_smarts = ['[$([H][O])]']


class Thiol(FunctionalGroup_):
    """
    Represents a thiol functional group.

    """

    def __init__(self, sulfur, hydrogen, atom, bonders, deleters):
        self._sulfur = sulfur
        self._hydrogen = hydrogen
        self._atom = atom
        atoms = (sulfur, hydrogen, atom)
        super().__init__(atoms, bonders, deleters)

    def get_sulfur(self):
        return self._sulfur.clone()

    def get_hydrogen(self):
        return self._hydrogen.clone()

    def get_atom(self):
        return self._atom.clone()


class ThiolFactory(SmartsFunctionalGroupFactory):
    _functional_group = Thiol
    _functional_group_smarts = '[*][S][H]'
    _bonder_smarts = ['[$([S][H])]']
    _deleter_smarts = ['[$([H][S])]']

    def get_functional_groups(self, molecule):
        for ids in self._get_ids(molecule):
            atoms = tuple(molecule.get_atoms(ids.atom_ids))
            yield Thioacid(
                carbon=atoms[1],
                oxygen=atoms[2],
                sulfur=atoms[3],
                hydrogen=atoms[4],
                atom=atoms[0],
                bonders=tuple(molecule.get_atoms(ids.bonder_ids)),
                deleters=tuple(molecule.get_atoms(ids.deleter_ids)),
            )


class Fluoro(FunctionalGroup_):
    """
    Represents a fluoro functional group.

    """

    def __init__(self, fluorine, atom, bonders, deleters):
        self._fluorine = fluorine
        self._atom = atom
        super().__init__((fluorine, atom), bonders, deleters)

    def get_fluorine(self):
        return self._fluorine.clone()

    def get_atom(self):
        return self._atom.clone()


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

    def __init__(self, bromine, atom, bonders, deleters):
        self._bromine = bromine
        self._atom = atom
        super().__init__((bromine, atom), bonders, deleters)

    def get_bromine(self):
        return self._bromine.clone()

    def get_atom(self):
        return self._atom.clone()


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

    def __init__(self, iodine, atom, bonders, deleters):
        self._iodine = iodine
        self._atom = atom
        super().__init__((iodine, atom), bonders, deleters)

    def get_iodine(self):
        return self._iodine.clone()

    def get_atom(self):
        return self._atom.clone()


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

    def __init__(
        self,
        carbon1,
        carbon2,
        hydrogen,
        atom,
        bonders,
        deleters,
    ):
        self._carbon1 = carbon1
        self._carbon2 = carbon2
        self._hydrogen = hydrogen
        self._atom = atom
        atoms = (carbon1, carbon2, hydrogen, atom)
        super().__init__(atoms, bonders, deleters)

    def get_carbon1(self):
        return self._carbon1.clone()

    def get_carbon2(self):
        return self._carbon2.clone()

    def get_hydrogen(self):
        return self._hydrogen.clone()

    def get_atom(self):
        return self._atom.clone()


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

    def __init__(
        self,
        carbon1,
        carbon2,
        hydrogen1,
        hydrogen2,
        atom,
        bonders,
        deleters,
    ):
        self._carbon1 = carbon1
        self._carbon2 = carbon2
        self._hydrogen1 = hydrogen1
        self._hydrogen2 = hydrogen2
        self._atom = atom
        atoms = (carbon1, carbon2, hydrogen1, hydrogen2, atom)
        super().__init__(atoms, bonders, deleters)

    def get_carbon1(self):
        return self._carbon1.clone()

    def get_carbon2(self):
        return self._carbon2.clone()

    def get_hydrogen1(self):
        return self._hydrogen1.clone()

    def get_hydrogen2(self):
        return self._hydrogen2.clone()

    def get_atom(self):
        return self._atom.clone()


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

    def __init__(
        self,
        boron,
        oxygen1,
        hydrogen1,
        oxygen2,
        hydrogen2,
        atom,
        bonders,
        deleters,
    ):
        self._boron = boron
        self._oxygen1 = oxygen1
        self._hydrogen1 = hydrogen1
        self._oxygen2 = oxygen2
        self._hydrogen2 = hydrogen2
        self._atom = atom
        atoms = (boron, oxygen1, hydrogen1, oxygen2, hydrogen2, atom)
        super().__init__(atoms, bonders, deleters)

    def get_boron(self):
        return self._boron.clone()

    def get_oxygen1(self):
        return self._oxygen1.clone()

    def get_hydrogen1(self):
        return self._hydrogen1.clone()

    def get_oxygen2(self):
        return self._oxygen2.clone()

    def get_hydrogen2(self):
        return self._hydrogen2.clone()

    def get_atom(self):
        return self._atom.clone()


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

    def __init__(
        self,
        atom1,
        oxygen1,
        hydrogen1,
        atom2,
        oxygen2,
        hydrogen2,
        atom3,
        atom4,
        bonders,
        deleters,
    ):
        self._atom1 = atom1
        self._oxygen1 = oxygen1
        self._hydrogen1 = hydrogen1
        self._atom2 = atom2
        self._oxygen2 = oxygen2
        self._hydrogen2 = hydrogen2
        self._atom3 = atom3
        self._atom4 = atom4
        atoms = (
            atom1,
            oxygen1,
            hydrogen1,
            atom2,
            oxygen2,
            hydrogen2,
            atom3,
            atom4,
        )
        super().__init__(atoms, bonders, deleters)

    def get_atom1(self):
        return self._atom1.clone()

    def get_oxygen1(self):
        return self._oxygen1.clone()

    def get_hydrogen1(self):
        return self._hydrogen1.clone()

    def get_atom2(self):
        return self._atom2.clone()

    def get_oxygen2(self):
        return self._oxygen2.clone()

    def get_hydrogen2(self):
        return self._hydrogen2.clone()

    def get_atom3(self):
        self._atom3.clone()

    def get_atom4(self):
        return self._atom4.clone()


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
