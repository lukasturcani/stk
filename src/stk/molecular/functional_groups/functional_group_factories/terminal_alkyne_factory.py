from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import Alkyne


class TerminalAlkyneFactory(SmartsFunctionalGroupFactory):
    """
    Creates :class:`.Alkyne` instances.

    """

    _functional_group_smarts = '[*][C]#[C][H]'
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
            yield Alkyne(
                atom1=atoms[0],
                carbon1=atoms[1],
                carbon2=atoms[2],
                atom2=atoms[3],
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
