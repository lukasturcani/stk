from .smarts_functional_group_factory import (
    SmartsFunctionalGroupFactory,
)
from ..functional_groups import PrimaryAmine


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
