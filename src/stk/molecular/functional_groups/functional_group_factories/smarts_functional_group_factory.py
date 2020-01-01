import rdkit.Chem.AllChem as rdkit

from .functional_group_factory import FunctionalGroupFactory


class SmartsFunctionalGroupFactory(FunctionalGroupFactory):
    """

    """

    def __init__(self, functional_group_smarts, bonders, deleters):
        self._functional_group_query = rdkit.MolFromSmarts(
            SMARTS=functional_group_smarts,
        )
        self._functional_group_smarts = functional_group_smarts
        self._bonders = bonders
        self._deleters = deleters

    def _get_atom_ids(self, molecule):
        rdkit_mol = molecule.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_mol)

        yield from rdkit_mol.GetSubstructMatches(
            query=self._functional_group_query,
        )

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._functional_group_smarts}, {self._bonders}, '
            f'{self._deleters})'
        )

    def __str__(self):
        return repr(self)
