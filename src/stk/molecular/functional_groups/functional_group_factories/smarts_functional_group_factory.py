import rdkit.Chem.AllChem as rdkit
from collections import Counter

from stk import FunctionalGroupFactory, flatten


class _FunctionalGroupIds:
    def __init__(self, atom_ids, bonder_ids, deleter_ids):
        self.atom_ids = atom_ids
        self.bonder_ids = bonder_ids
        self.deleter_ids = deleter_ids


class SmartsFunctionalGroupFactory(FunctionalGroupFactory):
    """

    """

    functional_group_smarts = None
    bonder_smarts = None
    deleter_smarts = None

    def __init__(self):
        self._functional_group = rdkit.MolFromSmarts(
            SMARTS=self.functional_group_smarts,
        )
        self._bonders = [
            (rdkit.MolFromSmarts(smarts), count)
            for smarts, count in Counter(self.bonder_smarts).items()
        ]
        self._deleters = [
            (rdkit.MolFromSmarts(smarts), count)
            for smarts, count in Counter(self.deleter_smarts).items()
        ]

    def _get_ids(self, molecule):
        rdkit_mol = molecule.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_mol)

        functional_groups = rdkit_mol.GetSubstructMatches(
            query=self._functional_group,
        )

        ids = zip(
            functional_groups,
            self._get_bonders(rdkit_mol, functional_groups),
            self._get_deleters(rdkit_mol, functional_groups),
        )
        for atom_ids, bonder_ids, deleter_ids in ids:
            yield _FunctionalGroupIds(
                atom_ids=atom_ids,
                bonder_ids=bonder_ids,
                deleter_ids=deleter_ids,
            )

    def _get_bonders(self, molecule, functional_groups):
        """

        """

        for bonder, count in self._bonders:
            matches = set(flatten(
                molecule.GetSubstructMatches(bonder)
            ))

            matched_bonders = [
                [aid for aid in fg if aid in matches]
                for fg in functional_groups
            ]

            for fg_id, fg in enumerate(functional_groups):
                yield matched_bonders[fg_id][:count]

    def _get_deleters(self, molecule, functional_groups):
        """

        """

        for deleter, count in self._deleters:
            matches = set(flatten(
                molecule.GetSubstructMatches(deleter)
            ))

            matched_deleters = [
                [aid for aid in fg if aid in matches]
                for fg in functional_groups
            ]

            for fg_id, fg in enumerate(functional_groups):
                yield matched_deleters[fg_id][:count]

    def __repr__(self):
        return f'{self.__class__.__name__}()'

    def __str__(self):
        return repr(self)
