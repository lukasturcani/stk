import rdkit.Chem.AllChem as rdkit
from collections import Counter

from .functional_group_factory import FunctionalGroupFactory
from stk import flatten


class _FunctionalGroupIds:
    def __init__(self, atom_ids, bonder_ids, deleter_ids):
        self.atom_ids = atom_ids
        self.bonder_ids = bonder_ids
        self.deleter_ids = deleter_ids


class SmartsFunctionalGroupFactory(FunctionalGroupFactory):
    """

    """

    _functional_group_smarts = None
    _bonder_smarts = None
    _deleter_smarts = None

    def __init_subclass__(cls):
        cls._set_queries(cls)

    @staticmethod
    def _set_queries(obj):
        obj._functional_group_query = rdkit.MolFromSmarts(
            SMARTS=obj._functional_group_smarts,
        )
        obj._bonder_queries = [
            (rdkit.MolFromSmarts(smarts), count)
            for smarts, count in Counter(obj._bonder_smarts).items()
        ]
        obj._deleter_queries = [
            (rdkit.MolFromSmarts(smarts), count)
            for smarts, count in Counter(obj._deleter_smarts).items()
        ]

    def _get_ids(self, molecule):
        rdkit_mol = molecule.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_mol)

        functional_groups = rdkit_mol.GetSubstructMatches(
            query=self._functional_group_query,
        )

        ids = zip(
            functional_groups,
            self._get_bonder_ids(rdkit_mol, functional_groups),
            self._get_deleter_ids(rdkit_mol, functional_groups),
        )
        for atom_ids, bonder_ids, deleter_ids in ids:
            yield _FunctionalGroupIds(
                atom_ids=atom_ids,
                bonder_ids=bonder_ids,
                deleter_ids=deleter_ids,
            )

    def _get_bonder_ids(self, molecule, functional_groups):
        """

        """

        bonders = [[] for _ in range(len(functional_groups))]

        for bonder, count in self._bonder_queries:
            matches = set(flatten(
                molecule.GetSubstructMatches(bonder)
            ))

            matched_bonders = [
                [aid for aid in fg if aid in matches]
                for fg in functional_groups
            ]
            for fg_id, bonder_ids in enumerate(matched_bonders):
                bonders[fg_id].extend(bonder_ids[:count])

        yield from bonders

    def _get_deleter_ids(self, molecule, functional_groups):
        """

        """

        deleters = [[] for _ in range(len(functional_groups))]

        for deleter, count in self._deleter_queries:
            matches = set(flatten(
                molecule.GetSubstructMatches(deleter)
            ))

            matched_deleters = [
                [aid for aid in fg if aid in matches]
                for fg in functional_groups
            ]

            for fg_id, deleter_ids in enumerate(matched_deleters):
                deleters[fg_id].extend(deleter_ids[:count])

        yield from deleters

    def __repr__(self):
        return f'{self.__class__.__name__}()'

    def __str__(self):
        return repr(self)
