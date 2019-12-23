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

    _functional_group = None
    _functional_group_smarts = None
    _bonder_smarts = None
    _deleter_smarts = None

    def __init_subclass__(cls):
        cls._functional_group_query = rdkit.MolFromSmarts(
            SMARTS=cls._functional_group_smarts,
        )
        cls._bonder_queries = [
            (rdkit.MolFromSmarts(smarts), count)
            for smarts, count in Counter(cls._bonder_smarts).items()
        ]
        cls._deleter_queries = [
            (rdkit.MolFromSmarts(smarts), count)
            for smarts, count in Counter(cls._deleter_smarts).items()
        ]

    def get_functional_groups(self, molecule):
        for ids in self._get_ids(molecule):
            atoms = tuple(molecule.get_atoms(ids.atom_ids))
            bonder_ids = set(ids.bonder_ids)
            deleter_ids = set(ids.deleter_ids)
            yield self._functional_group(
                atoms=atoms,
                bonders=tuple(
                    a for a in atoms if a.id in bonder_ids
                ),
                deleters=tuple(
                    a for a in atoms if a.id in deleter_ids
                ),
            )

    def _get_ids(self, molecule):
        rdkit_mol = molecule.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_mol)

        functional_groups = rdkit_mol.GetSubstructMatches(
            query=self._functional_group_query,
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

        for bonder, count in self._bonder_queries:
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

        for deleter, count in self._deleter_queries:
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
