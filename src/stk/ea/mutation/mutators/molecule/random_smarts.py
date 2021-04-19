"""
Random Smarts
=====================

"""

import numpy as np


from .mutator import MoleculeMutator
from ...records import MutationRecord
from ....molecule_records import MoleculeRecord
from rdkit import Chem
from rdkit.Chem.rdchem import BondType
from .....molecular import BuildingBlock
from rdkit import Chem
import stk


class RandomSmarts(MoleculeMutator):
    """
    Substitutes functional groups within building blocks.

    This mutator takes a :class:`ConstructedMolecule` and substitutes
    the existing building blocks with ones containing
    new functionalities, as specified in the SMARTS strings.
    Atoms in functional groups cannot be adjusted as these are used to
    construct the :class: `ConstructedMolecule`, and will
    result in an error being raised.
    Only external functionalities can be replaced.
    """
    def __init__(
        self,
        query_smarts,
        replacement_smarts,
        is_replaceable,
        replacement_functional_groups,
        replacement_specifier='all',
        name='RandomSmarts',
        random_seed=None,
    ):
        """
        Initialize a :class:`.RandomSmarts` instance.

        Parameters
        ----------
        query_smarts : :class`string`
            SMARTS string to match on the :class: `BuildingBlock`

        replacement_smarts : :class`string`
            SMARTS string to replace those in `query_smarts` on the
            :class: `BuildingBlock`

        is_replaceable : :class:`callable`
            A function which takes a :class:`.BuildingBlock` and
            returns ``True`` or ``False``. This function is applied to
            every building block in the molecule being mutated.
            Building blocks which returned ``True`` are liable for
            substitution by one of the molecules in `building_blocks`.

        replacement_specifier : :class:`str`, optional
            Specifies the number of SMARTS replacements to make.
            If `all`, all matching SMARTS Swill be replacemed.
            If `one`, a single random SMARTS replacement will occur.

        replacement_functional_groups : :class:`list` of :class:`.FunctionalGroupFactory`
            The functional group factories used to assign functional groups within the generated molecule.

        name : :class:`str`, optional
            A name to help identify the mutator instance.

        random_seed : :class:`bool`, optional
            The random seed to use.

        """
        self._query_smarts = query_smarts
        self._replacement_smarts = replacement_smarts
        self._is_replaceable = is_replaceable
        self._name = name
        self._replacement_specifier = replacement_specifier
        self._generator = np.random.RandomState(random_seed)
        self._replacement_functional_groups = replacement_functional_groups

    def mutate(self, record):
        replaceable_building_blocks = tuple(filter(
            self._is_replaceable,
            record.get_molecule().get_building_blocks(),
        ))
        replaced_building_block = self._generator.choice(
            a=replaceable_building_blocks,
        )
        rdmol = replaced_building_block.to_rdkit_mol()
        query = Chem.MolFromSmarts(self._query_smarts)
        replacer_smarts = Chem.MolFromSmarts(self._replacement_smarts)
        if self._replacement_specifier == 'all':
            new_rdmol = Chem.rdmolops.ReplaceSubstructs(rdmol, query, replacer_smarts, replaceAll=True)[0]

        elif self._replacement_specifier == 'one':
            new_rdmols = Chem.rdmolops.ReplaceSubstructs(rdmol, query, replacer_smarts, replaceAll=True)
            new_rdmol = self._generator.choice(
                new_rdmols
            )
        else:
            raise RuntimeError('Invalid argument for replacement_specifier received')
        # Create new BuildingBlock
        replacement = BuildingBlock.init_from_rdkit_mol(new_rdmol, functional_groups=self._replacement_functional_groups)
        graph = record.get_topology_graph().with_building_blocks({
            replaced_building_block: replacement,
        })

        return MutationRecord(
            molecule_record=MoleculeRecord(graph),
            mutator_name=self._name,
        )
