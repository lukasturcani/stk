"""
Random SMARTS
=============

"""

import numpy as np


from .mutator import MoleculeMutator
from ...records import MutationRecord
from ....molecule_records import MoleculeRecord
import rdkit.Chem.AllChem as rdkit
from .....molecular import BuildingBlock


class RandomSmarts(MoleculeMutator):
    """
    Substitutes functional groups within building blocks.

    This mutator takes a :class:`.ConstructedMolecule` and substitutes
    the existing building blocks with ones containing
    new functionalities, as specified in the SMILES strings.
    Atoms in functional groups used to
    construct the :class: `.ConstructedMolecule` cannot be changed,
    and will result in an error being raised, unless the
    new functional groups are specified in `replacement_fgs`.

    Examples
    --------
    *Smarts Mutation*

    .. testcode:: constructed-molecule-mutation

        import stk

        # Create a molecule which is to be mutated.
        bb1 = stk.BuildingBlock(
            smiles='NCC(CF)CCN',
            functional_groups=[stk.PrimaryAminoFactory()]
        )
        bb2 = stk.BuildingBlock('O=CCCCC=O', [stk.AldehydeFactory()])
        polymer = stk.MoleculeRecord(
            topology_graph=stk.polymer.Linear((bb1, bb2), 'AB', 3),
        )

        # Create the mutator.

        def has_primary_amino_group(building_block):
            fg, = building_block.get_functional_groups(0)
            return type(fg) is stk.PrimaryAmino

        random_smarts = stk.RandomSmarts(
            query_smarts="F",
            replacement_smiles="Br",
            is_replaceable=has_primary_amino_group,
            replacement_fgs=[stk.PrimaryAminoFactory()],
        )

        # Mutate a molecule.
        mutation_record1 = random_smarts.mutate(polymer)

        # Create a mutator that will fail.

        random_smarts_failed = stk.RandomSmarts(
            query_smarts="N",
            replacement_smiles="Br",
            is_replaceable=has_primary_amino_group,
            replacement_fgs=[stk.PrimaryAminoFactory()],
        )
        # This will raise an exception`.
        mutation_record2 = random_smarts.mutate(polymer)

    """

    def __init__(
        self,
        query_smarts,
        replacement_smiles,
        is_replaceable,
        replacement_functional_groups,
        name='RandomSmarts',
        random_seed=None,
    ):
        """
        Initialize a :class:`.RandomSmarts` instance.

        Parameters
        ----------
        query_smarts : :class:`str`
            SMARTS string to match on the :class:`.BuildingBlock`

        replacement_smiles : :class`str`
            SMILES string to replace those in `query_smarts` on the
            :class: `.BuildingBlock`

        is_replaceable : :class:`callable`
            A function which takes a :class:`.BuildingBlock` and
            returns ``True`` or ``False``. This function is applied to
            every building block in the molecule being mutated.
            Building blocks which returned ``True`` are liable for
            substitution.

        replacement_functional_groups : :class:`iterable` of
        :class:`.FunctionalGroupFactory`
            The functional group factories functional groups to the
            mutated building block.

        name : :class:`str`, optional
            A name to help identify the mutator instance.

        random_seed : :class:`bool`, optional
            The random seed to use.

        """

        self._query_smarts = query_smarts
        self._replacement_smiles = replacement_smiles
        self._is_replaceable = is_replaceable
        self._name = name
        self._generator = np.random.RandomState(random_seed)
        self._replacement_functional_groups = tuple(
            replacement_functional_groups,
        )

    def mutate(self, record):
        replaceable_building_blocks = tuple(
            filter(
                self._is_replaceable,
                record.get_molecule().get_building_blocks(),
            ),
        )
        replaced_building_block = self._generator.choice(
            a=replaceable_building_blocks,
        )
        rdmol = replaced_building_block.to_rdkit_mol()
        query = rdkit.MolFromSmarts(self._query_smarts)
        replacer_smarts = rdkit.MolFromSmarts(self._replacement_smiles)
        new_rdmol = rdkit.rdmolops.ReplaceSubstructs(
            mol=rdmol,
            query=query,
            replacement=replacer_smarts,
            replaceAll=True,
        )[0]
        # Create new BuildingBlock
        replacement = BuildingBlock.init_from_rdkit_mol(
            molecule=new_rdmol,
            functional_groups=self._replacement_functional_groups,
        )
        graph = record.get_topology_graph().with_building_blocks({
            replaced_building_block: replacement,
        })

        return MutationRecord(
            molecule_record=MoleculeRecord(graph),
            mutator_name=self._name,
        )
