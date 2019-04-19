"""
Defines mutation operations via the :class:`Mutation` class.

.. _`adding mutation functions`:

Extending stk: Adding mutation functions.
-----------------------------------------

If a new mutation operation is to be added to ``stk`` it should be
added as a method in the :class:`Mutation` class defined in this
module. The only requirement is that the first argument is `macro_mol`,
excluding any `self` or `cls` arguments.

The naming requirement of `macro_mol` exists to help users identify
which arguments are handled automatically by ``stk`` and which they
need to define in the input file. The convention is that if the
mutation function takes an argument called  `macro_mol` it does not
have to be specified in the input file.

If the mutation function does not fit neatly into a single function,
make sure that any helper functions are private, i.e. that their names
start with a leading underscore.

"""

import logging
import numpy as np
import rdkit.Chem.AllChem as rdkit


logger = logging.getLogger(__name__)


class MutationError(Exception):
    """
    Used for errors which occuring during mutation operations.

    """

    ...


class Mutator:
    """
    Creates mutants.

    """

    def mutate(self, mol):
        """

        """

        raise NotImplementedError()


class RandomBuildingBlock(Mutator):
    """
    Substitutes random building blocks.

    This mutator takes a :class:`.MacroMolecule` and substitutes the
    building blocks with one chosen at random from a given set.

    Attributes
    ----------
    building_blocks : :class:`list` of :class:`.StructUnit`
        A group of molecules which are used to replace building blocks
        in molecules being mutated.

    key : :class:`function`
        A function which takes a :class:`.StructUnit` and returns
        ``True`` or ``False``. This function is applied to every
        building block of the molecule being mutated. Building blocks
        which returned ``True`` are liable for substition by one of the
        molecules in :attr:`building_blocks`.

    duplicate_building_blocks : :class:`bool`
        Indicates whether the building blocks used to construct the
        mutant must all be unique.

    Examples
    --------
    >>> # Create a molecule which is to be mutated.
    >>> bb1 = StructUnit2.smiles_init('NCCN', ['amine'])
    >>> bb2 = StructUnit2.smiles_init('O=CCC=O', ['aldehyde'])
    >>> polymer = Polymer([bb1, bb2], Linear('AB', [0, 0], n=3))

    >>> # Create molecules used to substitute building blocks.
    >>> mols = [
            StructUnit2.smiles_init('NC[Si]CCN', ['amine']),
            StructUnit2.smiles_init('NCCCCCCCN', ['amine']),
            StructUnit2.smiles_init('NC1CCCCC1N', ['amine'])
        ]
    >>> # Create the mutator.
    >>> random_bb = RandomBuildingBlock(
               mols=mols,
               key=lambda mol: mol.func_group_infos[0].name == 'amine'
        )
    >>> # Mutate a molecule.
    >>> mutant1 = random_bb.mutate(polymer)
    >>> # Mutate the molecule a second time.
    >>> mutant2 = random_bb.mutate(polymer)
    >>> # Mutate a mutant.
    >>> mutant3 = random_bb.mutate(mutant1)

    """

    def __init__(self,
                 building_blocks,
                 key,
                 duplicate_building_blocks=False):
        """
        Initializes a :class:`RandomBuildingBlock` instance.

        Parameters
        ----------
        building_blocks : :class:`list` of :class:`.StructUnit`
            A group of molecules which are used to replace building blocks
            in molecules being mutated.

        key : :class:`function`
            A function which takes a :class:`.StructUnit` and returns
            ``True`` or ``False``. This function is applied to every
            building block of the molecule being mutated. Building blocks
            which returned ``True`` are liable for substition by one of the
            molecules in :attr:`building_blocks`.

        duplicate_building_blocks : :class:`bool`
            Indicates whether the building blocks used to construct the
            mutant must all be unique.

        """

        self.building_blocks = building_blocks
        self.key = key
        self.duplicate_building_blocks = duplicate_building_blocks
        super().__init__()

    def mutate(self, mol):
        """
        Substitute a building block at random.

        Parameters
        ----------
        macro_mol : :class:`.MacroMolecule`
            The molecules which is to have its building block
            substituted.

        Returns
        -------
        :class:`.MacroMolecule`
            The produced mutant.

        """

        # Choose the building block which undergoes mutation.
        valid_bbs = [bb for bb in mol.building_blocks if self.key(bb)]
        chosen_bb = np.random.choice(valid_bbs)

        # If the mutant can have more than one of the same building
        # block, prevent only the building block getting replaced
        # from being used as a replacement.
        if self.duplicate_building_blocks:
            excluded_bbs = {chosen_bb}
        # If the mutant is to be composed of unique building blocks
        # only, prevent any building block already present in the
        # macro_mol from being used as a replacement.
        else:
            excluded_bbs = set(mol.building_blocks)
        # Make sure that the building block itself will not be picked.
        mols = [
            mol for mol in self.building_blocks
            if mol not in excluded_bbs
        ]

        # Choose a replacement building block.
        replacement = np.random.choice(mols)

        # Build the new MacroMolecule.
        new_bbs = [
            bb for bb in mol.building_blocks if bb is not chosen_bb
        ]
        new_bbs.append(replacement)
        return mol.__class__(new_bbs, mol.topology)


class SimilarBuildingBlock(Mutator):
    """
    Substitutes similar building blocks.

    This mutator takes a :class:`.MacroMolecule` and substitutes the
    building blocks with the most similar one from a given set.

    Attributes
    ----------
    building_blocks : :class:`list` of :class:`.StructUnit`
        A group of molecules which are used to replace building blocks
        in molecules being mutated.

    key : :class:`function`
        A function which takes a :class:`.StructUnit` and returns
        ``True`` or ``False``. This function is applied to every
        building block of the molecule being mutated. Building blocks
        which returned ``True`` are liable for substition by one of the
        molecules in :attr:`building_blocks`.

    duplicate_building_blocks : :class:`bool`
        Indicates whether the building blocks used to construct the
        mutant must all be unique.

    Examples
    --------
    >>> # Create a molecule which is to be mutated.
    >>> bb1 = StructUnit2.smiles_init('NCCN', ['amine'])
    >>> bb2 = StructUnit2.smiles_init('O=CCC=O', ['aldehyde'])
    >>> polymer = Polymer([bb1, bb2], Linear('AB', [0, 0], n=3))

    >>> # Create molecules used to substitute building blocks.
    >>> mols = [
            StructUnit2.smiles_init('NC[Si]CCN', ['amine']),
            StructUnit2.smiles_init('NCCCCCCCN', ['amine']),
            StructUnit2.smiles_init('NC1CCCCC1N', ['amine'])
        ]
    >>> # Create the mutator.
    >>> random_bb = RandomBuildingBlock(
               mols=mols,
               key=lambda mol: mol.func_group_infos[0].name == 'amine'
        )
    >>> # Mutate a molecule.
    >>> mutant1 = random_bb.mutate(polymer)
    >>> # Mutate the molecule a second time.
    >>> mutant2 = random_bb.mutate(polymer)
    >>> # Mutate a mutant.
    >>> mutant3 = random_bb.mutate(mutant1)

    """

    def __init__(self,
                 building_blocks,
                 key,
                 duplicate_building_blocks):
        """
        Initializes a :class:`SimilarBuildingBlock` instance.

        Parameters
        ----------
        building_blocks : :class:`list` of :class:`.StructUnit`
            A group of molecules which are used to replace building blocks
            in molecules being mutated.

        key : :class:`function`
            A function which takes a :class:`.StructUnit` and returns
            ``True`` or ``False``. This function is applied to every
            building block of the molecule being mutated. Building blocks
            which returned ``True`` are liable for substition by one of the
            molecules in :attr:`building_blocks`.

        duplicate_building_blocks : :class:`bool`
            Indicates whether the building blocks used to construct the
            mutant must all be unique.

        """

        self.building_blocks = building_blocks
        self.key = key
        self.duplicate_building_blocks = duplicate_building_blocks
        super().__init__()

    def mutate(self, mol):
        """
        Substitute a similar building block.

        Parameters
        ----------
        macro_mol : :class:`.MacroMolecule`
            The molecules which is to have its building block
            substituted.

        Returns
        -------
        :class:`.MacroMolecule`
            The produced mutant.

        """

        if not hasattr(self, '_similar_bb'):
            # This will map every macro_mol to a set of iterators which
            # yield the most similar building blocks. 1 iterator for
            # each previously chosen building block.
            self._similar_bb = {}
        if mol not in self._similar_bb:
            # Maps the macro_mol to a dict. The dict maps each
            # building block of the macro mol to an iterator.
            # The iterators yield the next most similar molecules in
            # `mols` to the building block.
            self._similar_bb[mol] = {}

        # Choose the building block which undergoes mutation.
        valid_bbs = [bb for bb in mol.building_blocks if self.key(bb)]
        chosen_bb = np.random.choice(valid_bbs)

        # Create a mapping from rdkit molecules to the StructUnits.
        mol_map = {
            struct_unit.inchi: struct_unit
            for struct_unit in self.building_blocks
        }

        # If the building block has not been chosen before, create an
        # iterator yielding similar molecules from `mols` for it.
        if chosen_bb not in self._similar_bb[mol]:
            rdkit_mols = (m.mol for m in self.building_blocks)
            self._similar_bb[mol][chosen_bb] = iter(
                      chosen_bb.similar_molecules(rdkit_mols))

        try:
            sim_mol = next(self._similar_bb[mol][chosen_bb])[-1]
        except StopIteration:
            rdkit_mols = (m.mol for m in self.building_blocks)
            self._similar_bb[mol][chosen_bb] = iter(
                      chosen_bb.similar_molecules(rdkit_mols))
            sim_mol = next(self._similar_bb[mol][chosen_bb])[-1]

        sim_mol_inchi = rdkit.MolToInchi(sim_mol)
        sim_struct_unit = mol_map[sim_mol_inchi]

        # If the most similar molecule in `mols` is itself, then take
        # the next most similar one.
        if sim_struct_unit is chosen_bb:
            sim_mol = next(self._similar_bb[mol][chosen_bb])[-1]
            sim_mol_inchi = rdkit.MolToInchi(sim_mol)
            sim_struct_unit = mol_map[sim_mol_inchi]

        # Build the new MacroMolecule.
        new_bbs = [bb for bb in mol.building_blocks if
                   bb is not chosen_bb]
        new_bbs.append(sim_struct_unit)
        return mol.__class__(new_bbs, mol.topology)


class RandomTopology(Mutator):
    """

        Changes `macro_mol` topology to a random one from `topologies`.

        A new instance of the same type as `macro_mol` is created. I.e.
        if `macro_mol` was a :class:`.Polymer` instance then a
        :class:`.Polymer` instance will be returned.

        Parameters
        ----------
        macro_mol : :class:`.MacroMolecule`
            The macromolecule which is to be mutated.

        topologies : :class:`list` of :class:`.Topology`
            This lists holds the topology instances from which one is
            selected at random to form a new molecule.

    """

    def __init__(self, topologies):
        """

        Changes `macro_mol` topology to a random one from `topologies`.

        A new instance of the same type as `macro_mol` is created. I.e.
        if `macro_mol` was a :class:`.Polymer` instance then a
        :class:`.Polymer` instance will be returned.

        Parameters
        ----------
        macro_mol : :class:`.MacroMolecule`
            The macromolecule which is to be mutated.

        topologies : :class:`list` of :class:`.Topology`
            This lists holds the topology instances from which one is
            selected at random to form a new molecule.

        """

        self.topologies = topologies
        super().__init__()

    def mutate(self, mol):
        """


        Returns
        -------
        :class:`.MacroMolecule`
            A molecule generated by initializing a new instance
            with all the same parameters as `macro_mol` except for the
            topology.

        """

        tops = [x for x in topologies if
                repr(x) != repr(macro_mol.topology)]
        topology = np.random.choice(tops)
        return macro_mol.__class__(macro_mol.building_blocks, topology)
