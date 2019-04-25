"""
Defines mutators.

Mutators are objects which mutate molecules. They inherit
:class:`Mutator` and define a method :meth:`~Mutator.mutate`. This
method must take a single molecule and return a mutant.

Examples of how mutators work can be seen the documentation of
the various :class:`Mutator` classes, for example
:class:`RandomBuildingBlock`, :class:`SimilarBuildingBlock`
or :class:`RandomMutation`.

.. _`adding mutators`:

Extending stk: Making new mutators.
-----------------------------------

Mutators must simple inherit the :class:`Mutator` class and define a
method called :meth:`~Mutatator.mutate`, which take a single molecule
and returns a mutant molecule. There are no requirements besides this.

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
        Mutates a molecule.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be mutated.

        Returns
        -------
        mol : :class:`.Molecule`
            The mutant molecule.

        """

        raise NotImplementedError()


class RandomMutation(Mutator):
    """
    Uses a random :class:`Mutator` to carry out mutations.

    Attributes
    ----------
    mutators : :class:`tuple` of :class:`Mutator`
        :class:`Mutator` objects which are used to carry out the
        mutations.

    weights : :class:`list` of :class:`float`
        The probability that each :class:`Mutator` will be selected
        to carry out a mutation.

    Examples
    --------
    .. code-block:: python

        # Create a molecule which is to be mutated.
        bb1 = StructUnit2.smiles_init('NCCN', ['amine'])
        bb2 = StructUnit3.smiles_init('O=CCC(=O)CC=O', ['aldehyde'])
        cage = Cage([bb1, bb2], FourPlusSix())


        # Create the first mutator.
        topologies = [
            TwoPlusThree(),
            EightPlusTwelve(),
            Dodecahedron()
        ]
        random_topology = RandomTopology(topologies)

        # Create the second and third mutator.
        building_blocks = [
            StructUnit2.smiles_init('NC[Si]CCN', ['amine']),
            StructUnit2.smiles_init('NCCCCCCCN', ['amine']),
            StructUnit2.smiles_init('NC1CCCCC1N', ['amine'])
        ]

        random_bb = RandomBuildingBlock(
            building_blocks=building_blocks,
            key=lambda mol: mol.func_group_infos[0].name == 'amine'
        )

        similar_bb = SimilarBuildingBlock(
            building_blocks=building_blocks,
            key=lambda mol: mol.func_group_infos[0].name == 'amine'
        )

        # Create the mutator used to carry out the mutations.
        random_mutator = RandomMutation(random_topology,
                                        random_bb,
                                        similar_bb)

        # Mutate a molecule, one of random_topology,
        # random_bb or similar_bb will be used.
        mutant1 = random_mutator.mutate(cage)

        # Mutate the molecule a second time, one of random_topology,
        # random_bb or similar_bb will be used.
        mutant2 = random_mutator.mutate(cage)

        # Mutate a mutant, one of random_topology,
        # random_bb or similar_bb will be used.
        mutant3 = random_mutator.mutate(mutant1)

    """

    def __init__(self, *mutators, weights=None):
        """
        Initializes a :class:`RandomMutation` instance.

        Parameters
        ----------
        *mutators : :class:`Mutator`
            :class:`Mutator` objects which are used to carry out the
            mutations.

        weights : :class:`list` of :class:`float`, optional
            The probability that each :class:`Mutator` will be selected
            to carry out a mutation.

        """

        self.mutators = mutators
        self.weights = weights

    def mutate(self, mol):
        """
        Create a mutant with a :class:`Mutator` in :attr:`mutators`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule which is to be mutated.

        Returns
        -------
        :class:`.Molecule`
            The produced mutant.

        """

        mutator = np.random.choice(self.mutators, p=self.weights)
        return mutator.mutate(mol)


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
    .. code-block:: python

        # Create a molecule which is to be mutated.
        bb1 = StructUnit2.smiles_init('NCCN', ['amine'])
        bb2 = StructUnit2.smiles_init('O=CCC=O', ['aldehyde'])
        polymer = Polymer([bb1, bb2], Linear('AB', [0, 0], n=3))

        # Create molecules used to substitute building blocks.
        building_blocks = [
            StructUnit2.smiles_init('NC[Si]CCN', ['amine']),
            StructUnit2.smiles_init('NCCCCCCCN', ['amine']),
            StructUnit2.smiles_init('NC1CCCCC1N', ['amine'])
        ]

        # Create the mutator.
        random_bb = RandomBuildingBlock(
            building_blocks=building_blocks,
            key=lambda mol: mol.func_group_infos[0].name == 'amine'
        )

        # Mutate a molecule.
        mutant1 = random_bb.mutate(polymer)

        # Mutate the molecule a second time.
        mutant2 = random_bb.mutate(polymer)

        # Mutate a mutant.
        mutant3 = random_bb.mutate(mutant1)

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
            A group of molecules which are used to replace building
            blocks in molecules being mutated.

        key : :class:`function`
            A function which takes a :class:`.StructUnit` and returns
            ``True`` or ``False``. This function is applied to every
            building block of the molecule being mutated. Building
            blocks which returned ``True`` are liable for substition by
            one of the molecules in :attr:`building_blocks`.

        duplicate_building_blocks : :class:`bool`, optional
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
            The molecule which is to have its building block
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

    _similar_bbs : :class:`dict`
        Maps a :class:`.MacroMolecule` to multiple :class:`iterator`,
        each of which yields the most similar molecules in
        :attr:`building_blocks`, in order.

    Examples
    --------
    .. code-block:: python

        # Create a molecule which is to be mutated.
        bb1 = StructUnit2.smiles_init('NCCN', ['amine'])
        bb2 = StructUnit2.smiles_init('O=CCC=O', ['aldehyde'])
        polymer = Polymer([bb1, bb2], Linear('AB', [0, 0], n=3))

        # Create molecules used to substitute building blocks.
        building_blocks = [
            StructUnit2.smiles_init('NC[Si]CCN', ['amine']),
            StructUnit2.smiles_init('NCCCCCCCN', ['amine']),
            StructUnit2.smiles_init('NC1CCCCC1N', ['amine'])
        ]

        # Create the mutator.
        similar_bb = SimilarBuildingBlock(
            building_blocks=building_blocks,
            key=lambda mol: mol.func_group_infos[0].name == 'amine'
        )

        # Mutate a molecule.
        mutant1 = random_bb.mutate(polymer)

        # Mutate the molecule a second time.
        mutant2 = random_bb.mutate(polymer)

        # Mutate a mutant.
        mutant3 = random_bb.mutate(mutant1)

    """

    def __init__(self,
                 building_blocks,
                 key,
                 duplicate_building_blocks):
        """
        Initializes a :class:`RandomBuildingBlock` instance.

        Parameters
        ----------
        building_blocks : :class:`list` of :class:`.StructUnit`
            A group of molecules which are used to replace building
            blocks in molecules being mutated.

        key : :class:`function`
            A function which takes a :class:`.StructUnit` and returns
            ``True`` or ``False``. This function is applied to every
            building block of the molecule being mutated. Building
            blocks which returned ``True`` are liable for substition by
            one of the molecules in :attr:`building_blocks`.

        duplicate_building_blocks : :class:`bool`, optional
            Indicates whether the building blocks used to construct the
            mutant must all be unique.

        """

        self.building_blocks = building_blocks
        self.key = key
        self.duplicate_building_blocks = duplicate_building_blocks
        self._similar_bbs = {}
        super().__init__()

    def mutate(self, mol):
        """
        Substitute a similar building block.

        Parameters
        ----------
        macro_mol : :class:`.MacroMolecule`
            The molecule which is to have its building block
            substituted.

        Returns
        -------
        :class:`.MacroMolecule`
            The produced mutant.

        """

        if mol not in self._similar_bbs:
            # Maps the macro_mol to a dict. The dict maps each
            # building block of the macro mol to an iterator.
            # The iterators yield the next most similar molecules in
            # `mols` to the building block.
            self._similar_bbs[mol] = {}

        # Choose the building block which undergoes mutation.
        valid_bbs = [bb for bb in mol.building_blocks if self.key(bb)]
        chosen_bb = np.random.choice(valid_bbs)

        # Create a mapping from inchis to the StructUnits.
        mol_map = {
            struct_unit.inchi: struct_unit
            for struct_unit in self.building_blocks
        }

        # If the building block has not been chosen before, create an
        # iterator yielding similar molecules from `mols` for it.
        if chosen_bb not in self._similar_bbs[mol]:
            rdkit_mols = (m.mol for m in self.building_blocks)
            self._similar_bbs[mol][chosen_bb] = iter(
                      chosen_bb.similar_molecules(rdkit_mols)
            )

        try:
            _, sim_mol = next(self._similar_bbs[mol][chosen_bb])
        except StopIteration:
            rdkit_mols = (m.mol for m in self.building_blocks)
            self._similar_bbs[mol][chosen_bb] = iter(
                chosen_bb.similar_molecules(rdkit_mols)
            )
            _, sim_mol = next(self._similar_bbs[mol][chosen_bb])

        sim_mol_inchi = rdkit.MolToInchi(sim_mol)
        sim_struct_unit = mol_map[sim_mol_inchi]

        # If the most similar molecule in `mols` is itself, then take
        # the next most similar one.
        if sim_struct_unit is chosen_bb:
            _, sim_mol = next(self._similar_bbs[mol][chosen_bb])
            sim_mol_inchi = rdkit.MolToInchi(sim_mol)
            sim_struct_unit = mol_map[sim_mol_inchi]

        # Build the new MacroMolecule.
        new_bbs = [
            bb for bb in mol.building_blocks if bb is not chosen_bb
        ]
        new_bbs.append(sim_struct_unit)
        return mol.__class__(new_bbs, mol.topology)


class RandomTopology(Mutator):
    """
    Changes topologies at random.

    Parameters
    ----------
    topologies : :class:`list` of :class:`.Topology`
        This :class:`list` holds the topology instances from which one
        is selected at random to form a mutant.

    Examples
    --------
    .. code-block:: python

        # Create a molecule which is to be mutated.
        bb1 = StructUnit2.smiles_init('NCCN', ['amine'])
        bb2 = StructUnit3.smiles_init('O=CCC(=O)CC=O', ['aldehyde'])
        cage = Cage([bb1, bb2], FourPlusSix())

        # Create topologies used for substition.
        topologies = [
            TwoPlusThree(),
            EightPlusTwelve(),
            Dodecahedron()
        ]

        # Create the mutator.
        random_topology = RandomTopology(topologies)

        # Mutate a molecule.
        mutant1 = random_topology.mutate(cage)

        # Mutate the molecule a second time.
        mutant2 = random_topology.mutate(cage)

        # Mutate a mutant.
        mutant3 = random_topology.mutate(mutant1)

    """

    def __init__(self, topologies):
        """
        Initializes a :class:`RandomTopology` instance.

        Parameters
        ----------
        topologies : :class:`list` of :class:`.Topology`
            This lists holds the topology instances from which one is
            selected at random to form a new molecule.

        """

        self.topologies = topologies
        super().__init__()

    def mutate(self, mol):
        """
        Substitute a random topology.

        Parameters
        ----------
        mol : :class:`.MacroMolecule`
            The molecule which is to have its topology substituted.

        Returns
        -------
        :class:`.MacroMolecule`
            The produced mutant.

        """

        tops = [
            x for x in self.topologies if repr(x) != repr(mol.topology)
        ]
        topology = np.random.choice(tops)
        return mol.__class__(mol.building_blocks, topology)
