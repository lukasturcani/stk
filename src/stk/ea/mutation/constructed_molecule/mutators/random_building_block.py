from .mutator import Mutator


class RandomBuildingBlock(Mutator):
    """
    Substitutes random building blocks.

    This mutator takes a :class:`.ConstructedMolecule` and substitutes
    the building blocks with one chosen at random from a given set.

    Examples
    --------
    .. code-block:: python

        import stk

        # Create a molecule which is to be mutated.
        bb1 = stk.BuildingBlock('NCCN', ['amine'])
        bb2 = stk.BuildingBlock('O=CCC=O', ['aldehyde'])
        polymer = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2],
            topology_graphs=stk.polymer.Linear('AB', [0, 0], n=3)
        )

        # Create molecules used to substitute building blocks.
        building_blocks = [
            stk.BuildingBlock('NC[Si]CCN', ['amine']),
            stk.BuildingBlock('NCCCCCCCN', ['amine']),
            stk.BuildingBlock('NC1CCCCC1N', ['amine'])
        ]

        # Create the mutator.
        random_bb = stk.RandomBuildingBlock(
            building_blocks=building_blocks,
            key=lambda mol: mol.func_groups[0].fg_type.name == 'amine'
        )

        # Mutate a molecule.
        mutant1 = random_bb.mutate(polymer)

        # Mutate the molecule a second time.
        mutant2 = random_bb.mutate(polymer)

        # Mutate a mutant.
        mutant3 = random_bb.mutate(mutant1)

    """

    def __init__(
        self,
        building_blocks,
        key,
        duplicate_building_blocks=False,
        random_seed=None,
        use_cache=False
    ):
        """
        Initialize a :class:`RandomBuildingBlock` instance.

        Parameters
        ----------
        building_blocks : :class:`tuple` of :class:`.Molecule`
            A group of molecules which are used to replace building
            blocks in molecules being mutated.

        key : :class:`function`
            A function which takes a :class:`.Molecule` and
            returns ``True`` or ``False``. This function is applied to
            every building block in the molecule being mutated.
            Building blocks which returned ``True`` are liable for
            substition by one of the molecules in `building_blocks`.

        duplicate_building_blocks : :class:`bool`, optional
            Indicates whether the building blocks used to construct the
            mutant must all be unique.

        random_seed : :class:`bool`, optional
            The random seed to use.

        use_cache : :class:`bool`, optional
            Toggles use of the molecular cache.

        """

        self._building_blocks = building_blocks
        self._key = key
        self._duplicate_building_blocks = duplicate_building_blocks
        self._generator = np.random.RandomState(random_seed)
        super().__init__(use_cache=use_cache)

    def _mutate(self, mol):
        """
        Return a mutant of `mol`.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The molecule to be mutated.

        Returns
        -------
        mol : :class:`.ConstructedMolecule`
            The mutant.

        """

        # Choose the building block which undergoes mutation.
        valid_bbs = [
            bb for bb in mol.building_block_vertices if self._key(bb)
        ]
        chosen_bb = self._generator.choice(valid_bbs)

        # If the mutant can have more than one of the same building
        # block, prevent only the building block getting replaced
        # from being used as a replacement.
        if self._duplicate_building_blocks:
            excluded_bbs = {chosen_bb}
        # If the mutant is to be composed of unique building blocks
        # only, prevent any building block already present in the
        # mol from being used as a replacement.
        else:
            excluded_bbs = set(mol.building_block_vertices.keys())
        # Make sure that the excluded_bbs will not be picked.
        mols = [
            mol for mol in self._building_blocks
            if mol not in excluded_bbs
        ]

        # Choose a replacement building block.
        replacement = self._generator.choice(mols)

        # Build the new ConstructedMolecule.
        new_bbs = [
            bb for bb in mol.building_block_vertices
            if bb is not chosen_bb
        ]
        new_bbs.append(replacement)
        return mol.__class__(
            building_blocks=new_bbs,
            topology_graph=mol.topology_graph,
            use_cache=self._use_cache
        )
