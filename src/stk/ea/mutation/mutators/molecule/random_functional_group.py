"""
Random Functional Group
=====================

"""

import numpy as np


from .mutator import MoleculeMutator
from ...records import MutationRecord
from ....molecule_records import MoleculeRecord

class RandomFunctionalGroup(MoleculeMutator):
    """
    Substitutes functional groups within building blocks.

    This mutator takes a :class:`ConstructedMolecule` and substitutes
    the existing building blocks with ones containing
    new functional groups from a given set. 
    """
    def __init__(
        self,
        functional_groups,
        is_replaceable_fg,
        is_replaceable_bb,
        replacement_count='all',
        name='RandomFunctionalGroup',
        random_seed=None,
    ):
        """
        Initialize a :class:`.RandomFunctionalGroup` instance.

        Parameters
        ----------
        functional_groups : :class:`tuple` of :class`GenericFunctionalGroup`
            A group of :class:`GenericFunctionalGroups` which are used 
            to replace the :class:`GenericFunctionalGroups` in the existing building blocks.

        is_replaceable : :class:`callable`
            A function which takes a :class:`.BuildingBlock` and
            returns ``True`` or ``False``. This function is applied to
            every building block in the molecule being mutated.
            Building blocks which returned ``True`` are liable for
            substitution by one of the molecules in `building_blocks`.

        name : :class:`str`, optional
            A name to help identify the mutator instance.

        random_seed : :class:`bool`, optional
            The random seed to use.

        """
