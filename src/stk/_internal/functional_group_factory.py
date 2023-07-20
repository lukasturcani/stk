from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass

import rdkit.Chem.AllChem as rdkit

from stk._internal.functional_group import FunctionalGroup


@dataclass(frozen=True, slots=True)
class FunctionalGroupFactory:
    """
    Creates functional groups.

    The purpose of a functional group factory is to create
    :class:`.FunctionalGroup` instances. It allows the user to avoid
    creating :class:`.FunctionalGroup` instances manually.
    """

    smarts: str
    bonders: Sequence[int]
    deleters: Sequence[int]

    def get_functional_groups(
        self,
        molecule: rdkit.Mol,
    ) -> Iterator[FunctionalGroup]:
        """
        Yield functional groups in `molecule`.

        .. note::

            This function expects `molecule` to be sanitized.

        Parameters:
            molecule:
                The molecule, whose functional groups are to be found.
        Yields:
            A functional group in `molecule`.
        """
        for atom_ids in molecule.GetSubstructMatches(
            query=rdkit.MolFromSmarts(self.smarts),
        ):
            yield FunctionalGroup(
                bonders=[atom_ids[i] for i in self.bonders],
                deleters=[atom_ids[i] for i in self.deleters],
            )


def bromo(
    bonders: Iterable[int] = (0,),
    deleters: Iterable[int] = (1,),
) -> FunctionalGroupFactory:
    return FunctionalGroupFactory("[*][Br]", tuple(bonders), tuple(deleters))
