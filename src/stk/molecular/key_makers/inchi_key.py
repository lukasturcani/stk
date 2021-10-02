"""
InChIKey
========

"""

from .molecule import MoleculeKeyMaker
from .utilities import get_inchi_key

__all__ = (
    'InchiKey',
)


class InchiKey(MoleculeKeyMaker):
    """
    Used to get the InChIKey of molecules.

    Examples:

        *Adding InChIKey to a Molecule's JSON*

        You want to use the InChIKey as part of a JSON representation
        of a molecule

        .. testcode:: adding-inchikey-to-a-molecules-json

            import stk

            jsonizer = stk.MoleculeJsonizer(
                key_makers=(stk.InchiKey(), ),
            )
            # Get the JSON representation, including an InChIKey.
            json = jsonizer.to_json(stk.BuildingBlock('NCCN'))

        .. testcode:: adding-inchikey-to-a-molecules-json
            :hide:

            assert (
                json['molecule']['InChIKey']
                == 'PIICEJLVQHRZGT-UHFFFAOYSA-N'
            )
            assert (
                json['molecule']['InChIKey']
                == 'PIICEJLVQHRZGT-UHFFFAOYSA-N'
            )

    """

    def __init__(self) -> None:
        """
        Initialize an :class:`InchiKey` instance.

        """

        MoleculeKeyMaker.__init__(
            self=self,
            key_name='InChIKey',
            get_key=get_inchi_key,
        )

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return 'InchiKey()'
