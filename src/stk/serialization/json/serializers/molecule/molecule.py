"""
Molecule JSONizer
=================

"""

from stk.molecular import InchiKey
from .utilities import atom_to_json, bond_to_json


class MoleculeJsonizer:
    """
    Creates JSON representations of :class:`.Molecule` instances.

    See Also
    --------
    :class:`.ConstructedMoleculeJsonizer`

    Examples
    --------
    You want to create a JSON representation of a molecule

    .. code-block:: python

        import stk

        jsonizer = stk.MoleculeJsonizer()
        json = jsonizer.to_json(stk.BuildingBlock('NCCN'))

    The apart from atoms and bonds, the JSON representation holds
    additional fields, one for each
    :class:`.MoleculeKeyMaker` provided to the initializer

    .. code-block:: python

        import stk

        jsonizer = stk.MoleculeJsonizer(
            key_makers=(
                stk.Inchi(),
                stk.InchiKey()
            ),
        )
        json = jsonizer.to_json(stk.BuildingBlock('NCCN'))

    In this case, ``json`` will have the form

    .. code-block:: python

        {
            # A tuple of JSON atom representations.
            'atoms': (...),

            # A tuple of JSON bond representations.
            'bonds': (...),

            'InChI': 'The InChI of the molecule',
            'InChIKey': 'The InChIKey of the molecule',
        }

    For every :class:`.MoleculeKeyMaker` provided to `key_makers`,
    a new key will be added to the JSON representation, with its name
    given by :meth:`.MoleculeKeyMaker.get_key_name` and the value
    given by :meth:`.MoleculeKeyMaker.get_key`.

    """

    def __init__(
        self,
        key_makers=(InchiKey(), ),
    ):
        """
        Initialize a :class:`.MoleculeJsonizer` instance.

        Parameters
        ----------
        key_makers : :class:`tuple` of \
                :class:`.MoleculeKeyMaker`
            Used to make the keys of molecules, which should be
            included in their JSON representations.

        """

        self._key_makers = key_makers

    def to_json(self, molecule):
        """
        Return a JSON representation of `molecule`.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            The molecule to serialize.

        Returns
        -------
        :class:`dict`
            A JSON representation of `molecule`.

        """

        json = {
            'atoms': tuple(map(atom_to_json, molecule.get_atoms())),
            'bonds': tuple(map(bond_to_json, molecule.get_bonds())),
        }
        for key_maker in self._key_makers:
            json[key_maker.get_key_name()] = (
                key_maker.get_key(molecule)
            )
        return json
