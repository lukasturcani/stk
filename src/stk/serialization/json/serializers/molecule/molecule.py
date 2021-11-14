"""
Molecule JSONizer
=================

"""

from stk.molecular import InchiKey

from .utilities import atom_to_json, bond_to_json


class MoleculeJsonizer:
    """
    Abstract base class for creating JSONs of molecules.

    See Also
    --------
    :class:`.ConstructedMoleculeJsonizer`

    Notes
    -----
    You might notice that the public methods of this abstract base
    class are implemented. These are just default implementations,
    which can be safely ignored or overridden, when implementing
    subclasses. However, the default implementation can be used
    directly, if it suits your needs.

    Examples
    --------
    *Converting a Molecule to JSON*

    You want to create a JSON representation of a molecule

    .. testcode:: converting-a-molecule-to-json

        import stk

        jsonizer = stk.MoleculeJsonizer()
        json = jsonizer.to_json(stk.BuildingBlock('NCCN'))

    *Adding Additional Molecular Keys*

    Apart from atoms, bonds and the position matrix, the JSON
    representation holds additional fields, one for each
    :class:`.MoleculeKeyMaker` provided to the initializer

    .. testcode:: adding-additional-molecular-keys

        import stk

        jsonizer = stk.MoleculeJsonizer(
            key_makers=(
                stk.Inchi(),
                stk.InchiKey(),
            ),
        )
        json = jsonizer.to_json(stk.BuildingBlock('NCCN'))

    In this case, ``json`` will look something like

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
            included in their JSON representations. Keys allow
            molecular data to reference itself when split across
            multiple JSONs.

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
            'a': tuple(map(atom_to_json, molecule.get_atoms())),
            'b': tuple(map(bond_to_json, molecule.get_bonds())),
        }
        position_matrix = {
            'm': molecule.get_position_matrix().tolist(),
        }
        for key_maker in self._key_makers:
            key_name = key_maker.get_key_name()
            key = key_maker.get_key(molecule)
            json[key_name] = key
            position_matrix[key_name] = key
        return {
            'molecule': json,
            'matrix': position_matrix,
        }

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'{self.__class__.__name__}({self._key_makers!r})'
