"""
Constructed Molecule Key
========================

"""

from ..molecule_keys import InchiKey
from ..building_block_keys import BuildingBlockKey


class ConstructedMoleculeKey:
    """
    An abstract base class for :class:`.ConstructedMolecule` keys.

    """

    def __init__(
        self,
        name='ConstructedMoleculeKey',
        molecule_key=InchiKey(),
        building_block_key=BuildingBlockKey(),
        topology_graph_keys=None,
    ):
        """
        Initialize a :class:`.ConstrcutedMoleculeKey` instance.

        Parameters
        ----------
        name : :class:`str`, optional

        molecule_key : :class:`.MoleculeKey`, optional

        building_block_key : :class:`.BuildingBlockKey`, optional

        topology_graph_keys : :class:`dict`, optional

        """

        if topology_graph_keys is None:
            topology_graph_keys = {}

        self._name = name
        self._molecule_key = molecule_key
        self._building_block_key = building_block_key
        self._topology_graph_keys = (
            self._get_default_topology_graph_keys()
        )
        self._topology_graph_key.update(topology_graph_keys)

    def get_name(self):
        """
        Get the name of the key.

        Returns
        -------
        :class:`str`
            The name.

        """

        return self._name

    def get_key(self, constructed_molecule):
        """
        Get the key of `constructed_molecule`.

        Parameters
        ----------
        constructed_molecule : :class:`.ConstructedMolecule`
            The constructed molecule for which a key is wanted.

        Returns
        -------
        :class:`object`

        """

        return (
        )

    @staticmethod
    def _get_default_topology_graph_keys():
        """
        Get the default topology graph key for each topology graph.

        Returns
        -------
        :class:`dict`

        """

        return {
        }
