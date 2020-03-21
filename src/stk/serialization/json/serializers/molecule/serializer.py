"""
JSON Molecule Serializer
========================

"""

from stk.molecular import (
    Molecule,
    BuildingBlock,
    ConstructedMolecule,
)
from ..functional_group import FunctionalGroupJsonizer
from ..topology_graph import TopologyGraphJsonizer
from .utilities import (
    MoleculeJsonizer,
    BuildingBlockJsonizer,
    ConstructedMoleculeJsonizer,
    get_inchi_key,
    GetBuildingBlockKey,
    GetConstructedMoleculeKey,
)
from ....molecule_serializer import MoleculeSerializer


class JsonMoleculeSerializer(MoleculeSerializer):
    """
    Serializes :class:`.Molecule` instances into JSON.

    Examples
    --------
    *Usage*

    *Using Different Serializers*

    *Adding Additional Serializers*


    """

    def __init__(self, jsonizers=None):
        """
        Initialize a :class:`.JsonMoleculeSerializer` instance.

        Parameters
        ----------
        jsonizers : :class:`dict`, optional
            Maps :class:`.Molecule` and each of its subclasses to the
            :class:`callable`, which should be used to create
            the JSON representation for that class.

        """

        if jsonizers is None:
            jsonizers = self.get_default_jsonizers()

        self._jsonizers = jsonizers

    def serialize(self, molecule):
        return {
            cls.__name__: jsonizer.to_json(molecule)
            for cls, jsonizer in self._jsonizers.items()
            if isinstance(molecule, cls)
        }

    @staticmethod
    def get_default_serializers():
        """
        Return the default serializers.

        Returns
        -------
        :class:`dict`
            Maps :class:`.Molecule` and each of its subclasses to the
            default :class:`callable`, which should be used to create
            the JSON representation for that class.

        """

        return {
            Molecule: MoleculeJsonizer(
                molecule_key=get_inchi_key,
            ),
            BuildingBlock: BuildingBlockJsonizer(
                molecule_key=get_inchi_key,
                building_block_key=GetBuildingBlockKey(),
                functional_group_jsonizer=FunctionalGroupJsonizer(),
            ),
            ConstructedMolecule: ConstructedMoleculeJsonizer(
                molecule_key=get_inchi_key,
                building_block_key=GetBuildingBlockKey(),
                constructed_molecule_key=GetConstructedMoleculeKey(),
                topology_graph_jsonizer=TopologyGraphJsonizer(),
            ),
        }
