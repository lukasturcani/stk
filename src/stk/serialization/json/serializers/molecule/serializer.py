from stk.molecular import (
    Molecule,
    BuildingBlock,
    ConstructedMolecule,
)
from .utilities import (
    molecule_to_json,
    _BuildingBlockSerializer,
    _ConstructedMoleculeSerializer,
    molecule_key,
    building_block_key,
)
from ....molecule_serializer import MoleculeSerializer


class JsonMoleculeSerializer(MoleculeSerializer):
    """

    Examples
    --------
    *Usage*

    *Using Different Serializers*

    *Adding Additional Serializers*


    """

    def __init__(self, serializers=None):
        """
        Initialize a :class:`.JsonMoleculeSerializer` instance.

        Parameters
        ----------
        serializers : :class:`dict`, optional
            Maps :class:`.Molecule` and each of its subclasses to the
            :class:`callable`, which should be used to create
            the JSON representation for that class.

        """

        if serializers is None:
            serializers = self.get_default_serializers()

        self._serializers = serializers

    def serialize(self, molecule):
        return {
            cls.__name__: serializer.serialize(molecule)
            for cls, serializer in self._serializers.items()
            if isinstance(molecule, cls)
        }

    @classmethod
    def get_default_serializers(self):
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
            Molecule: molecule_to_json,
            BuildingBlock: _BuildingBlockSerializer(
                molecule_key=molecule_key,
                functional_group_serializer=...,
            ),
            ConstructedMolecule: _ConstructedMoleculeSerializer(
                molecule_key=molecule_key,
                building_block_key=building_block_key,
                topology_graph_serializer=...,
            ),
        }
