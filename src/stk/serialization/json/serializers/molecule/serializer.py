from stk.molecular import (
    Molecule,
    BuildingBlock,
    ConstructedMolecule,
)
from .molecule import _MoleculeSerializer
from .building_block import _BuildingBlockSerializer
from .constructed_molecule import _ConstructedMoleculeSerializer
from .utilities import (
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
        serializers : :class:`dict`

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
            default serializer, which should be used to create the JSON
            representation for that class.

        """

        return {
            Molecule: _MoleculeSerializer(),
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
