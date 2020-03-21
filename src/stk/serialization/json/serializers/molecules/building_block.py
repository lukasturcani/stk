from .utilities import molecule_to_json
from ..functional_group import FunctionalGroupSerializer


class _BuildingBlockSerializer:
    """
    Serializes :class:`.BuildingBlock` instances into JSON.

    Notes
    -----
    This class is an implementation detail of
    :class:`.MoleculeSerializer`, use that class directly for your
    serialization needs.

    """

    def __init__(
        self,
        functional_group_serializer=FunctionalGroupSerializer(),
    ):
        """
        Initialize a :class:`._BuildingBlockSerializer` instance.

        Parameters
        ----------
        functional_group_serializer : :class:`.FunctionalGroupSerializer`
            Used to serialize the functional groups of the serialized
            building blocks.

        """

        return super().__init__()

    def serialize(self, molecule):
        """
        Serialize `molecule`.

        Parameters
        -----------
        molecule : :class:`.BuildingBlock`
            The molecule to serialize.

        Returns
        -------
        :class:`dict`
            A JSON representation of `molecule`.

        """

        return {
            'molecule': molecule_to_json(molecule),
            'functional_groups': tuple(map(
                self._functional_group_serializer.serialize,
                molecule.get_functional_groups(),
            )),
            'placer_ids': tuple(molecule.get_placer_ids()),
        }

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._functional_group_serializer!r})'
        )
