class _BuildingBlockSerializer:
    """
    Serializes :class:`.BuildingBlock` instances into JSON.

    Notes
    -----
    This class is an implementation detail of
    :class:`.MoleculeSerializer`, use that class directly for your
    serialization needs.

    """

    def __init__(self, molecule_key, functional_group_serializer):
        """
        Initialize a :class:`._BuildingBlockSerializer` instance.

        Parameters
        ----------
        molecule_key : :class:`callable`
            Takes a single parameter, `molecule`, and returns
            a key used for referencing that molecule. The parameter
            requires a :class:`.Molecule` instance.

        functional_group_serializer : \
                :class:`.FunctionalGroupSerializer`
            Used to serialize the functional groups of the serialized
            building blocks.

        """

        self._molecule_key = molecule_key
        self._functional_group_serializer = functional_group_serializer

    def serialize(self, building_block):
        """
        Serialize `molecule`.

        Parameters
        -----------
        building_block : :class:`.BuildingBlock`
            The building block to serialize.

        Returns
        -------
        :class:`dict`
            A JSON representation of `molecule`.

        """

        return {
            'molecule': self._molecule_key(building_block),
            'functional_groups': tuple(map(
                self._functional_group_serializer.serialize,
                building_block.get_functional_groups(),
            )),
            'placer_ids': tuple(building_block.get_placer_ids()),
        }

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._functional_group_serializer!r})'
        )
