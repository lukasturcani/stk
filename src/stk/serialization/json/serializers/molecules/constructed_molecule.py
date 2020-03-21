class _ConstructedMoleculeSerializer:
    """
    Serializes :class:`.ConstructedMolecule` instances into JSON.

    Notes
    -----
    This class is an implementation detail of
    :class:`.MoleculeSerializer`, use that class directly for your
    serialization needs.

    """

    def __init__(
        self,
        molecule_key,
        building_block_key,
        topology_graph_serializer,
    ):
        """
        Initializes a :class:`._ConstructedMoleculeSerializer`.

        Parameters
        ----------
        molecule_key : :class:`callable`
            Takes a single parameter, `molecule`, and returns
            a key used for referencing that molecule. The parameter
            requires a :class:`.Molecule` instance.

        building_block_key : :class:`callable`
            Takes a single parameter, `building_block`, and returns
            a key used for referencing that building block. The
            parameter requires a :class:`.BuildingBlock` instance.

        topology_graph_serializer : \
                :class:`.TopologyGraphSerializer`
            Used to serialize the topology graph of the serialized
            constructed molecules.

        """

        self._molecule_key = molecule_key
        self._building_block_key = building_block_key
        self._topology_graph_serializer = topology_graph_serializer

    def serialize(self, constructed_molecule):
        """
        Serialize : :class:`molecule`.

        Parameters
        ----------
        constructed_molecule : :class:`.ConstructedMolecule`
            The constructed molecule to serialize.

        Returns
        -------
        :class:`dict`
            A JSON representation of `molecule`.

        """

        return {
            'molecule': self._molecule_key(constructed_molecule),
            'topology_graph':
                self._topology_graph_serializer.serialize(
                    # Access to the topology graph is not part of the
                    # public interface of a ConstructedMolecule so
                    # private attribute access must be used.
                    topology_graph=(
                        constructed_molecule._topology_graph
                    ),
                ),
            'atom_infos': tuple(map(
                self._atom_info_to_json,
                constructed_molecule.get_atom_infos(),
            )),
            'bond_infos': tuple(map(
                self._bond_info_to_json,
                enumerate(constructed_molecule.get_bond_infos()),
            )),
        }

    def _atom_info_to_json(self, atom_info):
        """
        Return a JSON representation of `atom_info`.

        Parameters
        ----------
        atom_info : :class:`.AtomInfo`
            The atom info to serialize.

        Returns
        -------
        :class:`dict`
            A JSON representation of `atom_info`.

        """

        return {
            'atom': atom_info.get_atom().get_id(),
            'building_block': self._building_block_key(
                building_block=atom_info.get_building_block(),
            ),
            'building_block_id': atom_info.get_building_block_id(),
        }

    def _bond_info_to_json(self, bond_info_data):
        """
        Return a JSON representation of `bond_info_data`.

        Parameters
        ----------
        bond_info_data : :class:`tuple`
            The first element of the :class:`tuple` is the id of the
            bond and second element of the :class:`tuple` is the
            bond info of the bond.

        Returns
        -------
        :class:`dict`
            A JSON representation of `bond_info_data`.

        """

        bond_id, bond_info = bond_info_data
        return {
            'bond': bond_id,
            'building_block': self._building_block_key(
                building_block=bond_info.get_building_block(),
            ),
            'building_block_id': bond_info.get_building_block_id(),
        }
