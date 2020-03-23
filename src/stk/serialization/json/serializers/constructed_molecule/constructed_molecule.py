from stk.molecular import InchiKey, ConstructedMoleculeKeyMaker


class ConstructedMoleculeJsonizer:
    """
    Serializes :class:`.ConstructedMolecule` instances into JSON.

    """

    def __init__(
        self,
        key_makers=(
            InchiKey(),
            ConstructedMoleculeKeyMaker(),
        ),
    ):
        """
        Initializes a :class:`.ConstructedMoleculeJsonizer`.

        Parameters
        ----------
        key_makers : :class:`tuple` of \
                :class:`.MoleculeKeyMaker` and \
                :class:`.ConstructedMoleculeKeyMaker`
            Used to make the keys of molecules, which should be
            included in their JSON representations.

        """

        self._key_makers = key_makers

    def to_json(self, molecule):
        """
        Serialize `molecule`.

        Parameters
        ----------
        molecule : :class:`.ConstructedMolecule`
            The constructed molecule to serialize.

        Returns
        -------
        :class:`dict`
            A JSON representation of `molecule`.

        """

        def get_keys(building_block):
            return {
                key_maker.get_key_name(): key_maker.get_key(molecule)
                for key_maker in self._key_makers
            }

        return {
            'building_blocks': tuple(map(
                get_keys,
                molecule.get_building_blocks(),
            )),
            'atom_infos': tuple(map(
                self._atom_info_to_json,
                molecule.get_atom_infos(),
            )),
            'bond_infos': tuple(map(
                self._bond_info_to_json,
                enumerate(molecule.get_bond_infos()),
            )),
            'num_building_blocks': {
                index: molecule.get_num_building_block(building_block)
                for index, building_block
                in enumerate(molecule.get_building_blocks())
            },
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
            'BuildingBlock': self._building_block_key(
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
            'BuildingBlock': self._building_block_key(
                building_block=bond_info.get_building_block(),
            ),
            'building_block_id': bond_info.get_building_block_id(),
        }
