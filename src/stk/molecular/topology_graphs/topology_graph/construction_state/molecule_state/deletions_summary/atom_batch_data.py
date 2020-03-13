from .....atoms import AtomInfo


class _AtomBatchData:

    __slots__ = ['_atoms', '_atom_infos', '_atom_map']

    def __init__(
        self,
        atoms,
        num_atoms,
        building_block,
        building_block_id,
    ):
        self._atoms = _atoms = []
        self._atom_infos = atom_infos = []
        self._atom_map = atom_map = {}

        for id_, atom in enumerate(atoms, num_atoms):
            _atoms.append(atom.with_id(id_))
            atom_map[atom.get_id()] = _atoms[-1]
            atom_infos.append(
                AtomInfo(
                    atom=_atoms[-1],
                    building_block=building_block,
                    building_block_id=building_block_id,
                ),
            )

    def get_atoms(self):
        yield from self._atoms

    def get_atom_infos(self):
        yield from self._atom_infos

    def get_atom_map(self):
        return dict(self._atom_map)
