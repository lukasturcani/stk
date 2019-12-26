

class TestGetPlaneNormal:
    def case1(atom_ids, normal):
        molecule = stk.BuildingBlock('NCCN')
        coords = molecule.get_position_matrix()
        coords[[1, 9], 2] = 0
        molecule = molecule.with_position_matrix(coords)
        return molecule, atom_ids, normal

    def case2(atom_ids, normal):
        molecule = stk.BuildingBlock('NCCN')
        coords = molecule.get_position_matrix()
        coords[:, 2] = 0
        molecule = molecule.with_position_matrix(coords)
        return molecule, atom_ids, normal

    @pytest.mark.parametrize(
        'molecule,atom_ids,normal',
        [
            case1(atom_ids=(1, 9), normal=[0, 0, 1]),
            case2(atom_ids=None, normal=[0, 0, 1]),
            case2(atom_ids=(1, 9), normal=[0, 0, 1]),
        ],
    )
    def test_1(self, molecule, atom_ids, normal):
        assert np.all(np.equal(
            molecule.get_plane_normal(atom_ids),
            normal,
        ))

    @pytest.mark.parametrize(
        'molecule,atom_ids,normal',
        [
            case1(atom_ids=None, normal=[0, 0, 1])
        ],
    )
    def test_2(self, molecule, atom_ids, normal):
        assert not np.all(np.equal(
            molecule.get_plane_normal(atom_ids),
            normal,
        ))
