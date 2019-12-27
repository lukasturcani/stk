


class TestGetBonderPlaneNormal:
    def case1():
        building_block = stk.BuildingBlock(
            smiles='BrCC(CBr)(CBr)CBr',
            functional_groups=['bromine'],
        )
        bonder_ids = list(building_block.get_bonder_ids())
        coords = building_block.get_position_matrix()
        coords[bonder_ids[0]] = [1, 1, 0]
        coords[bonder_ids[1]] = [0, 0, 0.5]
        coords[bonder_ids[2]] = [0, 0, -0.5]
        coords[bonder_ids[3]] = [1, -1, 0]
        building_block.set_position_matrix(coords)
        return building_block, None, [0, 0, -1]

    @pytest.mark.parametrize(
        argnames=(
            'building_block',
            'fg_ids',
            'expected_bonder_plane_normal',
        ),
        argvalues=(
            case1(),
        ),
    )
    def test(
        self,
        building_block,
        fg_ids,
        expected_bonder_plane_normal,
    ):
        bonder_plane_normal = building_block.get_bonder_plane_normal(
            fg_ids=fg_ids,
        )
        assert np.allclose(
            a=bonder_plane_normal,
            b=expected_bonder_plane_normal,
            atol=1e-32,
        )
