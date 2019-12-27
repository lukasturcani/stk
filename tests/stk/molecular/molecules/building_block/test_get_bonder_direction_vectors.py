
class TestGetBonderDirectionVectors:
    def case1():
        building_block = stk.BuildingBlock(
            smiles='BrCC(CBr)(CBr)CBr',
            functional_groups=['bromine'],
        )
        coords = building_block.get_position_matrix()
        bonders = enumerate(building_block.get_bonder_ids())
        for fg_id, bonder_id in bonders:
            coords[bonder_id] = [fg_id, fg_id, fg_id]
        building_block.set_position_matrix(coords)
        expected_vectors = {
            (1, 0): [-1, -1, -1],
            (2, 0): [-2, -2, -2],
            (3, 0): [-3, -3, -3],
            (2, 1): [-1, -1, -1],
            (3, 1): [-2, -2, -2],
            (3, 2): [-1, -1, -1],
        }
        return building_block, None, expected_vectors

    @pytest.mark.parametrize(
        argnames=(
            'building_block',
            'fg_ids',
            'expected_vectors',
        ),
        argvalues=(
            case1(),
        ),
    )
    def test(
        self,
        building_block,
        fg_ids,
        expected_vectors,
    ):
        direction_vectors = (
            building_block.get_bonder_direction_vectors(fg_ids=fg_ids)
        )
        for fg1, fg2, direction in direction_vectors:
            assert np.allclose(
                a=direction,
                b=expected_vectors[(fg1, fg2)],
                atol=1e-32,
            )
