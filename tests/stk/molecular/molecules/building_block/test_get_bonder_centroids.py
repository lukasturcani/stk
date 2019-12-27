

class TestGetBonderCentroids:
    def case1():
        building_block = stk.BuildingBlock('BrCCBr', ['bromine'])
        position_matrix = np.zeros((len(building_block.atoms), 3))
        bonder1 = [1.0, 0., 0.]
        bonder2 = [10., 0., 0.]
        position_matrix[1, :] = bonder1
        position_matrix[2, :] = bonder2
        building_block.set_position_matrix(position_matrix)
        fg_ids = (0, 1)
        return building_block, fg_ids, [bonder1, bonder2]

    @pytest.mark.parametrize(
        argnames=(
            'building_block',
            'fg_ids',
            'expected_bonder_centroids',
        ),
        argvalues=(
            case1(),
        )
    )
    def test(
        self,
        building_block,
        fg_ids,
        expected_bonder_centroids,
    ):
        bonder_centroids = building_block.get_bonder_centroids(
            fg_ids=fg_ids,
        )
        bonder_centroids = it.zip_longest(
            bonder_centroids,
            expected_bonder_centroids,
        )
        for centroid, expected_bonder_centroid in bonder_centroids:
            assert np.allclose(
                a=centroid,
                b=expected_bonder_centroid,
                atol=1e-32,
            )
