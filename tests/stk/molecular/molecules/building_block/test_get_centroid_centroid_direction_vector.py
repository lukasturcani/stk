import pytest
import stk
import numpy as np


class _TestGetCentroidCentroidDirectionVector:
    def case1():
        building_block = stk.BuildingBlock(
            smiles='BrCC(CBr)(CBr)CBr',
            functional_groups=['bromine'],
        )
        coords = np.zeros((len(building_block.atoms), 3))
        for bonder_id in building_block.get_bonder_ids():
            coords[bonder_id] = [1, 0, 0]
        building_block.set_position_matrix(coords)
        return building_block, None, [-1, 0, 0]

    @pytest.mark.parametrize(
        argnames=(
            'building_block',
            'fg_ids',
            'expected_vector',
        ),
        argvalues=(
            # case1(),
        )
    )
    def test(
        self,
        building_block,
        fg_ids,
        expected_vector,
    ):
        direction = (
            building_block.get_centroid_centroid_direction_vector(
                fg_ids=fg_ids,
            )
        )
        assert np.allclose(
            a=stk.normalize_vector(direction),
            b=expected_vector,
            atol=1e-32,
        )
