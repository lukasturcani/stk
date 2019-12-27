import stk
import pytest


class _TestBonderDistances:
    def case1():
        building_block = stk.BuildingBlock(
            smiles='BrCC(CBr)(CBr)CBr',
            functional_groups=['bromine'],
        )
        coords = building_block.get_position_matrix()
        bonders = enumerate(building_block.get_bonder_ids())
        for fg, bonder_id in bonders:
            coords[bonder_id] = [fg, 0, 0]
        building_block.set_position_matrix(coords)
        distances = {
            (0, 1): 1,
            (0, 2): 2,
            (0, 3): 3,
            (1, 2): 1,
            (1, 3): 2,
            (2, 3): 1,
        }
        return building_block, None, distances

    @pytest.mark.parametrize(
        argnames=(
            'building_block',
            'fg_ids',
            'expected_distances',
        ),
        argvalues=(
            # case1(),
        ),
    )
    def test(
        self,
        building_block,
        fg_ids,
        expected_distances,
    ):
        distances = building_block.get_bonder_distances(
            fg_ids=fg_ids,
        )
        for fg1, fg2, distance in distances:
            assert (
                abs(expected_distances[(fg1, fg2)] - distance) < 1e-32
            )
