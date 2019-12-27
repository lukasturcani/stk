import stk
import pytest


class _BonderPlanePoint:
    def __init__(self, x, y, z, distance):
        self.x = x
        self.y = y
        self.z = z
        self.distance = distance


class TestGetBonderPlane:

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

        points = (
            _BonderPlanePoint(1, 1, 0, 0),
            _BonderPlanePoint(0, 0, 0.5, 0.5),
            _BonderPlanePoint(0, 0, -0.5, 0.5),
            _BonderPlanePoint(1, -1, 0, 0),
        )

        return building_block, None, points

    @pytest.mark.parametrize(
        argnames=(
            'building_block',
            'fg_ids',
            'points',
        ),
        argvalues=(
            # case1(),
        ),
    )
    def test(self, building_block, fg_ids, points):
        a, b, c, d = building_block.get_bonder_plane(
            fg_ids=fg_ids,
        )
        for point in points:
            product = a*point.x + b*point.y + c*point.z
            assert abs(point.distance - abs(product-d)) < 1e-16
