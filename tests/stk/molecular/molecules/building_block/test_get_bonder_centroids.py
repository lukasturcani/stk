import numpy as np
import itertools as it


def test_get_bonder_centroids(building_block, get_fg_ids):
    bonder_centroids = list(get_bonder_centroids(
        building_block=building_block,
        fg_ids=get_fg_ids(building_block),
    ))
    position_matrix = get_position_matrix(
        building_block=building_block,
        fg_ids=get_fg_ids(building_block),
        bonder_centroids=bonder_centroids,
    )
    building_block = building_block.with_position_matrix(
        position_matrix=position_matrix,
    )
    result = it.zip_longest(
        building_block.get_bonder_centroids(
            fg_ids=get_fg_ids(building_block),
        ),
        bonder_centroids,
    )
    for centroid, expected_centroid in result:
        assert np.allclose(
            a=centroid,
            b=expected_centroid,
            atol=1e-32,
        )


def get_bonder_centroids(building_block, fg_ids):
    generator = np.random.RandomState(4)
    for fg in building_block.get_functional_groups(fg_ids):
        yield generator.normal(scale=50, size=3)


def get_position_matrix(building_block, fg_ids, bonder_centroids):
    position_matrix = building_block.get_position_matrix()
    centroids = zip(
        building_block.get_functional_groups(fg_ids),
        bonder_centroids,
    )
    generator = np.random.RandomState(5)
    for fg, centroid in centroids:
        set_centroid(position_matrix, fg, centroid, generator)
    return position_matrix


def set_centroid(
    position_matrix,
    functional_group,
    centroid,
    generator,
):
    bonder_ids = tuple(functional_group.get_bonder_ids())
    position_matrix[bonder_ids, :] = centroid
    for i in range(0, len(bonder_ids)-1, 2):
        displacement = generator.normal(scale=10, size=3)
        position_matrix[bonder_ids[i]] += displacement
        position_matrix[bonder_ids[i+1]] -= displacement
