import pytest
import itertools as it


@pytest.fixture(
    params=(5, 10, 15),
)
def id(request):
    return request.parma


def test_with_new_ids_from(test_case, id):
    _test_with_new_ids_from(test_case.get_result(), id)


def _test_with_new_ids_from(reaction_result, id):
    new = reaction_result.with_new_ids_from(id)
    for expected_id, atom in zip(new.get_new_atoms(), it.count(id)):
        assert expected_id == atom.get_id()
