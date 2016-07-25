from ..classes import Population

def test_all_members():
    a = Population.init_empty()
    print(len(a))
    assert len(a) == 24