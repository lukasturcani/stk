def test_get_topology_graph(case_data):
    assert (
        case_data.constructed_molecule.get_topology_graph()
        is case_data.topology_graph
    )
