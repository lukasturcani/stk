

class TestWithStructureFromFile1:
    @pytest.fixture(params=[
        'molecule.mol',
        'molecule.xyz',
    ])
    def path(self, request, tmpdir):
        return os.path.join(tmpdir, request.param)

    def case1():
        conformer1 = stk.BuildingBlock('NCCN')
        conformer2 = stk.BuildingBlock('NCCN')
        conformer2 = conformer2.with_position_matrix(
            position_matrix=np.zeros((conformer2.get_num_atoms(), 3)),
        )
        return conformer1, conformer2

    def case2():
        bb = stk.BuildingBlock('BrCCBr', ['bromine'])
        topology_graph = stk.polymer.Linear('A', 3)
        conformer1 = stk.ConstructedMolecule(
            building_blocks=[bb],
            topology_graph=topology_graph,
        )
        conformer2 = stk.ConstructedMolecule(
            building_blocks=[bb],
            topology_graph=topology_graph,
        )
        conformer2.set_position_matrix(
            position_matrix=np.zeros((conformer2.get_num_atoms(), 3)),
        )
        return conformer1, conformer2

    @pytest.mark.parametrize(
        'get_conformers',
        [
            case1,
            # case2,
        ],
    )
    def test(self, get_conformers, path):
        conformer1, conformer2 = get_conformers()
        assert not np.allclose(
            a=conformer1.get_position_matrix(),
            b=conformer2.get_position_matrix(),
            atol=1e-4,
        )

        conformer2.write(path)
        # Keep a clone for immutability testing.
        clone = conformer1.clone()
        new = conformer1.with_structure_from_file(path)

        assert np.allclose(
            a=new.get_position_matrix(),
            b=conformer2.get_position_matrix(),
            atol=1e-4,
        )
        _test_unchanged(clone, conformer1)


class TestWithStructureFromFile2:
    @pytest.fixture
    def path(self, datadir, request):
        return str(datadir / request.param)

    def case1():
        return stk.BuildingBlock('NCCN'), 'NCCN.mae'

    @pytest.mark.parametrize(
        'molecule,path',
        [
            case1(),
        ],
        indirect=['path'],
    )
    def test(self, molecule, path):
        # Keep clone for immutability testing.
        clone = molecule.clone()
        new = molecule.with_structure_from_file(path)
        size_diff = abs(
            molecule.get_maximum_diameter()
            - new.get_maximum_diameter()
        )
        assert size_diff > 1
        _test_unchanged(clone, molecule)
