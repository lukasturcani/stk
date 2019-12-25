import numpy as np


def test_get_maximum_diameter(
    molecule,
    get_atom_ids,
    maximum_diameter,
):
    position_matrix = get_position_matrix(
        molecule=molecule,
        atom_ids=get_atom_ids(molecule),
        maximum_diameter=maximum_diameter,
    )
    molecule = molecule.with_position_matrix(position_matrix)
    assert np.allclose(
        a=maximum_diameter,
        b=molecule.get_maximum_diameter(get_atom_ids(molecule)),
        atol=1e-32,
    )


def get_position_matrix(molecule, atom_ids, maximum_diameter):
    """
    Create a position matrix with a specific `maximum_diameter`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule for which the position matrix is created.

    atom_ids : :class:`iterable` of :class:`int`
        The ids of atoms which are to have a specific maximum
        diameter. If ``None``, then all atoms are used.

    maximum_diameter : :class:`float`
        The desired maximum diameter.

    Returns
    -------
    :class:`numpy.ndarray`
        The position matrix for `molecule`.

    """

    if atom_ids is None:
        atom_ids = range(molecule.get_num_atoms())
    elif not isinstance(molecule, (list, tuple)):
        atom_ids = tuple(atom_ids)

    position_matrix = molecule.get_position_matrix()



class TestGetMaximumDiameter:
    def case1():
        molecule = stk.BuildingBlock('NCCN')
        coords = np.array([
            [i, 0, 0] for i in range(molecule.get_num_atoms())
        ])
        molecule = molecule.with_position_matrix(coords)
        return molecule, None, molecule.get_num_atoms()-1

    def case2(atom_ids, maximum_diameter):
        molecule = stk.BuildingBlock('NCCN')
        coords = np.zeros((molecule.get_num_atoms(), 3))
        coords[[1]] = [0, -50, 0]
        coords[[9]] = [0, 50, 0]
        molecule = molecule.with_position_matrix(coords)
        return molecule, atom_ids, maximum_diameter

    @pytest.mark.parametrize(
        'molecule,atom_ids,maximum_diameter',
        [
            case1(),
            case2(atom_ids=None, maximum_diameter=100),
            case2(atom_ids=(1, 9), maximum_diameter=100),
            case2(atom_ids=(1, 0), maximum_diameter=50),
            case2(atom_ids=(0, 9), maximum_diameter=50),
            case2(atom_ids=(0, 2, 3, 4), maximum_diameter=0),
        ],
    )
    def test(
        self,
        molecule,
        atom_ids,
        maximum_diameter
    ):
        assert (
            molecule.get_maximum_diameter(atom_ids) == maximum_diameter
        )
