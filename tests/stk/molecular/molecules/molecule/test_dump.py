from .utilities import has_same_structure, is_equivalent_molecule


def test_dump(molecule, tmpdir):
    path = str(tmpdir / 'molecule.dump')
    molecule.dump(path)
    new = molecule.__class__.load(path)
    is_equivalent_molecule(molecule, new)
    has_same_structure(molecule, new)
