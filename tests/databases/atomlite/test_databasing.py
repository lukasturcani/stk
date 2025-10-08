import pathlib

import atomlite
import numpy as np

import stk

from .case_data import CaseData


def test_atomlite(molecule: CaseData) -> None:
    """Test :class:`.AtomliteDatabase`."""
    path = pathlib.Path(__file__).resolve().parent / "test.db"

    if path.exists():
        path.unlink()

    database = atomlite.Database(path)
    assert database.num_entries() == 0
    for mol, prop in zip(
        molecule.molecules,
        molecule.property_dicts,
        strict=True,
    ):
        key = stk.Smiles().get_key(mol)

        database.add_entries(
            atomlite.Entry(
                key=key,
                molecule=mol.to_atomlite(),
                properties=prop,
            )
        )

        entry = database.get_entry(key)

        assert np.allclose(
            mol.get_position_matrix(),
            mol.with_structure_from_atomlite(
                database=database,
                key=key,
            ).get_position_matrix(),
        )

        assert entry.molecule == atomlite.json_from_rdkit(mol.to_rdkit_mol())
        assert entry.molecule == mol.to_atomlite()
        assert (
            stk.Smiles().get_key(
                mol.init_from_atomlite(database=database, key=key)
            )
            == key
        )

        assert database.get_property(key, path="$.1") == prop["1"]
        assert isinstance(database.get_property(key, path="$.1"), int)
        assert database.get_property(key, path="$.3") == prop["3"]
        assert isinstance(database.get_property(key, path="$.3"), str)
        if "2" in prop:
            assert database.get_property(key=key, path="$.2") == prop["2"]
            assert isinstance(database.get_property(key, path="$.2"), dict)

    assert database.num_entries() == molecule.expected_count

    path.unlink()
