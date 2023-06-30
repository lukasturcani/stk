import atomlite
import stk


def test_serde(
    constructed_molecule: stk.ConstructedMolecule,
) -> None:
    db = atomlite.Database(":memory:")
    db.add_entries(
        entries=constructed_molecule.to_atomlite_entry(
            key="test-molecule",
            properties={"a": "b"},
        ),
    )
    retrieved = stk.ConstructedMolecule.init_from_atomlite_entry(
        entry=next(db.get_entries("test-molecule")),
    )
