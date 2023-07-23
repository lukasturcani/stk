import pathlib

import rdkit.Chem.AllChem as rdkit
import stk


def test_four_plus_six(tmp_path: pathlib.Path) -> None:
    graph = stk.cage.four_plus_six()
    cage = graph.build(
        primary=stk.cage.BuildingBlock.from_smiles(
            smiles="BrCc1cc(CBr)cc(CBr)c1",
            functional_groups=stk.bromo(),
        ),
        secondary=stk.cage.BuildingBlock.from_smiles(
            smiles="BrCC(CC)CBr",
            functional_groups=stk.bromo(),
        ),
    )
    rdkit.MolToMolFile(cage.to_rdkit(), str(tmp_path / "cage.mol"))
