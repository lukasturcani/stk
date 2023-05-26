import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        lambda: CaseData(
            key_maker=stk.Inchi(),
            molecule=stk.BuildingBlock("NCCN"),
            key_name="InChI",
            key="InChI=1S/C2H8N2/c3-1-2-4/h1-4H2",
        ),
        lambda: CaseData(
            key_maker=stk.InchiKey(),
            molecule=stk.BuildingBlock("NCCN"),
            key_name="InChIKey",
            key="PIICEJLVQHRZGT-UHFFFAOYSA-N",
        ),
        lambda: CaseData(
            key_maker=stk.Smiles(),
            molecule=stk.BuildingBlock("NCCN"),
            key_name="SMILES",
            key="NCCN",
        ),
        lambda: CaseData(
            key_maker=stk.Smiles(),
            molecule=stk.BuildingBlock("C(N)CN"),
            key_name="SMILES",
            key="NCCN",
        ),
        lambda: CaseData(
            key_maker=stk.Smiles(),
            molecule=stk.BuildingBlock(
                "C(#Cc1cccc2cnccc12)c1ccc2[nH]c3ccc"
                "(C#Cc4cccc5ccncc45)cc3c2c1"
            ),
            key_name="SMILES",
            key=(
                "C(#Cc1cccc2cnccc12)c1ccc2[nH]c3ccc"
                "(C#Cc4cccc5ccncc45)cc3c2c1"
            ),
        ),
        lambda: CaseData(
            key_maker=stk.Smiles(),
            molecule=stk.BuildingBlock(
                "C(#Cc1cccc2cnccc21)c1ccc2[nH]c3ccc"
                "(C#Cc4cccc5ccncc54)cc3c2c1"
            ),
            key_name="SMILES",
            key=(
                "C(#Cc1cccc2cnccc12)c1ccc2[nH]c3ccc"
                "(C#Cc4cccc5ccncc45)cc3c2c1"
            ),
        ),
        lambda: CaseData(
            key_maker=stk.Smiles(),
            molecule=stk.BuildingBlock("C[C@H](O)c1ccccc1"),
            key_name="SMILES",
            key="C[C@H](O)c1ccccc1",
        ),
        lambda: CaseData(
            key_maker=stk.MoleculeKeyMaker(
                key_name="NumAtoms",
                get_key=lambda molecule: molecule.get_num_atoms(),
            ),
            molecule=stk.BuildingBlock("NCCN"),
            key_name="NumAtoms",
            key=12,
        ),
    ),
)
def case_data(request) -> CaseData:
    return request.param()
