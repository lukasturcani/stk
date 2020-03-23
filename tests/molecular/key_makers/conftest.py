import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            key_maker=stk.Inchi(),
            molecule=stk.BuildingBlock('NCCN'),
            key_name='InChI',
            key='InChI=1S/C2H8N2/c3-1-2-4/h1-4H2',
        ),
        CaseData(
            key_maker=stk.InchiKey(),
            molecule=stk.BuildingBlock('NCCN'),
            key_name='InChIKey',
            key='PIICEJLVQHRZGT-UHFFFAOYSA-N',
        ),
        CaseData(
            key_maker=stk.MoleculeKeyMaker(
                key_name='NumAtoms',
                get_key=lambda molecule: molecule.get_num_atoms(),
            ),
            molecule=stk.BuildingBlock('NCCN'),
            key_name='NumAtoms',
            key=12,
        ),
        CaseData(
            key_maker=stk.ConstructedMoleculeKeyMaker(),
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='NCCN',
                            functional_groups=[
                                stk.PrimaryAminoFactory()
                            ],
                        ),
                    ),
                    repeating_unit='A',
                    num_repeating_units=2,
                ),
            ),
            key_name='ConstructedMoleculeKey',
            key=str((
                'KCEQCQJHUVIRJO-UHFFFAOYSA-N',
                ('PIICEJLVQHRZGT-UHFFFAOYSA-N', ),
                (2, ),
            )),
        ),
        CaseData(
            key_maker=stk.ConstructedMoleculeKeyMaker(
                key_name='MyKey',
                molecule_key_maker=stk.Inchi(),
            ),
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='NCCN',
                            functional_groups=[
                                stk.PrimaryAminoFactory()
                            ],
                        ),
                    ),
                    repeating_unit='A',
                    num_repeating_units=2,
                ),
            ),
            key_name='MyKey',
            key=str((
                'InChI=1S/C4H14N4/c5-1-3-7-8-4-2-6/h7-8H,1-6H2',
                ('InChI=1S/C2H8N2/c3-1-2-4/h1-4H2', ),
                (2, ),
            )),
        ),
    ),
)
def case_data(request):
    return request.param
