import stk


class CaseData:
    """
    A test case.

    Attributes
    ----------
    molecule : :class:`.Molecule`
        The molecule to be tested.

    smiles : :class:`str`
        The canonical smiles for the molecule.

    """

    def __init__(self, molecule, smiles):
        self.molecule = molecule
        self.smiles = smiles

    @classmethod
    def init_constructed_molecule(
        cls,
        building_blocks,
        topology_graph,
        building_block_vertices,
        smiles,
    ):
        return cls(
            molecule=stk.ConstructedMolecule(
                building_blocks=building_blocks,
                topology_graph=topology_graph,
                building_block_vertices={
                    building_blocks[building_block_index]:
                        tuple(topology_graph.get_vertices(vertex_ids))
                    for building_block_index, vertex_ids
                    in building_block_vertices.items()
                },
            ),
            smiles=smiles,
        )
