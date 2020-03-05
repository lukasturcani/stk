import stk


class CaseData:
    def __init__(
        self,
        building_blocks,
        topology_graph,
        num_new_atoms,
        num_new_bonds,
    ):
        self.building_block_vertices = (
            topology_graph.get_building_block_vertices(
                building_blocks=building_blocks,
            )
        )
        self.constructed_molecule = stk.ConstructedMolecule(
            building_blocks=building_blocks,
            topology_graph=topology_graph,
            building_block_vertices=self.building_block_vertices,
        )
        self.num_new_atoms = num_new_atoms
        self.num_new_bonds = num_new_bonds
        self.building_blocks_counts = {
            building_block: len(vertices)
            for building_block, vertices
            in self.building_block_vertices.items()
        }
        self.building_blocks = building_blocks
        self.topology_graph = topology_graph
