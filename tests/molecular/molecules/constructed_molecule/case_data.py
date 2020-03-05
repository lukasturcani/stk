class CaseData:
    def __init__(
        self,
        constructed_molecule,
        num_new_atoms,
        num_new_bonds,
        building_block_counts,
        building_block_vertices,
        building_blocks,
        topology_graph,
    ):
        self.constructed_molecule = constructed_molecule
        self.num_new_atoms = num_new_atoms
        self.num_new_bonds = num_new_bonds
        self.building_block_counts = building_block_counts
        self.building_block_vertices = building_block_vertices
        self.building_blocks = building_blocks
        self.topology_graph = topology_graph
