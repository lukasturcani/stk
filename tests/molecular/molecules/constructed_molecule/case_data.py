import stk


class CaseData:
    def __init__(
        self,
        topology_graph,
        num_new_atoms,
        num_new_bonds,
    ):
        self.constructed_molecule = stk.ConstructedMolecule(
            topology_graph=topology_graph,
        )
        self.num_new_atoms = num_new_atoms
        self.num_new_bonds = num_new_bonds
        self.topology_graph = topology_graph
