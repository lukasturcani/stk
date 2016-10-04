import numpy as np
import rdkit.Chem as chem

from .base import VertexOnlyCageTopology
from ..base import Vertex

class OnePlusOne(VertexOnlyCageTopology):
    positions_A = [Vertex(25,0,0),
                   Vertex(-25,0,0)]
    a,b = positions_A
    connections = [(a,b)]
    
    a.edge_plane_normal = lambda : np.array([1,0,0])
    b.edge_plane_normal = lambda : np.array([-1,0,0])
    
    a.edge_centroid = lambda : np.array([0,0,0])
    b.edge_centroid = lambda : np.array([0,0,0])
    
    n_windows = 3
    n_window_types = 1    
    
    def join_mols(self):
        
        editable_mol = chem.EditableMol(self.macro_mol.heavy_mol)
        
        for position in self.positions_A:
            other_position = next(x for x in self.positions_A if 
                                x is not position)
            
            position.atom_position_pairs = [(atom, other_position) for 
                                    atom in position.heavy_ids]
            
            
            for atom_id, vertex in position.atom_position_pairs:
                # Get all the distances between the atom and the heavy
                # atoms on the vertex. Store this information on the 
                # vertex.
                for atom2_id in vertex.heavy_ids:
                    distance = self.macro_mol.atom_distance('heavy', 
                                                            atom_id, 
                                                            atom2_id)
                    position.distances.append((distance, 
                                             atom_id, atom2_id))

        
        paired = set()        
        for position in self.positions_A:
            for _, atom1_id, atom2_id in sorted(position.distances):
                if atom1_id in paired or atom2_id in paired:
                    continue            

                bond_type = self.determine_bond_type(atom1_id, atom2_id)
                # Add the bond.                
                editable_mol.AddBond(atom1_id, atom2_id, bond_type)
                self.bonds_made += 1
                paired.add(atom1_id)
                paired.add(atom2_id)
                
        self.macro_mol.heavy_mol = editable_mol.GetMol()     

class TwoPlusTwo(VertexOnlyCageTopology):
    positions_A = [Vertex(50,0,-50/np.sqrt(2)), 
                Vertex(-50,0,-50/np.sqrt(2)), 
                Vertex(0,50,50/np.sqrt(2)), 
                Vertex(0,-50,50/np.sqrt(2))]
                
    a,b,c,d = positions_A

    for x in positions_A:
        old_normal = x.edge_plane_normal
        x.edge_plane_normal = lambda a=old_normal: np.multiply(a(), -1)

    
    connections = [(a,b), (a,c), (a,d),
                   (b,c), (b,d),
                   (c,d)]
                   
    n_windows = 4
    n_window_types = 1
                   
class FourPlusFour(VertexOnlyCageTopology):
    positions_A = [Vertex(-50, 50, -50), 
                Vertex(-50, -50, -50), 
                Vertex(50, 50, -50), 
                Vertex(50, -50, -50),

                Vertex(-50, 50, 50), 
                Vertex(-50, -50, 50), 
                Vertex(50, 50, 50), 
                Vertex(50, -50, 50)]
                
    a,b,c,d, e,f,g,h = positions_A
    
    connections = [(a,b), (a,c), (a,e), (b,d), (b,f), (c,g), (c,d),
                   (d,h), (e,g), (e,f), (f,h), (g,h)]
    
    n_windows = 6
    n_window_types = 1
    
    
    
    