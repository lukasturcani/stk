import numpy as np

from .base import VertexOnlyCageTopology
from ..base import Vertex

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