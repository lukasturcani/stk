from .base import CageTopology
from ..base import Vertex

class SixPlusEight(CageTopology):
    
    positions_A = [Vertex(-50, 50, 0), 
                    Vertex(-50, -50, 0), 
                    Vertex(50, 50, 0), 
                    Vertex(50, -50, 0),
    
                    Vertex(0, 0, 50), 
                    Vertex(0, 0, -50)]

    a,b,c,d,e,f = positions_A

    positions_B = [Vertex.vertex_init(a,e,b),
                 Vertex.vertex_init(b,e,d),
                Vertex.vertex_init(e,d,c),
                Vertex.vertex_init(e,c,a),
    
                Vertex.vertex_init(a,f,b),
                Vertex.vertex_init(f,b,d),
                Vertex.vertex_init(d,f,c),
                Vertex.vertex_init(c,f,a)]
                
    n_windows = 12
    n_window_types = 1