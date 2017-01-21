import numpy as np
import itertools

from .base import CageTopology, Vertex, Edge

class TwoPlusThree(CageTopology):
    positions_A = [Vertex(0,0,20), Vertex(0,0,-20)]
    a,b = positions_A
    
    positions_B = [Edge(a, b),
                  Edge(a, b),
                  Edge(a,b)]
    
    alpha, beta, gamma = positions_B
  
    b.edge_plane_normal = lambda a=a: np.multiply(a.edge_plane_normal(), 
                                                  -1)
         
    alpha.coord = np.array([-20,
                              -10*np.sqrt(3),
                                0])

    beta.coord = np.array([20,
                              -10*np.sqrt(3),
                                0])

    gamma.coord = np.array([0,
                              10*np.sqrt(3),
                                0])
                                
    n_windows = 3
    n_window_types = 1

class FourPlusSix(CageTopology):
    """
    Defines the tetrahedral, 4+6, topology.

    This is a topology of cages where 4 building-blocks* are placed on
    vertices and 6 linkers are placed on the edges between them. This
    class defines functions which place these molecules in the correct
    positions within an rdkit instance. The rdkit instance is stored in 
    the `heavy_mol` attribute of a ``Cage`` instance.
        
    """
    
    
    # Vertices of a tetrahdron so that origin is at the origin. Source:
    # http://tinyurl.com/lc262h8.
    positions_A = [Vertex(100,0,-100/np.sqrt(2)), 
                Vertex(-100,0,-100/np.sqrt(2)), 
                Vertex(0,100,100/np.sqrt(2)), 
                Vertex(0,-100,100/np.sqrt(2))]
        
    positions_B = [Edge(v1,v2) for v1, v2 in 
                itertools.combinations(positions_A, 2)]  
    
    n_windows = 4
    n_window_types = 1

class SixPlusNine(CageTopology):
    
    # source: http://eusebeia.dyndns.org/4d/prism3
    positions_A = [Vertex(-50,-50/np.sqrt(3),-50), 
                Vertex(-50,-50/np.sqrt(3),50),
                Vertex(50,-50/np.sqrt(3),-50),
                Vertex(50,-50/np.sqrt(3),50),
                Vertex(0, 100/np.sqrt(3),-50),
                Vertex(0, 100/np.sqrt(3), 50)]
           
    a,b,c, d,e,f = positions_A          
           
    positions_B = [Edge(a, b),
             Edge(a, c),
             Edge(c, d),
             Edge(b, d),
             Edge(a, e),
             Edge(c, e),
             Edge(b, f),
             Edge(d, f),
             Edge(e, f)]
    
    n_windows = 5
    n_window_types = 1

class EightPlusTwelve(CageTopology):
    """
    Defines a square-like topology.    
    
    """
    
    positions_A = [Vertex(-50, 50, -50), 
                Vertex(-50, -50, -50), 
                Vertex(50, 50, -50), 
                Vertex(50, -50, -50),

                Vertex(-50, 50, 50), 
                Vertex(-50, -50, 50), 
                Vertex(50, 50, 50), 
                Vertex(50, -50, 50)]
    
    a,b,c,d, e,f,g,h = positions_A
    
    positions_B = [Edge(a, c), 
             Edge(a, b),
             Edge(b, d),
             Edge(c, d),
             
             Edge(e, g), 
             Edge(e, f),
             Edge(f, h),
             Edge(g, h),


             Edge(a, e), 
             Edge(b, f),
             Edge(c, g),
             Edge(d, h)]  
    
    n_windows = 6  
    n_window_types = 1

class Dodecahedron(CageTopology):
    
    # Source: http://tinyurl.com/h2dl949
    phi = (1 + np.sqrt(5))/2
    x = 50
    positions_A = [Vertex(x* phi, 0.0, x/phi), 
                Vertex(x*-phi, 0.0, x/phi), 
                Vertex(x*-phi, 0.0, x/-phi),
                Vertex(x* phi, 0.0, x/-phi), 
                
                Vertex(x/ phi, x* phi, 0.0), 
                Vertex(x/ phi, x*-phi, 0.0),
                Vertex(x/-phi, x*-phi, 0.0), 
                Vertex(x/-phi, x* phi, 0.0), 
                Vertex(0.0, x/ phi, x* phi),
                Vertex(0.0, x/ phi, x*-phi), 
                Vertex(0.0, x/-phi, x*-phi), 
                Vertex(0.0, x/-phi, x* phi),
  
                Vertex( x, x, x), 
                Vertex( x,-x, x), 
                Vertex(-x,-x, x), 
                Vertex(-x, x, x), 
                Vertex(-x, x,-x), 
                Vertex( x, x,-x),
                Vertex( x,-x,-x), 
                Vertex(-x,-x,-x)]

    A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T = positions_A
    positions_B = [Edge(A,N), Edge(A,M), Edge(A,D), Edge(B,O), Edge(B,P),
             Edge(B,C), Edge(C,T), Edge(C,Q), Edge(D,S), Edge(D,R),
             Edge(E,M), Edge(E,H), Edge(E,R), Edge(F,G), Edge(F,S),
             Edge(F,N), Edge(G,O), Edge(G,T), Edge(H,P), Edge(H,Q),
             Edge(I,L), Edge(I,M), Edge(I,P), Edge(J,K), Edge(J,R),
             Edge(J,Q), Edge(K,S), Edge(K,T), Edge(L,O), Edge(L,N)]    
    
    n_windows = 12
    n_window_types = 1
