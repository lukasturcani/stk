import itertools
from collections import deque
from scipy.spatial.distance import euclidean
import numpy as np
import rdkit.Chem as chem 

from ..base import Topology
from ....convenience_functions import LazyAttr, atom_vdw_radii
from ....pyWindow import window_sizes
from ...molecular import FGInfo

class CageTopology(Topology):
    """
    A topology class which cage topologies should inherit.
        
    Attributes
    ----------
    In addition to all the attributes defined within ``Topology`` this
    class has the following attributes:

    pair_up : function object (default = pair_up_edges_with_vertices)
        This is the function which pairs up molecules placed using the
        ``Vertex`` and ``Edge`` classes. This should be how cage
        topologies should be defined.    
    
    """

    @property
    def windows(self):
        """
        
        Returns
        -------
        None : NoneType
            If the function for finding windows and their sizes
            found fewer than the required number of windows or
            if it failed for some other reason.
            
        list of floats
            Each float in the list represents the size of a window in
            the cage. If the window finding function found more than
            the expected number of windows, only the largest n windows
            are returned. Where n is the number of expected windows.
        
        """
        
        all_windows = window_sizes(self.macro_mol.prist_mol_file)

        if all_windows is None:          
            return None
        
        all_windows = sorted(all_windows, reverse=True)[:self.n_windows]
        for i, x in enumerate(all_windows):
            if x > 500:
                all_windows[i] = 500
                
        return all_windows
   
    def cavity_size(self):
        """
        Returns the size of the cage cavity.

        Returns
        -------
        float
            The size of the cage cavity.        
        
        """
        
        center_of_mass = self.macro_mol.center_of_mass('prist')
        min_dist = min((euclidean(coord, center_of_mass) -
        atom_vdw_radii[self.macro_mol.atom_symbol('prist', atom_id)]) 
                                 for atom_id, coord in 
                                 self.macro_mol.all_atom_coords('prist'))
        return 2 * abs(min_dist)    

    def window_difference(self, default=500):
        """
        The total difference between all cage size.
        
        Every combination of windows is considered and all the size
        differences are summed and returned.

        Parameters
        ----------
        default : float or int (default = 500)
            The number returned if no windows were found in the cage.
        
        Returns
        -------
        default : float or int
            If the `windows` attribute is ``None``. This happens when
            the window finidning algorithm fails. In these cases the
            `default` value is returned.
            
        float
            The total difference of window size when considering every
            combination of windows.
                           
        """
        
        if self.windows is None:
            return default
        if len(self.windows) < self.n_windows:
            return default
        if 500 in self.windows:
            return default
    
        # Cluster the windows into groups so that only size differences
        # between windows of the same type are taken into account. To do
        # this, first sort the windows by size. If two windows types are
        # present split the windows at the two groups at the point where
        # the window sizes have the biggest difference. If there are
        # three types split it at the two biggest differences and so on.
                
        windows = np.array(self.windows)
        
        diffs = list(abs(np.ediff1d(windows)))
        sorted_diffs = sorted(diffs, reverse=True)
        
        # Get indices of where the list should be split.
        split = []
        for x in range(self.n_window_types-1):
            i = diffs.index(sorted_diffs[x]) + 1
            split.append(i)
    
        # Get the sub-lists.
        og = list(windows)
        clusters = []
        for i in sorted(split, reverse=True):
            clusters.append(og[i:])            
            og = og[:i]

        if self.n_window_types == 1:
            clusters.append(og)

        # After this sum the differences in each group and then sum the
        # group totals.
        diff_sums = []
        for cluster in clusters:
            diff_sum = sum(abs(w1 - w2) for w1, w2 in 
                                    itertools.combinations(cluster, 2))

            diff_num = sum(1 for _ in 
                itertools.combinations(cluster, 2))
            
            diff_sums.append(diff_sum / diff_num)
            
        return sum(diff_sums)
        
class VertexOnlyCageTopology(CageTopology): 
    
    def __init__(self, macro_mol, random_placement=False):
        Topology.__init__(self, macro_mol)
        self.random_placement = random_placement
        self.connect()
        
    def place_mols(self):
        
        self.macro_mol.heavy_mol = chem.Mol()        
        
        if self.random_placement:
            return self.place_mols_random()
        return self.place_mols_assigned()
        
    def place_mols_random(self):
        for position in self.positions_A:
            bb = np.random.choice(self.macro_mol.building_blocks)
            n_bb = len(bb.find_functional_group_atoms())
            
            self.macro_mol.heavy_mol = chem.CombineMols(
                                        self.macro_mol.heavy_mol,
                                        position.place_mol(bb))
            self.bb_counter.update([bb])                                        
                                        
            heavy_ids = deque(maxlen=n_bb)
            for atom in self.macro_mol.heavy_mol.GetAtoms():
                if atom.GetAtomicNum() in FGInfo.heavy_atomic_nums:
                    heavy_ids.append(atom.GetIdx())
            
            position.heavy_ids = sorted(heavy_ids)
            self.pair_heavy_ids_with_connected(position)

    @classmethod
    def connect(cls):
        if getattr(cls, 'connected', False):
            return
        
        for v1, v2 in cls.connections:
            v1.connected.append(v2)
            v2.connected.append(v1)                                    
        cls.connected = True
        
    def place_mols_assigned(self):
        pass