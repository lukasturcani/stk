import itertools
from collections import deque
from scipy.spatial.distance import euclidean
import numpy as np
import rdkit.Chem as chem 

from ..base import Topology
from ....convenience_tools import atom_vdw_radii
from ....pyWindow import window_sizes
from ...fg_info import FGInfo

class WindowError(Exception):
    def __init__(self, message):
        self.message = message

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

    def place_mols(self):
        """
        Places all building block molecules on correct coordinates.

        The building block molecules are placed in their appropriate 
        positions based on the topology. It does not join them. 

        
        Modifies
        --------
        self.macro_mol.mol
            An rdkit instance of the macromolecule with disconnected
            building blocks is placed in this attribute.
            
        Returns
        -------
        None : NoneType
        
        """
        
        self.macro_mol.mol = chem.Mol()
        
        # Get the StructUnit instances of the building blocks.
        bb1, bb2 = self.macro_mol.building_blocks
        # Get the number of functional groups in each building block.
        n_fg1 = len(bb1.functional_group_atoms())
        n_fg2 = len(bb2.functional_group_atoms())
        
        # Depending on the number of functional groups, assigned a 
        # building block to be either a linker or a building-block*. 
        if n_fg1 < n_fg2:
            lk = bb1
            n_lk = n_fg1
            bb = bb2
            n_bb = n_fg2
        else:
            lk = bb2
            n_lk = n_fg2
            bb = bb1
            n_bb = n_fg1
        
        # This loop places all building-blocks* on the points at 
        # `positions_A`. It then pairs all atoms which form a new bond
        # with the positions to which they will be bonding. It also
        # counts the nubmer of building-blocks* which make up the 
        # structure.
        for position in self.positions_A:
            self.macro_mol.mol = chem.CombineMols(
                                        self.macro_mol.mol, 
                                        position.place_mol(bb))
            # Update the counter each time a building-block* is added.
            self.bb_counter.update([bb])                            
            
            # Get ids of atoms which form new bonds.
            bonder_ids = deque(maxlen=n_bb)
            for atom in self.macro_mol.mol.GetAtoms():
                if (atom.HasProp('on_react') and 
                                atom.GetProp('on_react') == 'bond'):
                    bonder_ids.append(atom.GetIdx())
            
            # Save the ids of atoms which form new bonds and pair them
            # up with positions.
            position.bonder_ids = sorted(bonder_ids)
            self.pair_bonders_with_positions(position)

        # This loop places all linkers on the points at `positions_B`. 
        # It then saves all atoms which form a new bond to the position
        # they are found at. It also counts the number of linkers which 
        # make up the structure.
        for position in self.positions_B:
            self.macro_mol.mol = chem.CombineMols(
                                        self.macro_mol.mol, 
                                        position.place_mol(lk))
            # Update the counter each time a linker is added.
            self.bb_counter.update([lk])
            
            # Get ids of atoms which form new bonds.
            bonder_ids = deque(maxlen=n_lk)
            for atom in self.macro_mol.mol.GetAtoms():
                if (atom.HasProp('on_react') and 
                                atom.GetProp('on_react') == 'bond'):
                    bonder_ids.append(atom.GetIdx())
                    
            # Save the ids of atoms which form new bonds. 
            position.bonder_ids = list(bonder_ids)

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

        # If pyWindow failed, return ``None``.
        if all_windows is None:          
            return None
        
        all_windows = sorted(all_windows, reverse=True)[:self.n_windows]
        for i, x in enumerate(all_windows):
            # Return ``None`` when pyWindow fucks up and outputs a
            # mistakenly large window size.
            if x > 500:
                return None
                
        return all_windows
        
    def cavity_size(self):
        """
        Returns the diameter of the cage cavity.

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

    def window_difference(self):
        """
        The total difference in all window sizes.
        
        Every combination of windows is considered and all the size
        differences are summed and returned. Only differences between
        windows of the same type are considered.
        
        Consider a triangular-based prism cage topology. Such a cage 
        will have triangular windows and square windows. You only want 
        to compare the triangulars with other triangular windows and 
        squares only with other squares.
        
        Returns
        -------
        float
            The total difference of window size when considering every
            combination of windows of the same type.
            
        None : NoneType
            If not all windows were found.
            
                       
        Raises
        ------
        WindowError
            When the number of found windows is less than the number of 
            expected windows. Likely due to a collapsed cage.
            
        """
        

        if self.windows is None or len(self.windows) < self.n_windows:
            return None

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
    
    def __init__(self, macro_mol, random_placement=True):
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