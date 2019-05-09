"""
Defines classes for cage molecules.

The difference between :class:`CageStructUnit` and :class:`Cage`
is that :class:`CageStructUnit` is for cages loaded directly into
``stk``, while :class:`Cage` is for cages constructed by ``stk``.

"""

import numpy as np
import itertools as it
import warnings
import rdkit.Chem.AllChem as rdkit
import pywindow
from scipy.optimize import minimize
from sklearn.metrics.pairwise import euclidean_distances

from .macro_molecule import MacroMolecule
from .struct_unit import StructUnit
from ...utilities import atom_vdw_radii


class CageBase:
    """
    A base class for cage molecules.

    """

    def _cavity_size(self, origin, conformer):
        """
        Calculates diameter of the molecule from `origin`.

        The cavity is measured by finding the atom nearest to
        `origin`, correcting for van der Waals diameter and multiplying
        by -2.

        This function should not be used directly. Use
        :meth:`cavity_size` instead, which finds the optimal value of
        `origin` to use.

        Parameters
        ----------
        origin : :class:`numpy.ndarray`
            Holds the x, y and z coordinate of the position from which
            the cavity is measured.

        conformer : :class:`int`
            The id of the conformer to be used.

        Returns
        -------
        :class:`float`
            The (negative) diameter of the molecules cavity.

        """

        atom_vdw = np.array([atom_vdw_radii[x.GetSymbol()] for x
                            in self.mol.GetAtoms()])
        pos_mat = self.mol.GetConformer(conformer).GetPositions()
        distances = euclidean_distances(pos_mat, np.array([origin]))
        distances = distances.flatten() - atom_vdw
        return -2*min(distances)

    def cavity_size(self, conformer=-1):
        """
        Calculates the diameter of the molecule's cavity.

        Parameters
        ----------
        conformer : :class:`int`, optional
            The id of the conformer to be used.

        Returns
        -------
        :class:`float`
            The diameter of the molecule's cavity in Angstroms.

        """

        # This function uses _cavity_size() to calculate the cavity
        # size. _cavity_size() finds the closest atom to `origin` to
        # get its value of the cavity.

        # What this function does is finds the value of `origin` which
        # causes _cavity_size() to calculate the largest possible
        # cavity.
        ref = self.center_of_mass(conformer)
        icavity = 0.5*self._cavity_size(ref, conformer)
        bounds = [(coord+icavity, coord-icavity) for coord in ref]
        cavity_origin = minimize(
                            lambda x: self._cavity_size(x, conformer),
                            x0=ref,
                            bounds=bounds).x
        cavity = -self._cavity_size(cavity_origin, conformer)
        return 0 if cavity < 0 else cavity


class CageStructUnit(CageBase, StructUnit):
    """
    Represents molecular cages loaded into ``stk``.

    """


class Cage(CageBase, MacroMolecule):
    """
    Represents molecular cages constructed by ``stk``.

    """

    def window_difference(self, conformer=-1):
        """
        The average difference across window sizes.

        First take the average window difference between all windows
        of the same type. Then take the average of window differences
        across window types.

        Consider a triangular-based prism cage topology. Such a
        cage will have triangular windows and square windows. You
        only want to compare the triangulars with other
        triangular windows and squares only with other squares.

        Parameters
        ---------
        conformer : :class:`int`, optional
            The id of the conformer to use.

        Returns
        -------
        :class:`float`
            The total difference of window size when considering
            every combination of windows of the same type.

        None : :class:`NoneType`
            If not all windows were found.

        """

        windows = self.windows(conformer)

        if windows is None or len(windows) < self.topology.n_windows:
            return None

        windows = np.array(windows)

        # Cluster the windows into groups so that only size
        # differences between windows of the same type are taken
        # into account. To do this, first sort the windows by
        # size. If two windows types are present split the
        # windows at the two groups at the point where the window
        # sizes have the biggest difference. If there are three
        # types split it at the two biggest differences and so
        # on.

        diffs = list(abs(np.ediff1d(windows)))
        sorted_diffs = sorted(diffs, reverse=True)

        # Get indices of where the list should be split.
        split = []
        for x in range(self.topology.n_window_types-1):
            i = diffs.index(sorted_diffs[x]) + 1
            split.append(i)

        # Get the sub-lists.
        og = list(windows)
        clusters = []
        for i in sorted(split, reverse=True):
            clusters.append(og[i:])
            og = og[:i]

        if self.topology.n_window_types == 1:
            clusters.append(og)

        # After this sum the differences in each group and then
        # sum the group totals.
        diff_sums = []
        for cluster in clusters:
            diff_sum = sum(abs(w1 - w2) for w1, w2 in
                           it.combinations(cluster, 2))

            diff_num = sum(1 for _ in it.combinations(cluster, 2))

            diff_sums.append(diff_sum / diff_num)

        return np.mean(diff_sums)

    def window_variance(self, conformer=-1):
        """
        The variance in window sizes.

        For cages where multiple window types are present, the
        window variance within each window type is calculated first and
        the returned value is the mean of these variances.

        Parameters
        ---------
        conformer : :class:`int`, optional
            The id of the conformer to use.

        Returns
        -------
        :class:`float`
            The variance in window sizes.

        """

        windows = self.windows(conformer)

        if windows is None or len(windows) < self.topology.n_windows:
            return None

        windows = np.array(windows)

        # Cluster the windows into groups so that only size
        # differences between windows of the same type are taken
        # into account. To do this, first sort the windows by
        # size. If two windows types are present split the
        # windows at the two groups at the point where the window
        # sizes have the biggest difference. If there are three
        # types split it at the two biggest differences and so
        # on.

        diffs = list(abs(np.ediff1d(windows)))
        sorted_diffs = sorted(diffs, reverse=True)

        # Get indices of where the list should be split.
        split = []
        for x in range(self.topology.n_window_types-1):
            i = diffs.index(sorted_diffs[x]) + 1
            split.append(i)

        # Get the sub-lists.
        og = list(windows)
        clusters = []
        for i in sorted(split, reverse=True):
            clusters.append(og[i:])
            og = og[:i]

        if self.topology.n_window_types == 1:
            clusters.append(og)

        variances = [np.var(c) for c in clusters]
        return np.mean(variances)

    def windows(self, conformer=-1):
        """
        Returns window sizes found by ``pywindow``.

        Parameters
        ----------
        conformer : :class:`int`, optional
            The id of the conformer to use.

        Returns
        -------
        None : :class:`NoneType`
            If the function for finding windows and their sizes
            found fewer than the required number of windows or
            if it failed for some other reason.

        :class:`list` of :class:`float`
            Each :class:`float` represents the size of a
            window in the cage. If the window finding function
            found more than the expected number of windows, only
            the largest ``n`` windows are returned. Where ``n`` is the
            number of expected windows for that cage topology.

        """

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # Load an RDKit molecule object to pywindow.

            # As pywindow doesnt support multiple conformers, first
            # make an rdkit molecule holding only the desired
            # conformer.
            new_mol = rdkit.Mol(self.mol)
            new_mol.RemoveAllConformers()
            new_mol.AddConformer(self.mol.GetConformer(conformer))

            loader = pywindow.molecular.Molecule.load_rdkit_mol
            pw_molecule = loader(new_mol)
            # Find windows and get a single array with windows' sizes.
            all_windows = pw_molecule.calculate_windows()

        # If pywindow failed, return ``None``.
        if all_windows is None:
            return None

        # pywindow sometimes detects super large windows by accident,
        # filter them out first.
        valid_windows = sorted((w for w in all_windows if w < 1e6),
                               reverse=True)
        return valid_windows[:self.topology.n_windows]


class CageComplex(MacroMolecule):
    """
    Represents a cage-guest system.

    """

    pass
