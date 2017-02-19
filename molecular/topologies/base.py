"""
Defines the base ``Topology`` type.

Extending MMEA: Adding new topologies
-------------------------------------
> Cages
To add a new cage topology a new class should be created, named
after the topology. This class should inhertic the ``CageTopology``
class. This will give access to various methods which are necessary
for dealing with any cage molecule. See the documenation of
``CageTopology`` for more details.

The new class will only need to have five class attributes added:
    1) a list called `vertices`
    2) a list called `edges`
    3) `n_windows` which holds the number of windows the cage
       topology has
    4) `n_window_types` which holds the number of different window
       types. For example, if `n_window_types` is 2 then the
       topology will have two kinds of windows, each with a
       different expected size even in a perfectly symmetrical case.
       Windows of the same type are expected to be of the same size.

The `vertices` list holds instances of the class ``Vertex``. Each
instance represents a vertex of a cage and needs to be initialized
with the coordinates of that vertex. Vertices of a cage are where
building-blocks* of cages are placed.

The `edges` list holds instances of the class ``Edge``. Each
instance represents an edge of a cage and needs to be initialized
with two instnaces of the ``Vertex`` class. The ``Vertex`` instances
should be held in the `vertices` list mentioned above. These are
the two vertices which the edge connects. Linkers of cages are
placed on edges. The edge instances automatically derive their
positions from the vertices supplied during initialization.

The vertices need to be positioned such that the center of the
topology is at the origin.


"""



from collections import Counter
import rdkit
import rdkit.Chem as chem

from ..fg_info import double_bond_combs

class Topology:
    """
    Builds macromolecules.

    More accurately, child classes of ``Topology`` class take care of
    building macromolecules.

    This class directly defines any operations and attributes that are
    needed by any topology during assembly. However, this class is not
    used directly by MMEA. It is intended to be inherited from. All
    macromolecules within MMEA have a `topology` attribute which hold
    an instance of a Topology child class. Child classes of Topology
    define operations specific to that one topology. For example, each
    child class must define a `join_mols()` method which creates bonds
    between the building blocks of a macromolecule. The way in which
    this is done will depend on what kind of macromolecules are being
    built. In addition, each child class must define methods which
    place the building blocks in approriate positions.

    """

    def build(self, macro_mol):
        """
        Assembles rdkit instances of macromolecules.

        Parameters
        ----------
        macro_mol : MacroMolecule
            The MacroMolecule instance which needs to be built.

        Modifies
        --------
        macro_mol.mol : rdkit.Chem.rdchem.Mol
            An rdkit instance of the assembled macromolecule is placed
            in this attribute.

        macro_mol.bb_counter : Counter
            The counter is updated with the amounts of building blocks
            used to make the macromolecule.

        macro_mol.bonds_made : int
            This count is updated with the number of bonds made during
            assembly.

        Returns
        -------
        None : NoneType

        """

        # When running ``build()`` in parallel, the atom tags are
        # cleared by the multiprocessing module. Make sure to readd the
        # tags before running ``build()``.
        for bb in macro_mol.building_blocks:
            bb.tag_atoms()

        self.place_mols(macro_mol)
        self.join_mols(macro_mol)
        self.del_atoms(macro_mol)

    def del_atoms(self, macro_mol):
        """
        Deletes the atoms which are lost during assembly.

        Modifies
        --------
        macro_mol.mol : rdkit.Chem.rdchem.Mol
            Atoms are removed from this molecule.

        Returns
        -------
        None : NoneType

        """

        mol = chem.EditableMol(macro_mol.mol)
        # Delete atoms with largest id last, so that when deleting
        # later atoms their ids do not change.
        for atom in reversed(macro_mol.mol.GetAtoms()):
            if atom.HasProp('del'):
                mol.RemoveAtom(atom.GetIdx())

        macro_mol.mol = mol.GetMol()

    @staticmethod
    def determine_bond_type(macro_mol, atom1_id, atom2_id):
        """
        Returns the bond order to be formed between the atoms.

        Some atoms will need to have a double bond created between
        them. This is defined in the `double_bond_combs` list. If the
        atom ids provided as paramters belong to functional grups found
        in this list, the rdkit double bond type will be returned.
        If not the rdkit single bond type will be returned. These types
        are needed when adding bonds using ``EditableMol`` instances.

        Parameters
        ----------
        macro_mol : MacroMolecule
            The macromolecule being assembled.

        atom1_id : int
            The id number of the first atom.

        atom2_id : int
            The id number of the second atom.

        Returns
        -------
        rdkit.Chem.rdchem.BondType.SINGLE
            If the atoms don't belong to functional groups which form
            a double bond.

        rdkit.Chem.rdchem.BondType.DOUBLE
            If the atoms belong to functional groups which form a
            double bond.

        """

        # Get the functional groups of the of the atoms whose atom ids
        # were supplied as arguments. If the groups form a tuple in
        # `double_bond_combs` return a rdkit double bond type. If they
        # do not, return a rdkit single bond type.

        atom1 = macro_mol.mol.GetAtomWithIdx(atom1_id)
        atom1_grp = atom1.GetProp('fg')
        atom2 = macro_mol.mol.GetAtomWithIdx(atom2_id)
        atom2_grp= atom2.GetProp('fg')

        double_bond_present = ((atom1_grp, atom2_grp) == tup or
                               (atom2_grp, atom1_grp) == tup for
                               tup in double_bond_combs)

        if any(double_bond_present):
            return rdkit.Chem.rdchem.BondType.DOUBLE
        else:
            return rdkit.Chem.rdchem.BondType.SINGLE
