from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import rdkit.Chem.AllChem as rdkit

from stk._internal.atom import Atom
from stk._internal.bonds import DativeBond, IntegerBond
from stk._internal.functional_group_factories.functional_group_factory import (
    FunctionalGroupFactory,
)
from stk._internal.functional_groups.functional_group import FunctionalGroup
from stk._internal.math import (
    get_acute_vector,
    get_centroid,
    get_orthogonal_vector,
    get_plane_normal,
)
from stk._internal.utilities.utilities import remake


@dataclass(frozen=True, slots=True)
class RotationAnchor:
    axis: npt.NDArray[np.float32]
    target: npt.NDArray[np.float32] | None


@dataclass(frozen=True, slots=True)
class BuildingBlock:
    """
    Represents a building block of a :class:`.ConstructedMolecule`.

    A :class:`BuildingBlock` can represent either an entire molecule or
    a molecular fragments used to construct a
    :class:`.ConstructedMolecule`. The building block uses
    :class:`.FunctionalGroup` instances to identify which atoms are
    modified during construction.

    See Also:
        * :class:`.Atom`: Represents atoms of a building block.
        * :class:`.Bond`: Represents bonds of a building block.
        * :class:`.FunctionalGroup`: Represents functional
          groups of a building block.
    """

    atoms: Sequence[Atom]
    integer_bonds: Sequence[IntegerBond]
    dative_bonds: Sequence[DativeBond]
    functional_groups: Sequence[FunctionalGroup]
    position_matrix: npt.NDArray[np.float32]
    position_anchor: npt.NDArray[np.float32]
    rotation_anchor: RotationAnchor | None

    @staticmethod
    def from_smiles(
        smiles: str,
        functional_groups: FunctionalGroup
        | FunctionalGroupFactory
        | Iterable[FunctionalGroup | FunctionalGroupFactory] = (),
    ) -> "BuildingBlock":
        molecule = rdkit.AddHs(rdkit.MolFromSmiles(smiles))
        rdkit.Kekulize(molecule)
        return BuildingBlock.from_rdkit(molecule, functional_groups)

    @staticmethod
    def from_rdkit(
        molecule: rdkit.Mol,
        functional_groups: FunctionalGroup
        | FunctionalGroupFactory
        | Iterable[FunctionalGroup | FunctionalGroupFactory] = (),
    ) -> "BuildingBlock":
        if isinstance(
            functional_groups,
            FunctionalGroup | FunctionalGroupFactory,
        ):
            functional_groups = (functional_groups,)

        molecule = remake(molecule)
        atoms = []
        for atom in molecule.GetAtoms():
            atoms.append(
                Atom(
                    id=atom.GetIdx(),
                    atomic_number=atom.GetAtomicNum(),
                    charge=atom.GetFormalCharge(),
                )
            )
        integer_bonds = []
        dative_bonds = []
        for bond in molecule.GetBonds():
            atom1 = atoms[bond.GetBeginAtomIdx()]
            atom2 = atoms[bond.GetEndAtomIdx()]
            match bond.GetBondType():
                case (
                    rdkit.BondType.SINGLE
                    | rdkit.BondType.DOUBLE
                    | rdkit.BondType.TRIPLE
                    | rdkit.BondType.QUADRUPLE
                    | rdkit.Bondype.QUINTUPLE
                    | rdkit.BondType.HEXTUPLE
                ):
                    integer_bonds.append(
                        IntegerBond(
                            atom1=atom1,
                            atom2=atom2,
                            order=int(bond.GetBondTypeAsDouble()),
                        )
                    )
                case rdkit.BondType.DATIVE:
                    dative_bonds.append(DativeBond(atom1, atom2, 1))

        normalized_fgs = []
        for functional_group in functional_groups:
            match functional_group:
                case FunctionalGroup():
                    normalized_fgs.append(functional_group)
                case FunctionalGroupFactory():
                    normalized_fgs.extend(
                        functional_group.get_functional_groups(molecule)
                    )

        position_matrix = np.array(
            molecule.GetConformer().GetPositions(),
            dtype=np.float32,
        )
        centroid = get_centroid(position_matrix)

        rotation_axis = None
        if len(normalized_fgs) == 2:
            centroid1, centroid2 = fg_centroids = [
                get_centroid(position_matrix[fg.bonders])
                for fg in normalized_fgs
            ]
            fg_centroid = get_centroid(fg_centroids)
            rotation_axis = get_acute_vector(
                reference=centroid - fg_centroid,
                vector=get_orthogonal_vector(centroid2 - centroid1),
            )

        if len(normalized_fgs) >= 3:
            fg_centroids = [
                get_centroid(position_matrix[fg.bonders])
                for fg in normalized_fgs
            ]
            rotation_axis = get_plane_normal(fg_centroids)
            fg_centroid = get_centroid(fg_centroids)
            fg_centroid_to_centroid = centroid - fg_centroid
            rotation_axis = get_acute_vector(
                reference=fg_centroid_to_centroid,
                vector=rotation_axis,
            )

        rotation_target = None
        if normalized_fgs and rotation_axis is not None:
            target_fg = normalized_fgs[0]
            rotation_target = get_centroid(position_matrix[target_fg.bonders])

        return BuildingBlock(
            atoms=atoms,
            integer_bonds=integer_bonds,
            dative_bonds=dative_bonds,
            functional_groups=normalized_fgs,
            position_matrix=position_matrix,
            position_anchor=centroid,
            rotation_anchor=RotationAnchor(rotation_axis, rotation_target)
            if rotation_axis is not None
            else None,
        )
