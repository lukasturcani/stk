from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import rdkit.Chem.AllChem as rdkit

from stk._internal.atom import Atom
from stk._internal.bonds import DativeBond, IntegerBond
from stk._internal.functional_group import FunctionalGroup
from stk._internal.functional_group_factory import FunctionalGroupFactory
from stk._internal.math import (
    get_acute_vector,
    get_centroid,
    get_orthogonal_vector,
    get_plane_normal,
)
from stk._internal.utilities.utilities import remake


class EmbedError(Exception):
    pass


class NumFunctionalGroupsError(Exception):
    pass


BuildingBlock = ""


def to_rdkit(
    atoms: Sequence[Atom],
    integer_bonds: Sequence[IntegerBond],
    dative_bonds: Sequence[DativeBond],
    position_matrix: npt.NDArray[np.float32],
) -> rdkit.Mol:
    mol = rdkit.EditableMol(rdkit.Mol())
    for atom in atoms:
        rdkit_atom = rdkit.Atom(atom.atomic_number)
        rdkit_atom.SetFormalCharge(atom.charge)
        mol.AddAtom(rdkit_atom)

    for ibond in integer_bonds:
        mol.AddBond(
            beginAtomIdx=ibond.atom1.id,
            endAtomIdx=ibond.atom2.id,
            order=rdkit.BondType(ibond.order),
        )
    for dbond in dative_bonds:
        mol.AddBond(dbond.atom1.id, dbond.atom2.id, rdkit.BondType.DATIVE)

    mol = mol.GetMol()
    conformer = rdkit.Conformer(len(atoms))
    for atom_id, atom_coord in enumerate(position_matrix):
        conformer.SetAtomPosition(atom_id, atom_coord.astype(np.float64))
        mol.GetAtomWithIdx(atom_id).SetNoImplicit(True)
    mol.AddConformer(conformer)
    return mol


@dataclass(frozen=True, slots=True)
class Molecule:
    atoms: Sequence[Atom]
    integer_bonds: Sequence[IntegerBond]
    dative_bonds: Sequence[DativeBond]

    @staticmethod
    def from_rdkit(molecule: rdkit.Mol) -> "Molecule":
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
        return Molecule(atoms, integer_bonds, dative_bonds)


@dataclass(frozen=True, slots=True)
class RotationBuildingBlock:
    atoms: Sequence[Atom]
    integer_bonds: Sequence[IntegerBond]
    dative_bonds: Sequence[DativeBond]
    functional_groups: Sequence[FunctionalGroup]
    position_matrix: npt.NDArray[np.float32]
    position_anchor: npt.NDArray[np.float32]
    rotation_anchor_axis: npt.NDArray[np.float32]
    rotation_anchor_target: npt.NDArray[np.float32]

    @staticmethod
    def from_smiles(
        smiles: str,
        functional_groups: FunctionalGroup
        | FunctionalGroupFactory
        | Iterable[FunctionalGroup | FunctionalGroupFactory] = (),
        position_matrix: npt.NDArray[np.float32] | None = None,
    ) -> "RotationBuildingBlock":
        molecule = rdkit.AddHs(rdkit.MolFromSmiles(smiles))
        if position_matrix is None:
            params = rdkit.ETKDGv3()
            params.randomSeed = 4
            if rdkit.EmbedMolecule(molecule, params) == -1:
                raise EmbedError("failed to embed building block")
        rdkit.Kekulize(molecule)
        return RotationBuildingBlock.from_rdkit(molecule, functional_groups)

    @staticmethod
    def from_rdkit(
        molecule: rdkit.Mol,
        functional_groups: FunctionalGroup
        | FunctionalGroupFactory
        | Iterable[FunctionalGroup | FunctionalGroupFactory] = (),
    ) -> "RotationBuildingBlock":
        if isinstance(
            functional_groups,
            FunctionalGroup | FunctionalGroupFactory,
        ):
            functional_groups = (functional_groups,)

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

        match len(normalized_fgs):
            case 0 | 1:
                raise NumFunctionalGroupsError(
                    "needs 2 or more functional groups"
                )

            case 2:
                centroid1, centroid2 = fg_centroids = [
                    get_centroid(position_matrix[fg.bonders])
                    for fg in normalized_fgs
                ]
                fg_centroid = get_centroid(fg_centroids)
                rotation_axis = get_acute_vector(
                    reference=centroid - fg_centroid,
                    vector=get_orthogonal_vector(centroid2 - centroid1),
                )
            case _:
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

        target_fg = normalized_fgs[0]
        rotation_target = get_centroid(position_matrix[target_fg.bonders])

        stk_molecule = Molecule.from_rdkit(molecule)
        return RotationBuildingBlock(
            atoms=stk_molecule.atoms,
            integer_bonds=stk_molecule.integer_bonds,
            dative_bonds=stk_molecule.dative_bonds,
            functional_groups=normalized_fgs,
            position_matrix=position_matrix,
            position_anchor=centroid,
            rotation_anchor_axis=rotation_axis,
            rotation_anchor_target=rotation_target,
        )

    def to_rdkit(self) -> rdkit.Mol:
        mol = rdkit.EditableMol(rdkit.Mol())
        for atom in self.atoms:
            rdkit_atom = rdkit.Atom(atom.atomic_number)
            rdkit_atom.SetFormalCharge(atom.charge)
            mol.AddAtom(rdkit_atom)

        for ibond in self.integer_bonds:
            mol.AddBond(
                beginAtomIdx=ibond.atom1.id,
                endAtomIdx=ibond.atom2.id,
                order=rdkit.BondType(ibond.order),
            )
        for dbond in self.dative_bonds:
            mol.AddBond(dbond.atom1.id, dbond.atom2.id, rdkit.BondType.DATIVE)

        mol = mol.GetMol()
        conformer = rdkit.Conformer(len(self.atoms))
        for atom_id, atom_coord in enumerate(self.position_matrix):
            conformer.SetAtomPosition(atom_id, atom_coord.astype(np.float64))
            mol.GetAtomWithIdx(atom_id).SetNoImplicit(True)
        mol.AddConformer(conformer)
        return mol
