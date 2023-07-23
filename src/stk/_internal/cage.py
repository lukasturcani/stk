import itertools
from collections import defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, replace
from typing import TypeAlias

import numpy as np
import numpy.typing as npt
import rdkit.Chem.AllChem as rdkit

from stk._internal.atom import Atom
from stk._internal.bonds import DativeBond, IntegerBond
from stk._internal.building_block import RotationBuildingBlock, to_rdkit
from stk._internal.math import normalize_vector
from stk._internal.vertices import RotationVertex


@dataclass(frozen=True, slots=True)
class Cage:
    atoms: Sequence[Atom]
    integer_bonds: Sequence[IntegerBond]
    dative_bonds: Sequence[DativeBond]
    position_matrix: npt.NDArray[np.float32]

    def to_rdkit(self) -> rdkit.Mol:
        return to_rdkit(
            atoms=self.atoms,
            integer_bonds=self.integer_bonds,
            dative_bonds=self.dative_bonds,
            position_matrix=self.position_matrix,
        )


Assignment: TypeAlias = tuple[RotationVertex, RotationBuildingBlock]


@dataclass(frozen=True, slots=True)
class TwoGraph:
    primary_vertices: Sequence[RotationVertex]
    secondary_vertices: Sequence[RotationVertex]

    @staticmethod
    def from_primary_vertex_positions(
        positions: Sequence[Sequence[float]],
        primary_alignments: dict[int, int] | None = None,
        secondary_alignments: dict[int, int] | None = None,
    ) -> "TwoGraph":
        if primary_alignments is None:
            primary_alignments = {}
        if secondary_alignments is None:
            secondary_alignments = {}

        primary_positions = [np.array(p, dtype=np.float32) for p in positions]
        secondary_positions: list[npt.NDArray[np.float32]] = []
        edges = []
        primary_neighbors = defaultdict(list)
        secondary_neighbors = defaultdict(list)
        for (i1, p1), (i2, p2) in itertools.combinations(
            enumerate(primary_positions),
            2,
        ):
            primary_neighbors[i1].append(len(secondary_positions))
            primary_neighbors[i2].append(len(secondary_positions))
            secondary_neighbors[len(secondary_positions)].append(i1)
            secondary_neighbors[len(secondary_positions)].append(i2)
            edges.append((i1, len(secondary_positions)))
            edges.append((i2, len(secondary_positions)))
            secondary_positions.append((p1 + p2) / 2)

        primary_vertices = []
        for i, position in enumerate(primary_positions):
            anchor_neighbor = primary_neighbors[i][
                primary_alignments.get(i, 0)
            ]
            neighbor_position = secondary_positions[anchor_neighbor]
            primary_vertices.append(
                RotationVertex(
                    position=position,
                    rotation_axis=normalize_vector(position),
                    rotation_target=neighbor_position - position,
                )
            )

        secondary_vertices = []
        for i, position in enumerate(secondary_positions):
            anchor_neighbor = secondary_neighbors[i][
                secondary_alignments.get(i, 0)
            ]
            neighbor_position = primary_positions[anchor_neighbor]
            secondary_vertices.append(
                RotationVertex(
                    position=position,
                    rotation_axis=normalize_vector(position),
                    rotation_target=neighbor_position - position,
                )
            )

        return TwoGraph(
            primary_vertices=primary_vertices,
            secondary_vertices=secondary_vertices,
        )

    def build(
        self,
        primary: RotationBuildingBlock,
        secondary: RotationBuildingBlock,
    ) -> Cage:
        return build(
            assignments=itertools.chain(
                ((vertex, primary) for vertex in self.primary_vertices),
                ((vertex, secondary) for vertex in self.secondary_vertices),
            ),
        )


def four_plus_six(
    primary_alignments: dict[int, int] | None = None,
    secondary_alignments: dict[int, int] | None = None,
) -> TwoGraph:
    return TwoGraph.from_primary_vertex_positions(
        positions=[
            [0, 0, np.sqrt(6) / 2],
            [-1, -np.sqrt(3) / 3, -np.sqrt(6) / 6],
            [1, -np.sqrt(3) / 3, -np.sqrt(6) / 6],
            [0, 2 * np.sqrt(3) / 3, -np.sqrt(6) / 6],
        ],
        primary_alignments=primary_alignments,
        secondary_alignments=secondary_alignments,
    )


def build(
    assignments: Iterable[Assignment],
) -> Cage:
    atoms: list[Atom] = []
    integer_bonds: list[IntegerBond] = []
    dative_bonds: list[DativeBond] = []
    position_matrix = []

    for vertex, building_block in assignments:
        new_atoms = [
            replace(atom, id=atom.id + len(atoms))
            for atom in building_block.atoms
        ]
        atoms.extend(new_atoms)
        integer_bonds.extend(
            replace(
                bond,
                atom1=new_atoms[bond.atom1.id],
                atom2=new_atoms[bond.atom2.id],
            )
            for bond in building_block.integer_bonds
        )
        dative_bonds.extend(
            replace(
                bond,
                atom1=new_atoms[bond.atom1.id],
                atom2=new_atoms[bond.atom2.id],
            )
            for bond in building_block.dative_bonds
        )
        position_matrix.append(
            replace(vertex, position=vertex.position * 10).place(
                matrix=building_block.position_matrix,
                position_anchor=building_block.position_anchor,
                rotation_anchor_axis=building_block.rotation_anchor_axis,
                rotation_anchor_target=building_block.rotation_anchor_target,
            )
        )

    return Cage(
        atoms=atoms,
        integer_bonds=integer_bonds,
        dative_bonds=dative_bonds,
        position_matrix=np.vstack(position_matrix, dtype=np.float32),
    )
