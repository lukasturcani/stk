#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.
import stk
import numpy as np


def main():
    print('testing script')
    metal = stk.BuildingBlock(
        '[Pd+2]',
        functional_groups=None,
        embed=False
    )
    metal.set_position_matrix(np.array([0, 0, 0]))
    metal_coord_info = {
        0: {
            'atom_ids': 0,
            'bonder_ids': 0,
            'deleter_ids': None
        },
        1: {
            'atom_ids': 0,
            'bonder_ids': 0,
            'deleter_ids': None
        },
        2: {
            'atom_ids': 0,
            'bonder_ids': 0,
            'deleter_ids': None
        },
        3: {
            'atom_ids': 0,
            'bonder_ids': 0,
            'deleter_ids': None
        },
    }
    metal = stk.assign_metal_fgs(
        building_block=metal,
        coordination_info=metal_coord_info
    )
    ligand = stk.BuildingBlock(
        'c1cc(-c2ccc(-c3ccncc3)cc2)ccn1',
        functional_groups=['pyridine_N_metal']
    )
    print(metal)
    print(ligand)

    sqpl = stk.metal_complex.SquarePlanar()
    pdl2_sqpl_complex = stk.ConstructedMolecule(
        building_blocks=[metal, ligand],
        topology_graph=sqpl,
        building_block_vertices={
            metal: sqpl.vertices[0],
            ligand: sqpl.vertices[1:2]
        }
    )
    print(pdl2_sqpl_complex)
    pdl2_sqpl_complex.write('metal_complex.mol')
    pdl2_sqpl_complex.write('metal_complex.pdb')

    print('testing script done')


if __name__ == "__main__":
    main()
