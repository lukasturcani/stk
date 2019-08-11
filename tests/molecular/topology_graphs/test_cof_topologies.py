# import stk
# import os
# from os.path import join
#
# from ...._test_utilities import _test_dump_and_load
#
#
# test_dir = 'cof_topology_tests_output'
# if not os.path.exists(test_dir):
#     os.mkdir(test_dir)
#
#
# def test_alignments(amine2_alt3, aldehyde3_alt3, aldehyde4_alt1):
#     for i in range(4):
#         v0 = stk.cof.Kagome.vertices[0]
#         v = stk.cof.Kagome.vertices[-1]
#         cof = stk.ConstructedMolecule(
#             building_blocks=[amine2_alt3, aldehyde4_alt1],
#             topology_graph=stk.cof.Kagome(
#                 lattice_size=(3, 3, 1),
#                 vertex_alignments={
#                     v0: v0.edges[i],
#                     v: v.edges[i % 2]
#                 }
#             )
#         )
#         cof.write(join(test_dir, f'aligning_{i}_{i%2}.sdf'))
#
#
# def test_multi_bb():
#     assert False
#
#
# def _test_construction(cof, num_expected_bbs):
#     path = join(
#         test_dir, f'{cof.topology_graph.__class__.__name__}.mol'
#     )
#     cof.write(path)
#
#     is_periodic = any(
#         any(d != 0 for d in bond.periodicity) for bond in cof.bonds
#     )
#     for bb in cof.get_building_blocks():
#         assert cof.building_block_counter[bb] == num_expected_bbs[bb]
#         # This test only holds true when each building block is
#         # involved in every construction bond and the cof is
#         # periodic.
#         if len(num_expected_bbs) < 3 and is_periodic:
#             assert (
#                 len(cof.construction_bonds) ==
#                 cof.building_block_counter[bb] * len(bb.func_groups)
#             )
#
#     num_deleters = sum(
#         len(fg.deleters)*cof.building_block_counter[bb]
#         for bb in cof.get_building_blocks() for fg in bb.func_groups
#     )
#     num_bb_atoms = sum(
#         len(bb.atoms)*cof.building_block_counter[bb]
#         for bb in cof.get_building_blocks()
#     )
#     num_bb_bonds = sum(
#         len(bb.bonds)*cof.building_block_counter[bb]
#         for bb in cof.get_building_blocks()
#     )
#     # Check that the correct number of bonds got made.
#     assert (
#         len(cof.construction_bonds) == len(cof.topology_graph.edges)
#     )
#     # Check correct total number of atoms.
#     if is_periodic:
#         assert len(cof.atoms) == num_bb_atoms - num_deleters
#     # Check correct total number of bonds.
#     if is_periodic:
#         assert (
#             len(cof.bonds) ==
#             num_bb_bonds + len(cof.construction_bonds) - num_deleters
#         )
#
#
# def test_topologies(
#     tmp_honeycomb,
#     tmp_periodic_honeycomb,
#     tmp_kagome,
#     tmp_periodic_kagome,
#     tmp_hexagonal,
#     tmp_periodic_hexagonal,
#     tmp_square,
#     tmp_periodic_square,
#     tmp_linkerless_honeycomb,
#     tmp_periodic_linkerless_honeycomb
# ):
#     cofs = (
#         (tmp_honeycomb, 3*9, 2*9),
#         (tmp_periodic_honeycomb, 3*9, 2*9)
#         (tmp_kagome, 6*9, 3*9),
#         (tmp_periodic_kagome, 6*9, 3*9),
#         (tmp_hexagonal, 12*9, 4*9),
#         (tmp_periodic_hexagonal, 12*9, 4*9)
#         (tmp_square, 2*9, 1*9),
#         (tmp_periodic_square, 2*9, 1*9)
#         (tmp_linkerless_honeycomb, 1*9, 1*9)
#         (tmp_periodic_linkerless_honeycomb, 1*9, 1*9)
#     )
#
#     for cof, num_linkers, num_building_blocks in cofs:
#         linker, building_block = sorted(
#             cof.get_building_blocks(),
#             key=lambda bb: len(bb.func_groups)
#         )
#         num_expected_bbs = {
#             linker: num_linkers,
#             building_block: num_building_blocks
#         }
#         _test_construction(cof, num_expected_bbs)
#         _test_dump_and_load(cof)
