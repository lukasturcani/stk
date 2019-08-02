def test_cage_complex(amine2, amine2_alt1, aldehyde3, chained_c60):
    c = stk.Cage([amine2, amine2_alt1, aldehyde3],
                 stk.FourPlusSix(bb_positions={
                     amine2: [5],
                     amine2_alt1: [0, 1, 2, 3, 4],
                     aldehyde3: [0, 1, 2, 3]
                 }))

    n = 3
    for i in range(n):
        cage_complex = stk.CageComplex(
            [c, chained_c60],
            stk.CageWithGuest(axis=[1, 0, 0],
                              angle=2*np.pi*i/n,
                              displacement=[2*i, 0, 0])
        )
        cage_complex.write(join(test_dir, f'cage_with_guest_{i}.mol'))
