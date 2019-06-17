import stk
from os.path import join
import numpy as np


def test_xtb_extractor():
    known_output_file = join('../data', 'xtb_energy.output')

    # Get properties from output_file.
    xtbext = stk.XTBExtractor(output_file=known_output_file)
    total_energy = xtbext.total_energy()
    homo_lumo_gap = xtbext.homo_lumo_gap()
    fermi_level = xtbext.fermi_level()
    homo_lumo_occ = xtbext.homo_lumo_occ()
    qonly_dipole_moment = xtbext.qonly_dipole_moment()
    full_dipole_moment = xtbext.full_dipole_moment()
    qonly_quadrupole_moment = xtbext.qonly_quadrupole_moment()
    qdip_quadrupole_moment = xtbext.qdip_quadrupole_moment()
    full_quadrupole_moment = xtbext.full_quadrupole_moment()
    total_free_energy = xtbext.total_free_energy()
    frequencies = xtbext.frequencies()

    assert np.isclose(
        total_energy, -76.323188311664, rtol=0, atol=1.e-8
    )
    assert np.isclose(
        homo_lumo_gap, 2.169914122867, rtol=0, atol=1.e-8
    )
    assert np.isclose(fermi_level, -8.4811, rtol=0, atol=1.e-8)
    assert homo_lumo_occ['HOMO'][0] == 70
    assert homo_lumo_occ['LUMO'][0] == 71
    assert np.isclose(
        homo_lumo_occ['HOMO'][2], -9.5661, rtol=0, atol=1.e-3
    )
    assert np.isclose(
        homo_lumo_occ['LUMO'][2], -7.3961, rtol=0, atol=1.e-3
    )
    assert np.isclose(
        total_free_energy, -75.905031656056, rtol=0, atol=1.e-8
    )
    assert np.allclose(
        qonly_dipole_moment, [0.041, -0.064, 1.209], rtol=0, atol=1.e-2
    )
    known_ = [-0.174, 0.224, 1.354, 3.517]
    assert np.allclose(
        full_dipole_moment, known_, rtol=0, atol=1.e-2
    )
    known_ = [3.746, 8.952, 4.757, -23.734, 0.600, -8.503]
    assert np.allclose(
        qonly_quadrupole_moment, known_, rtol=0, atol=1.e-2
    )
    known_ = [4.992, 11.616, 6.936, -22.340, 4.353, -11.927]
    assert np.allclose(
        qdip_quadrupole_moment, known_, rtol=0, atol=1.e-2
    )
    known_ = [5.134, 12.073, 7.577, -21.996, 4.439, -12.711]
    assert np.allclose(
        full_quadrupole_moment, known_, rtol=0, atol=1.e-2
    )
    known_ = [
        -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 2.93, 7.16, 9.99, 13.28, 17.23,
        21.94, 25.18, 32.49, 43.17, 53.79, 56.81, 68.99, 73.16, 82.13,
        84.81, 99.5, 108.7, 111.38, 142.18, 146.85, 175.72, 190.38,
        205.29, 229.83, 242.85, 251.96, 281.39, 283.88, 295.98, 310.92,
        322.62, 344.52, 361.61, 367.22, 384.87, 422.02, 440.56, 499.64,
        508.19, 522.54, 551.47, 578.47, 603.28, 606.42, 635.94, 673.35,
        713.33, 770.62, 785.53, 818.46, 834.27, 840.81, 853.73, 865.85,
        876.21, 888.31, 894.59, 901.26, 907.05, 941.82, 947.98, 951.35,
        952.66, 959.57, 974.27, 986.86, 1013.48, 1021.7, 1038.77,
        1041.33, 1056.1, 1065.43, 1067.97, 1074.06, 1079.15, 1082.84,
        1087.1, 1106.08, 1109.94, 1111.64, 1119.43, 1123.42, 1134.23,
        1145.88, 1153.22, 1155.18, 1158.01, 1168.0, 1178.86, 1192.58,
        1213.2, 1233.85, 1235.34, 1242.96, 1247.8, 1248.26, 1252.7,
        1253.31, 1260.26, 1267.08, 1280.28, 1314.27, 1316.0, 1317.92,
        1330.69, 1331.0, 1333.01, 1340.74, 1341.97, 1346.33, 1350.62,
        1356.74, 1373.62, 1378.21, 1402.28, 1446.54, 1449.61, 1454.18,
        1460.77, 1465.75, 1466.32, 1469.06, 1480.2, 1489.47, 1489.49,
        1495.76, 1552.2, 1729.76, 1737.41, 1761.18, 1764.79, 1768.34,
        1796.44, 2748.06, 2810.05, 2833.24, 2850.66, 2852.9, 2860.05,
        2868.69, 2905.77, 2909.07, 2932.37, 2932.85, 2933.49, 2934.21,
        2936.47, 2938.02, 2940.2, 2971.16, 2982.15, 2983.83, 2987.45,
        2989.03, 2991.81, 2992.93, 2994.95, 2995.59, 3002.24, 3003.05,
        3003.7, 3006.08, 3009.6, 3352.09, 3405.37
    ]
    assert np.allclose(
        frequencies, known_, rtol=0, atol=1.e-8
    )
