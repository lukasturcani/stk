import argparse
import logging
import pathlib
import shutil

import numpy as np
import pytest

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
logger = logging.getLogger(__name__)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "test_path",
        type=pathlib.Path,
        help="where are the tests?",
    )
    parser.add_argument(
        "output_path",
        type=pathlib.Path,
        help="where to save structures to",
    )

    return parser.parse_args()


def main() -> None:
    """Run script."""
    args = _parse_args()
    tests = args.test_path

    pos_mats = (
        tests / "molecular/molecules/molecule/fixtures/position_matrices/"
    )
    old_pos_mats = (
        tests / "molecular/molecules/molecule/fixtures/old_position_matrices/"
    )
    args.output_path.mkdir(parents=True, exist_ok=True)

    # Copy from current to old and remove all in current while maintaining
    # directory for tests.
    pos_mats.rename(old_pos_mats)
    pos_mats.mkdir(parents=True, exist_ok=False)

    # Run pytest, but only those that will make the position matrices.
    pytest.main(["-k test_get_plane_normal"])

    # Compare.
    new_numpy_files = list(sorted(pos_mats.glob("*.npy")))
    old_numpy_files = list(sorted(old_pos_mats.glob("*.npy")))
    assert len(new_numpy_files) == len(old_numpy_files)

    rmsds = []
    for new_file, old_file in zip(
        new_numpy_files,
        old_numpy_files,
        strict=True,
    ):
        assert new_file.stem == old_file.stem
        new_array = np.load(new_file)
        old_array = np.load(old_file)
        assert len(new_array) == len(old_array)

        rmsd = np.sqrt(np.mean((new_array - old_array) ** 2))

        rmsds.append(rmsd)
        if rmsd >= 0.1:
            logger.info(
                "very large rmsd for %s is %s", old_file.name, round(rmsd, 3)
            )
            output_file = args.output_path / f"{old_file.stem}_0.xyz"
            top_line = f"{len(old_array)}\n"
            comment_line = "\n"
            pos_lines = []
            for pos in old_array:
                pos_lines.append(f"C {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
            with output_file.open("w") as f:
                f.write(top_line)
                f.write(comment_line)
                f.writelines(pos_lines)

            new_output_file = args.output_path / f"{new_file.stem}_1.xyz"
            top_line = f"{len(new_array)}\n"
            comment_line = "\n"
            pos_lines = []
            for pos in new_array:
                pos_lines.append(f"C {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
            with new_output_file.open("w") as f:
                f.write(top_line)
                f.write(comment_line)
                f.writelines(pos_lines)

    logger.info("average rmsd is %s", round(np.mean(rmsds), 3))

    # Delete test and return back to normal.
    shutil.rmtree(pos_mats)
    old_pos_mats.rename(pos_mats)
    logger.info("do not forget to delete %s", args.output_path)


if __name__ == "__main__":
    main()
