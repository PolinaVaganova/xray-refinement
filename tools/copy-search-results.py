import shutil
from pathlib import Path
from typing import List


def get_pdb_codes(prefix: Path) -> List[str]:
    return sorted(
        filename.stem[:4].lower()
        for filename in (prefix / "structure_factors").glob("*-sf_P1.mtz")
    )


def main(search_results_dir: Path, output_dir: Path):
    pdb_codes = get_pdb_codes(search_results_dir)

    for pdb_code in pdb_codes:
        orig_pdb = search_results_dir / "db" / f"pdb{pdb_code}.ent"
        orig_mtz = search_results_dir / "structure_factors" / f"{pdb_code}-sf_P1.mtz"

        pdb_dir = output_dir / pdb_code

        if orig_mtz.exists() and orig_pdb.exists():
            pdb_dir.mkdir(exist_ok=True)
            shutil.copy(orig_pdb, pdb_dir / f"{pdb_code}.pdb")
            shutil.copy(orig_mtz, pdb_dir / f"{pdb_code}.mtz")


if __name__ == "__main__":
    main(
        search_results_dir=Path("~/bionmr/oleg/amber_xray_test_DB/").expanduser(),
        output_dir=Path.cwd() / "data" / "input",
    )
