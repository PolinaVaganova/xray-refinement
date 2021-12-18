from pathlib import Path
from typing import List
import shutil


def get_pdb_codes() -> List[str]:
    return [
        "149l",
        "171l",
        "1ae2",
        "1ae3",
        "1ail",
        "1ame",
        "1anu",
        "1bkl",
        "1bmg",
        "1dt4",
        "1du5",
        "1k40",
        "1kem",
        "1loz",
        "1mjc",
        "1o9h",
        "1q2y",
        "1rn7",
        "1smt",
        "1uue",
        "1w45",
        "1wdx",
        "1x6j",
        "1yib",
        "1z27",
        "2aak",
        "2eql",
        "2fht",
        "2fkl",
        "2j7i",
        "2jee",
        "2msi",
        "2nn4",
        "2nsb",
        "2o85",
        "2o87",
        "2o89",
        "2ont",
        "2qdo",
        "2snw",
        "3dvp",
        "3fis",
        "3hpm",
        "3k9p",
        "3le4",
        "3m3t",
        "3ndf",
        "3q2c",
        "3qd7",
        "3rd3",
        "3tsv",
        "3u4z",
        "3zq7",
        "3zy1",
        "3zye",
        "4bhc",
        "4c0m",
        "4c86",
        "4f17",
        "4f26",
        "4fis",
        "4hll",
        "4niq",
        "4oyc",
        "4qfq",
        "4ug3",
        "4wfw",
        "4x37",
        "5a7l",
        "5arj",
        "5ewr",
        "5f6a",
        "5fd7",
        "5h79",
        "5jqz",
        "5lhx",
        "5t8l",
        "5t8n",
        "5teo",
        "5tut",
        "5xbh",
        "6cyr",
        "6dz6",
        "6msi",
    ]


def main(prefix: Path):
    output_dir = Path(__file__).parent / "data" / "input"
    pdb_codes = get_pdb_codes()

    for pdb_code in pdb_codes:
        orig_d_path = (
            prefix / f"{pdb_code}_asu_5_rounds_d/4phenix_final_structure_bf.pdb"
        )
        orig_s_path = (
            prefix / f"{pdb_code}_asu_5_rounds_md1/4phenix_final_structure_bf_asu.pdb"
        )
        orig_mtz = prefix / f"{pdb_code}_uc_5_rounds_d/real.mtz"

        pdb_dir = output_dir / pdb_code
        pdb_dir.mkdir(exist_ok=True)
        shutil.copy(orig_d_path, pdb_dir / f"{pdb_code}.D-set.pdb")
        shutil.copy(orig_s_path, pdb_dir / f"{pdb_code}.S-set.pdb")
        shutil.copy(orig_mtz, pdb_dir / f"{pdb_code}.mtz")


if __name__ == "__main__":
    main(prefix=Path("/home/sergei/bionmr/oleg/structures2refine/"))
