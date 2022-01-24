from pathlib import Path

import gemmi

from arx.amber import create_topology_and_input
from arx.prepare import (
    add_missing_atoms,
    add_missing_b_factors,
    add_missing_occupancies,
    add_water,
    expand_crystallographic_symmetries,
    expand_non_crystallographic_symmetries,
    neutralize_with_ions,
    read_pdb,
    remove_alternative_conformations,
    remove_ligands_and_water,
    write_pdb,
)


class StructurePipeInput:
    DEBUG = True

    def __init__(self, st: gemmi.Structure, *, debug=DEBUG, seq=0):
        self.st: gemmi.Structure = st
        self.debug = debug
        self.seq = seq

    def do(self, func, *args, **kwargs):
        if self.DEBUG:
            write_pdb(self.st, f"{self.seq:02d}_{func.__name__}.before.pdb")
        result: gemmi.Structure = func(self.st, *args, **kwargs)
        if self.DEBUG:
            write_pdb(result, f"{self.seq:02d}_{func.__name__}.after.pdb")
        return StructurePipeInput(result, debug=self.debug, seq=self.seq + 1)


def main():
    input_structure = read_pdb("2msi.pdb")
    assert len(input_structure) == 1

    spce = read_pdb("spce.pdb")
    cl = read_pdb("cl.pdb")
    na = read_pdb("na.pdb")

    parm7_path = (Path().cwd() / "wbox.parm7").absolute()
    rst7_path = (Path().cwd() / "wbox.rst7").absolute()

    (
        StructurePipeInput(input_structure)
        .do(remove_ligands_and_water)
        .do(remove_alternative_conformations)
        .do(add_missing_atoms)
        .do(add_missing_b_factors, reference=input_structure)
        .do(add_missing_occupancies, reference=input_structure)
        .do(expand_non_crystallographic_symmetries)
        .do(expand_crystallographic_symmetries)
        .do(add_water, water=spce)
        .do(neutralize_with_ions, negative_ion=cl, positive_ion=na)
        .do(create_topology_and_input, parm7_path=parm7_path, rst7_path=rst7_path)
    )


if __name__ == "__main__":
    main()
