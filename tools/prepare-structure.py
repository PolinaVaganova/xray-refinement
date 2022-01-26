import logging
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

logger = logging.getLogger("prepare-structure")


class StructurePipeInput:
    DEBUG = True

    def __init__(self, st: gemmi.Structure, *, debug=DEBUG, seq=0):
        self.st: gemmi.Structure = st
        self.debug = debug
        self.seq = seq

    def do(self, func, *args, **kwargs):
        logger.info(f"{func.__name__}...")
        if self.DEBUG:
            write_pdb(self.st, self.fmt_debug_path(f"{func.__name__}.before"))
        result: gemmi.Structure = func(self.st, *args, **kwargs)
        if self.DEBUG:
            write_pdb(result, self.fmt_debug_path(f"{func.__name__}.after"))
        return StructurePipeInput(result, debug=self.debug, seq=self.seq + 1)

    def fmt_debug_path(self, name: str, ext="pdb") -> Path:
        debug_dir = Path.cwd() / "debug"
        debug_dir.mkdir(exist_ok=True)
        return debug_dir / f"{self.seq:02d}_{name}.{ext}"


def create_parm7_rst7_from(
    structure: gemmi.Structure,
    solvent: gemmi.Structure,
    positive_ion: gemmi.Structure,
    negative_ion: gemmi.Structure,
    parm7_path: Path,
    rst7_path: Path,
) -> gemmi.Structure:
    return (
        StructurePipeInput(structure)
        .do(remove_ligands_and_water)
        .do(remove_alternative_conformations)
        .do(add_missing_atoms)
        .do(add_missing_b_factors, reference=structure)
        .do(add_missing_occupancies, reference=structure)
        .do(expand_non_crystallographic_symmetries)
        .do(expand_crystallographic_symmetries)
        .do(add_water, water=solvent)
        .do(neutralize_with_ions, negative_ion=negative_ion, positive_ion=positive_ion)
        .do(create_topology_and_input, parm7_path=parm7_path, rst7_path=rst7_path)
    ).st


def main():
    logging.basicConfig()
    logger.setLevel(logging.INFO)

    data = Path.cwd() / "data"
    sm = data / "small-molecules"
    wbox = create_parm7_rst7_from(
        structure=read_pdb(data / "input" / "2msi" / "2msi.pdb"),
        positive_ion=read_pdb(sm / "na.pdb"),
        negative_ion=read_pdb(sm / "cl.pdb"),
        solvent=read_pdb(sm / "spce.pdb"),
        parm7_path=(Path().cwd() / "wbox.parm7").absolute(),
        rst7_path=(Path().cwd() / "wbox.rst7").absolute(),
    )
    write_pdb(wbox, (Path().cwd() / "wbox.pdb").absolute())


if __name__ == "__main__":
    main()
