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
)


class PipeInput:
    def __init__(self, arg):
        self.arg = arg

    def do(self, func, *args, **kwargs):
        print(func.__name__, "...")
        return PipeInput(func(self.arg, *args, **kwargs))


def main():
    input_structure = read_pdb("2msi.pdb")
    assert len(input_structure) == 1

    spce = read_pdb("spce.pdb")
    cl = read_pdb("cl.pdb")
    na = read_pdb("na.pdb")

    (
        PipeInput(input_structure)
        .do(remove_ligands_and_water)
        .do(remove_alternative_conformations)
        .do(add_missing_atoms)
        .do(add_missing_b_factors, reference=input_structure)
        .do(add_missing_occupancies, reference=input_structure)
        .do(expand_non_crystallographic_symmetries)
        .do(expand_crystallographic_symmetries)
        .do(add_water, water=spce)
        .do(neutralize_with_ions, negative_ion=cl, positive_ion=na)
    )


if __name__ == "__main__":
    main()
