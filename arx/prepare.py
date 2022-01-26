import os
import subprocess
import tempfile
import typing

import gemmi

from .utils import chdir, check_call

AnyPath = typing.Union[str, bytes, os.PathLike]


def add_ions(
    st: gemmi.Structure, ion: gemmi.Structure, count: int, *, n_protein_atoms: int
) -> gemmi.Structure:
    with tempfile.TemporaryDirectory() as tmpdir:
        with chdir(tmpdir):
            write_pdb(st, "input.pdb")
            write_pdb(ion, "ion.pdb")
            check_call(
                [
                    "AddToBox",
                    "-c",
                    "input.pdb",
                    "-a",
                    "ion.pdb",
                    "-na",
                    str(count),
                    "-o",
                    "result.pdb",
                    "-P",
                    str(n_protein_atoms),
                    "-RP",
                    str(3.0),
                    "-RW",
                    str(6.0),
                    "-G",
                    str(0.2),
                    "-V",
                    str(1),
                ]
            )
            return read_pdb("result.pdb")


def add_water(st: gemmi.Structure, water: gemmi.Structure) -> gemmi.Structure:
    n_non_water_atoms = st[0].count_atom_sites()
    count = estimate_water_molecules(st)
    with tempfile.TemporaryDirectory() as tmpdir:
        with chdir(tmpdir):
            write_pdb(st, "input.pdb")
            write_pdb(water, "water.pdb")
            check_call(
                [
                    "AddToBox",
                    "-c",
                    "input.pdb",
                    "-a",
                    "water.pdb",
                    "-na",
                    str(count),
                    "-o",
                    "result.pdb",
                    "-P",
                    str(n_non_water_atoms),
                    "-RP",
                    str(3.0),
                    "-RW",
                    str(3.0),
                    "-G",
                    str(0.2),
                    "-V",
                    str(1),
                ]
            )
            return read_pdb("result.pdb")


def calc_total_charge(st: gemmi.Structure) -> int:
    tleap_in = """
source oldff/leaprc.ff14SB

HOH = SPC
WAT = SPC

loadAmberParams frcmod.ionsjc_spce
loadAmberParams frcmod.ions1lm_1264_spce
loadAmberParams frcmod.spce

wbox = loadpdb input.pdb
setBox wbox vdw 1.0
saveamberparm wbox wbox.prmtop wbox.inpcrd
savepdb wbox wbox.dry.pdb
quit
"""

    with tempfile.TemporaryDirectory() as tmpdir:
        with chdir(tmpdir):
            write_pdb(st, "input.pdb")
            with open("tleap.in", "w") as f:
                f.write(tleap_in)

            with open("parmed.in", "w") as f:
                f.write("summary\n")
                f.write("quit\n")

            check_call(["tleap", "-s", "-f", "tleap.in"])

            parmed_output = subprocess.check_output(
                [
                    "parmed",
                    "wbox.prmtop",
                    "parmed.in",
                ]
            ).decode("utf-8")

    for line in parmed_output.split("\n"):
        if "Total charge (e-)" in line:
            return int(float(line.split(":")[-1].strip().split()[0]))

    raise RuntimeError("Can't find Total charge in parmed output")


def estimate_water_molecules(st: gemmi.Structure) -> int:
    psv = 0.74
    MW = estimate_weight(st)
    nau = 1
    # 100K->RT volume expansion coefficient
    expand = 1.0  # no expansion
    # ==============================================================================
    # Volume of protein
    VP = MW * nau * psv * 1e24 / 6.02e23
    V = st.cell.volume
    # print("Water content: %.3f" % (1 - VP / V))
    V *= expand
    # estimated water molecule numbers
    nwat = (V - VP) * 1e-24 * 1.00 / 18.0 * 6.02e23
    assert nwat > 10, f"Too low number of water molecules: {nwat}"
    return int(nwat)


def estimate_weight(st: gemmi.Structure) -> float:
    mass = 0.0
    aname_map = {"Na+": 22.989769, "Cl-": 35.453}
    element_map = {
        "P": 30.973762,
        "H": 1.00794,
        "C": 12.0107,
        "S": 32.065,
        "O": 15.999,
        "N": 14.0067,
    }
    for m in st:
        for c in m:
            for r in c:
                for a in r:
                    if a.name in aname_map:
                        mass += aname_map[a.name]
                    else:
                        mass += element_map[a.name[0]]

    return mass


def expand_non_crystallographic_symmetries(st: gemmi.Structure) -> gemmi.Structure:
    result = st.clone()
    result.expand_ncs(how=gemmi.HowToNameCopiedChain.AddNumber)
    return result


def expand_crystallographic_symmetries(st: gemmi.Structure) -> gemmi.Structure:
    g: gemmi.SpaceGroup = gemmi.find_spacegroup_by_name(st.spacegroup_hm)

    def apply_to_chain(chain, op, cell):
        for residue in chain:
            for atom in residue:
                fractional = cell.fractionalize(atom.pos)
                transformed = gemmi.Fractional(*op.apply_to_xyz(fractional.tolist()))
                cartesian = cell.orthogonalize(transformed)
                atom.pos = cartesian

    result = gemmi.Structure()

    for model in st:
        new_model = gemmi.Model(model.name)
        for op in g.operations():
            for chain in model:
                new_model.add_chain(chain, pos=-1)
                new_chain = new_model[-1]
                apply_to_chain(new_chain, op, cell=st.cell)
        result.add_model(new_model, pos=-1)
    result.cell = st.cell
    result.spacegroup_hm = st.spacegroup_hm
    return result


def neutralize_with_ions(
    st: gemmi.Structure, negative_ion: gemmi.Structure, positive_ion: gemmi.Structure
) -> gemmi.Structure:
    total_charge = calc_total_charge(st)
    if total_charge > 0:
        ion = negative_ion
    else:
        ion = positive_ion
    count = abs(total_charge)
    n_protein_atoms = st[0].count_atom_sites()
    return add_ions(st, ion, count=count, n_protein_atoms=n_protein_atoms)


def add_missing_atoms(st: gemmi.Structure) -> gemmi.Structure:
    with tempfile.TemporaryDirectory() as tmpdir:
        with chdir(tmpdir):
            input_pdb = "input.pdb"
            result_pdb = "result.pdb"
            write_pdb(st, input_pdb)
            check_call(
                [
                    "pdb4amber",
                    "-d",
                    "--reduce",
                    "--add-missing-atoms",
                    "-i",
                    input_pdb,
                    "-o",
                    result_pdb,
                ]
            )
            result = read_pdb(result_pdb)
    result.cell = st.cell
    result.spacegroup_hm = st.spacegroup_hm
    return result


def apply_to_atoms(
    st: gemmi.Structure, atom_mutator: typing.Callable[[gemmi.Atom], None]
) -> gemmi.Structure:
    result = st.clone()
    for model in result:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom_mutator(atom)
    return result


def set_b_factors_to(st: gemmi.Structure, value: float = 0) -> gemmi.Structure:
    def setter(atom: gemmi.Atom):
        atom.b_iso = value

    return apply_to_atoms(st, atom_mutator=setter)


def add_missing_b_factors(
    st: gemmi.Structure, reference: gemmi.Structure
) -> gemmi.Structure:
    index: gemmi.NeighborSearch = gemmi.NeighborSearch(reference, 3.0)
    index.populate()

    result = set_b_factors_to(st, 0)
    total_assigned = 0
    min_distance = 0
    max_distance = 1.0
    while True:
        n_assigned = 0
        for model in result:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        if atom.b_iso > 0:
                            continue
                        closest = index.find_neighbors(atom, min_distance, max_distance)
                        best_dist = 99
                        best_ref_at = None
                        for mark in closest:
                            ref_at = mark.to_cra(reference[0]).atom
                            dist = ref_at.pos.dist(atom.pos)
                            if dist < best_dist:
                                best_ref_at = ref_at
                                best_dist = dist
                        if best_ref_at is not None:
                            atom.b_iso = best_ref_at.b_iso
                            n_assigned += 1
        total_assigned += n_assigned
        min_distance = max_distance
        max_distance += 1

        if n_assigned == 0:
            break

    return result


def add_missing_occupancies(
    st: gemmi.Structure, reference: gemmi.Structure
) -> gemmi.Structure:
    def setter(atom: gemmi.Atom):
        atom.occ = 1.0

    return apply_to_atoms(st, setter)


def remove_alternative_conformations(st: gemmi.Structure) -> gemmi.Structure:
    result = st.clone()
    result.remove_alternative_conformations()
    return result


def remove_ligands_and_water(st: gemmi.Structure) -> gemmi.Structure:
    result = st.clone()
    result.remove_ligands_and_waters()
    return result


def remove_water(st: gemmi.Structure) -> gemmi.Structure:
    result = st.clone()
    result.remove_waters()
    return result


def read_pdb(path: AnyPath) -> gemmi.Structure:
    return gemmi.read_pdb(str(path), split_chain_on_ter=True)


def write_pdb(st: gemmi.Structure, path: AnyPath):
    return st.write_pdb(str(path), numbered_ter=False, ter_ignores_type=True)


def copy_coordinates(
    st: gemmi.Structure, reference: gemmi.Structure
) -> gemmi.Structure:
    result = st.clone()
    for (_, _, st_r, st_a), (_, _, ref_r, ref_a) in zip(
        iterate_over_atoms(result), iterate_over_atoms(reference)
    ):
        assert st_r.name == ref_r.name
        assert st_a.name == ref_a.name
        st_a.pos = ref_a.pos
    return result


def iterate_over_atoms(
    st: gemmi.Structure,
) -> typing.Iterator[typing.Tuple[gemmi.Model, gemmi.Chain, gemmi.Residue, gemmi.Atom]]:
    for model in st:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    yield model, chain, residue, atom
