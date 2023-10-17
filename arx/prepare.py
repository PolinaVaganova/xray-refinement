import os
import subprocess
import tempfile
import typing

import gemmi

from .utils import chdir, check_call

AnyPath = typing.Union[str, bytes, os.PathLike]

tleap_in = """
source leaprc.protein.ff14SB
source leaprc.DNA.OL15
source leaprc.RNA.OL3

HOH = SPC
WAT = SPC

loadAmberParams frcmod.ionsjc_spce
loadAmberParams frcmod.ions1lm_1264_spce
loadAmberParams frcmod.spce

wbox = loadpdb input.pdb
set default nocenter on
setBox wbox vdw 1.0
saveamberparm wbox wbox.prmtop wbox.inpcrd
savepdb wbox wbox.dry.pdb
quit
"""


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
    aname_map = dict([(n, gemmi.Element(n).weight) for n in ["Na+", "Cl-"]])
    element_map = dict([(n, gemmi.Element(n).weight) for n in ["P", "H", "C", "S", "O", "N"]])
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
    result.spacegroup_hm = 'P1'  # st.spacegroup_hm
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

    result = set_b_factors_to(st, -1)
    total_assigned = 0
    shell_step = 0.5
    min_distance = 0
    max_distance = shell_step
    is_first_pass = True
    while True:
        n_assigned = 0
        unassigned_exist = False
        for model in result:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        if atom.b_iso > 0:
                            continue
                        unassigned_exist = True
                        closest = index.find_neighbors(atom, min_distance, max_distance)
                        max_b_factor = -1
                        for mark in closest:
                            if is_first_pass:
                                ref_at = mark.to_cra(reference[0]).atom
                            else:
                                ref_at = mark.to_cra(model).atom
                            dist = ref_at.pos.dist(atom.pos)
                            if dist < max_distance:
                                max_b_factor = max(max_b_factor, ref_at.b_iso)
                        if max_b_factor > 0:
                            atom.b_iso = max_b_factor
                            n_assigned += 1
        if is_first_pass:
            index: gemmi.NeighborSearch = gemmi.NeighborSearch(result, 3.0)
            index.populate()
            is_first_pass = False
        total_assigned += n_assigned
        min_distance = max_distance
        max_distance += shell_step
        shell_step = 2.0

        if n_assigned == 0 and not unassigned_exist:
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


def remove_empty_chains(st: gemmi.Structure) -> gemmi.Structure:
    result = st.clone()
    result.remove_empty_chains()
    return result


def remove_hydrogens(
    st: gemmi.Structure, reference: gemmi.Structure
) -> gemmi.Structure:
    result = st.clone()
    for (_, _, st_r), (_, _, ref_r) in zip(
            iterate_over_residues(result), iterate_over_residues(reference)
    ):
        # condition for NAs
        if len(ref_r.name) < 3:
            unter_na = st_r.name
            unter_na = unter_na.replace('3', '').replace('5', '')
            condition = unter_na == ref_r.name
        # condition for histidines
        elif ref_r.name.startswith('HI'):
            condition = st_r.name.startswith('HI') and ref_r.name.startswith('HI')
        else:
            condition = st_r.name == ref_r.name
        assert condition
        to_del = []
        for h_atom in st_r:
            if h_atom.is_hydrogen():
                if h_atom.name not in [atom.name for atom in ref_r]:
                    to_del.append(h_atom)
        for h_atom in to_del[::-1]:
            st_r.remove_atom(h_atom.name, ' ', gemmi.Element('H'))
        if st_r.name.startswith('HI'):
            st_r.name = ref_r.name
        if len(ref_r.name) < 3:
            st_r.name = ref_r.name
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


def iterate_over_residues(
    st: gemmi.Structure,
) -> typing.Iterator[typing.Tuple[gemmi.Model, gemmi.Chain, gemmi.Residue]]:
    for model in st:
        for chain in model:
            for residue in chain:
                yield model, chain, residue


def copy_residue_names(
    st: gemmi.Structure, reference: gemmi.Structure
) -> gemmi.Structure:
    result = st.clone()
    for (_, _, st_r), (_, _, ref_r) in zip(
        iterate_over_residues(result), iterate_over_residues(reference)
    ):
        if st_r.name != ref_r.name:
            # very crude check that there was no accidental shift so far
            assert st_r.name[:2] == ref_r.name[:2]
            st_r.name = ref_r.name
    return result


def check_chain_and_residue_numbering(st: gemmi.Structure, ref: gemmi.Structure, strict: bool = True) -> bool:
    result = True
    for (_, st_c, st_r), (_, ref_c, ref_r) in zip(
            iterate_over_residues(st), iterate_over_residues(ref)
    ):
        if strict:
            res_names_equal = st_r.name == ref_r.name
        else:
            res_names_equal = st_r.name[:2] == ref_r.name[:2]
        result = result and st_c.name == ref_c.name and st_r.seqid.num == ref_r.seqid.num and res_names_equal
    return result
