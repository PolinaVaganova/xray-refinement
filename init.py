from __future__ import annotations
from amber_runner.MD import *
from typing import Optional, List, Union, Tuple, Iterator
import gemmi
from dataclasses import dataclass
import os
import re


@dataclass
class Config:
    draft_tleap_config_path: str
    production_tleap_config_path: str
    parmed_summary_in_path: str
    na_ion_pdb_path: str
    cl_ion_pdb_path: str
    water_pdb_path: str
    parmed_add_xray_parameters_in_config: str

    @classmethod
    def from_config_dir(cls, config_dir: Path):
        config = cls(
            parmed_summary_in_path=(config_dir / "parmed.summary.in").absolute().__str__(),
            draft_tleap_config_path=(config_dir / "draft.tleap.in.config").absolute().__str__(),
            production_tleap_config_path=(config_dir / "production.tleap.in.config").absolute().__str__(),
            na_ion_pdb_path=(config_dir / "na.pdb").absolute().__str__(),
            cl_ion_pdb_path=(config_dir / "cl.pdb").absolute().__str__(),
            water_pdb_path=(config_dir / "spce.pdb").absolute().__str__(),
            parmed_add_xray_parameters_in_config=(
                    config_dir / "parmed.add-xray-parameters.in.config").absolute().__str__()
        )
        return config


class PrepareHelper:
    @classmethod
    def estimate_weight(cls, st: gemmi.Structure):
        mass = 0.0
        aname_map = {"Na+": 22.989769, "Cl-": 35.453}
        element_map = {"P": 30.973762, "H": 1.00794, "C": 12.0107, "S": 32.065, "O": 15.999, "N": 14.0067}
        for m in st:
            for c in m:
                for r in c:
                    for a in r:
                        if a.name in aname_map:
                            mass += aname_map[a.name]
                        else:
                            mass += element_map[a.name[0]]

        return mass

    @classmethod
    def estimate_n_water_molecules(cls, st: gemmi.Structure):
        psv = 0.74
        MW = cls.estimate_weight(st)
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

    @classmethod
    def find_ss_bond_pairs(cls, st: gemmi.Structure) -> Iterator[Tuple[int, int]]:
        cyx_residues = []
        for chain in st[0]:
            for residue in chain:
                if residue.name == "CYX":
                    cyx_residues.append(residue)

        for i in range(len(cyx_residues)):
            ri = cyx_residues[i]
            sgi = ri["SG"][0]
            for j in range(i + 1, len(cyx_residues)):
                rj = cyx_residues[j]
                sgj = rj["SG"][0]
                if sgi.pos.dist(sgj.pos) < 2.3:
                    yield ri.seqid.num, rj.seqid.num

    @classmethod
    def write_sf_dat_file(cls, mtz: gemmi.Mtz,
                          output_filename: Union[os.PathLike, str]
                          ):
        """Write .tab file for pmemd.arx
        :param mtz: mtz file with P1 symmetry
        :param output_filename: output .dat filename
        """

        assert mtz.spacegroup.hm == "P 1"

        def missing_column(label):
            raise RuntimeError(f"MTZ file missing {label} column")

        R_FREE_FLAG = None
        r_free_pattern_string = r"(r.*free)|(free.*r)"
        r_free_pattern = re.compile(r_free_pattern_string, flags=re.IGNORECASE)

        for column in mtz.columns:  # type: gemmi.Mtz.Column
            if r_free_pattern.search(column.label):
                R_FREE_FLAG = column

        if R_FREE_FLAG is None:
            raise RuntimeError(f"MTZ file missing R-FREE-FLAG column (pattern: `{r_free_pattern_string}`)"
                               f"\nPresent columns: {[column.label for column in mtz.columns]}"
                               )

        H, K, L, FOBS, SIGMA_FOBS = [
            mtz.column_with_label(label) or missing_column(label)
            for label in ("H", "K", "L", "FOBS", "SIGFOBS")
        ]

        n_positive_r_flags = sum(R_FREE_FLAG)
        flag_is_one = n_positive_r_flags > len(R_FREE_FLAG) / 2

        with open(output_filename, "w") as output:
            output.write(f'{len(H)} 0\n')

            for h, k, l, fobs, sigma, r_flag in zip(H, K, L, FOBS, SIGMA_FOBS, R_FREE_FLAG):
                r = r_flag if flag_is_one else 1 - r_flag
                output.write(f"{h:3.0f} {k:3.0f} {l:3.0f} {fobs:15.8e} {sigma:15.8e} {r:1.0f}\n")

    @classmethod
    def get_b_factors(cls, st: gemmi.Structure):
        return [
            at.b_iso
            for res in cls.get_residues(st)
            for at in res
        ]

    @classmethod
    def get_residues(cls, st: gemmi.Structure):
        assert len(st) == 1, "Structure MUST have one MODEL"
        for chain in st[0]:
            for residue in chain:
                yield residue

    @classmethod
    def copy_b_factors(cls, src: gemmi.Structure, dst: gemmi.Structure):
        import warnings
        target_residues = list(cls.get_residues(dst))
        reference_residues = list(cls.get_residues(src))

        assert len(target_residues) % len(reference_residues) == 0
        prev_at = gemmi.Atom()
        prev_res = None
        prev_b_iso = 0
        for i in range(len(target_residues)):
            tgt_res = target_residues[i]
            ref_res = reference_residues[i % len(reference_residues)]
            if tgt_res.name != ref_res.name:
                warnings.warn(f"Residues don't match: {tgt_res} ~ {ref_res}")
            for at in tgt_res:  # type: gemmi.Atom
                try:
                    at.b_iso = ref_res[at.name][0].b_iso
                    prev_at = at
                    prev_res = ref_res
                    prev_b_iso = at.b_iso
                except RuntimeError as e:
                    warnings.warn(f"B-factor {prev_b_iso:.2f} assigned "
                                  f"to {tgt_res}.{at.name} "
                                  f"from {prev_res}.{prev_at.name} "
                                  f"({e})")
                    at.b_iso = prev_b_iso


class Prepare(Step, PrepareHelper):

    def __init__(self, name: str,
                 input_pdb_path: str,
                 input_mtz_path: str,
                 config: Config,
                 ):
        super().__init__(name)

        self.charge = None
        self.input_pdb_path = input_pdb_path
        self.input_mtz_path = input_mtz_path

        self.config: Config = config

        self.n_polymer_residues = None
        self.b_factors = None

    @property
    def wbox_prmtop(self):
        return str(self.step_dir / "wbox.prmtop")

    @property
    def wbox_tmp_prmtop_path(self):
        return str(self.step_dir / "wbox.tmp.prmtop")

    @property
    def wbox_xray_prmtop_path(self):
        return str(self.step_dir / "wbox.xray.prmtop")

    @property
    def structure_factors_dat(self):
        return str(self.step_dir / "sf.dat")

    @property
    def parmed_output(self):
        return str(self.step_dir / "parmed.out.txt")

    @property
    def parmed_add_xray_parameters_in(self):
        return str(self.step_dir / "parmed.add-xray-parameters.in")

    @property
    def input_pdb_protonated_path(self):
        return str(self.step_dir / "input.1.protonated.pdb")

    @property
    def input_pdb_reduced_path(self):
        return str(self.step_dir / "input.2.reduced.pdb")

    @property
    def input_pdb_p1_path(self):
        return str(self.step_dir / "input.3.p1.pdb")

    @property
    def input_pdb_neutral(self):
        return str(self.step_dir / "input.4.neutral.pdb")

    @property
    def production_tleap_input_pdb_path(self):
        return str(self.step_dir / "input.5.solvated.pdb")

    @property
    def draft_tleap_in_path(self):
        return str(self.step_dir / "draft.tleap.in")

    @property
    def production_tleap_in_path(self):
        return str(self.step_dir / "production.tleap.in")

    @property
    def wbox_dry_pdb_path(self):
        return str(self.step_dir / "wbox.dry.pdb")

    @property
    def wbox_pdb_path(self):
        return str(self.step_dir / "wbox.pdb")

    @property
    def wbox_inpcrd_path(self):
        return str(self.step_dir / "wbox.inpcrd")

    def check_occupancy(self, st: gemmi.Structure):
        import gemmi
        for model in st:
            for chain in model:
                for residue in chain:
                    for atom in residue:  # type: gemmi.Atom
                        if atom.occ < 1.0:
                            raise RuntimeError(f'Atom {atom} occupancy is less than one')

    def remove_water_and_alt_conformations(self) -> gemmi.Structure:
        import gemmi

        input_structure = gemmi.read_pdb(self.input_pdb_path, split_chain_on_ter=True)
        input_structure.remove_waters()
        input_structure.remove_alternative_conformations()

        return input_structure

    def protonate(self, input_structure: gemmi.Structure):
        import subprocess
        input_structure.write_pdb(self.input_pdb_protonated_path, numbered_ter=False, ter_ignores_type=True)
        subprocess.check_call([
            "conda", "run",
            "pdb4amber",
            "-d", "--reduce",
            "-i", self.input_pdb_protonated_path,
            '-o', self.input_pdb_reduced_path
        ])
        assert Path(self.input_pdb_reduced_path).exists()

    def build_unit_cell(self):
        import gemmi
        st = gemmi.read_pdb(self.input_pdb_reduced_path, split_chain_on_ter=True)
        g = gemmi.find_spacegroup_by_name(st.spacegroup_hm)

        def apply_to_chain(chain, op, cell):
            for residue in chain:
                for atom in residue:
                    fractional = cell.fractionalize(atom.pos)
                    transformed = gemmi.Fractional(*op.apply_to_xyz(fractional.tolist()))
                    cartesian = cell.orthogonalize(transformed)
                    atom.pos = cartesian

        identity = gemmi.Op("x,y,z")
        model = st[0]
        n_chains = len(model)
        for op in g.operations():
            if op != identity:
                for i in range(n_chains):
                    model.add_chain(model[i], pos=-1)
                    chain = model[-1]
                    apply_to_chain(chain, op, cell=st.cell)

        st.write_pdb(self.input_pdb_p1_path, numbered_ter=False, ter_ignores_type=True)

    def prepare_draft_simulation(self):
        import subprocess

        with open(self.config.draft_tleap_config_path) as config:
            tleap_in_config = config.read()
        tleap_in = tleap_in_config.format(
            step_dir=self.step_dir,
            input_pdb=self.input_pdb_p1_path
        )

        with open(self.draft_tleap_in_path, "w") as f:
            f.write(tleap_in)

        subprocess.check_call([
            'tleap', '-s',
            '-f', self.draft_tleap_in_path
        ])

        parmed_output = subprocess.check_output([
            'conda', 'run',
            'parmed',
            self.wbox_prmtop,
            self.config.parmed_summary_in_path,
        ]).decode('utf-8')

        with open(self.parmed_output, "w") as out:
            out.write(parmed_output)

        for line in parmed_output.split('\n'):
            if 'Total charge (e-)' in line:
                self.charge = int(float(line.split(':')[-1].strip().split()[0]))
                return

        raise RuntimeError("couldn't create a draft refinement unit cell")

    def prepare_b_factors(self):
        import gemmi
        input_pdb_protonated = gemmi.read_pdb(self.input_pdb_protonated_path, split_chain_on_ter=True)
        wbox_pdb = gemmi.read_pdb(self.wbox_dry_pdb_path, split_chain_on_ter=True)
        self.copy_b_factors(input_pdb_protonated, wbox_pdb)
        self.n_polymer_residues = len(list(self.get_residues(wbox_pdb)))
        wbox_pdb.write_pdb(self.wbox_dry_pdb_path, numbered_ter=False, ter_ignores_type=True)

    def prepare_structure_factors(self):
        mtz = gemmi.read_mtz_file(self.input_mtz_path)
        self.write_sf_dat_file(mtz=mtz,
                               output_filename=self.structure_factors_dat)

    def prepare_structure(self, input_structure: gemmi.Structure):
        import gemmi
        import shutil
        import subprocess
        st: gemmi.Structure = gemmi.read_pdb(self.wbox_dry_pdb_path, split_chain_on_ter=True)
        st.cell = input_structure.cell
        n_wat = self.estimate_n_water_molecules(st)

        n_ions = abs(self.charge)
        n_water_molecules = n_wat - 2 * n_ions

        st.spacegroup_hm = 'P 1'
        st.write_pdb(self.wbox_dry_pdb_path, numbered_ter=False, ter_ignores_type=True)

        n_protein_atoms = st[0].count_atom_sites()

        if self.charge == 0:
            shutil.copyfile(self.wbox_dry_pdb_path, self.input_pdb_neutral)
        else:
            ion_type = self.config.na_ion_pdb_path if self.charge < 0 else self.config.cl_ion_pdb_path
            subprocess.check_call([
                'AddToBox',
                '-c', self.wbox_dry_pdb_path,
                '-a', ion_type,
                '-na', str(n_ions),
                '-o', self.input_pdb_neutral,
                '-P', str(n_protein_atoms),
                '-RP', str(3.0),
                '-RW', str(6.0),
                '-G', str(0.2),
                '-V', str(1)
            ])
        subprocess.check_call([
            'AddToBox',
            '-c', self.input_pdb_neutral,
            '-a', self.config.water_pdb_path,
            '-na', str(n_water_molecules),
            '-o', self.production_tleap_input_pdb_path,
            '-P', str(n_protein_atoms + n_ions),
            '-RP', str(3.0),
            '-RW', str(3.0),
            '-G', str(0.2),
            '-V', str(1)
        ])

        ss_bonds = self.find_ss_bond_pairs(st)
        ss_bond_commands = "\n".join(f"bond wbox.{i}.SG wbox.{j}.SG" for i, j in ss_bonds)

        with open(self.config.production_tleap_config_path) as config:
            with open(self.production_tleap_in_path, "w") as rc:
                rc.write(config.read().format(
                    reference_pdb=self.production_tleap_input_pdb_path,
                    bonds_commands=ss_bond_commands,
                    step_dir=self.step_dir
                ))

        subprocess.check_call([
            'tleap', '-s',
            '-f', self.production_tleap_in_path
        ])

        subprocess.check_call([
            'ChBox',
            '-c', self.wbox_inpcrd_path,
            '-o', self.wbox_inpcrd_path,
            '-X', str(st.cell.a),
            '-Y', str(st.cell.b),
            '-Z', str(st.cell.c),
            '-al', str(st.cell.alpha),
            '-bt', str(st.cell.beta),
            '-gm', str(st.cell.gamma)
        ])

    def prepare_xray_prmtop(self):
        import subprocess
        with open(self.config.parmed_add_xray_parameters_in_config) as f:
            parmed_in = f.read().format(
                wbox_pdb_path=self.wbox_pdb_path,
                out_prmtop_path=self.wbox_tmp_prmtop_path
            )

        with open(self.parmed_add_xray_parameters_in, "w") as f:
            f.write(parmed_in)
        subprocess.check_call([
            'conda', 'run',
            'parmed',
            '--overwrite',
            self.wbox_prmtop,
            self.parmed_add_xray_parameters_in,
        ])

        subprocess.check_call([
            "add_xray",
            "-i", self.wbox_tmp_prmtop_path,
            "-o", self.wbox_xray_prmtop_path,
            "-scattering", "xray"
        ])

    def prepare_files_for_next_stages(self, md: RefinementProtocol):
        # Set global attributes
        md.sander.prmtop = self.wbox_xray_prmtop_path
        md.sander.inpcrd = self.wbox_inpcrd_path
        md.sander.refc = self.wbox_inpcrd_path

        # Configure minimize
        md.minimize.input.cntrl(
            imin=1, maxcyc=500, ncyc=150,
            ntb=1,
            ntr=0,
            cut=8.,
        )
        # Configure heating
        md.heating.input.cntrl(
            cut=8.,
            dt=0.002,
            imin=0,
            ioutfm=1,
            irest=0,
            nstlim=10000,
            ntb=1,
            ntc=2,
            ntf=2,
            ntpr=50,
            ntr=1,
            ntt=1,
            ntwr=1000,
            ntwx=200,
            ntx=1,
            temp0=293.0,
            tempi=0.0,
        )
        md.heating.input.pin(
            AmberInput.GroupSelection(
                title="Keep protein fixed with weak restraints",
                weight=10.0,
                residue_id_ranges=[(1, self.n_polymer_residues)]
            )
        )

        # Configure evolution
        md.evolution.input.cntrl(
            imin=0,
            irest=1, ntx=5, iwrap=1,
            ntb=1,
            ntt=3, gamma_ln=3., ig=-1,
            tempi=298.0, temp0=298.0,
            ntp=0, pres0=1.0, taup=2.0,
            cut=8.0,
            ntr=0,
            ntc=2, ntf=2,
            nstlim=5000, nscm=100, dt=0.002,
            ntpr=100, ntwx=100, ntwr=5000,
            ioutfm=1
        )

        md.evolution.input._get("xray")(
            spacegroup_name='P1',
            pdb_infile=self.wbox_dry_pdb_path,
            pdb_read_coordinates=False,
            reflection_infile=self.structure_factors_dat,
            xray_weight_initial=0.0,
            xray_weight_final=1.0,
            target='ml',
            bulk_solvent_model='afonine-2013'
        )

        # Configure cool
        md.cooling.input.cntrl(
            imin=0,
            ntx=5,
            irest=1,
            iwrap=1,
            nstlim=5000,
            dt=0.002,
            ntf=2,
            ntc=2,
            tempi=298.0,
            temp0=0.0,
            ntpr=100,
            ntwx=100,
            cut=8.0,
            ntb=1,
            ntp=0,
            ntt=3,
            gamma_ln=2.0,
            nscm=200,
            nmropt=1,
        )

        md.cooling.input._get("xray")(
            spacegroup_name='P1',
            pdb_infile=self.wbox_dry_pdb_path,
            pdb_read_coordinates=False,
            reflection_infile=self.structure_factors_dat,
            xray_weight_initial=1.0,
            xray_weight_final=1.0,
            target='ml',
            bulk_solvent_model='afonine-2013'
        )

        for istep1, istep2, value1, value2 in [
            (0, 500, 293.0, 262.5),
            (501, 625, 262.5, 262.5),
            (626, 1125, 262.5, 225.0),
            (1125, 1250, 225.0, 225.0),
            (1251, 1750, 225.0, 187.5),
            (1751, 1875, 187.5, 187.5),
            (1876, 2375, 187.5, 150.0),
            (2376, 2500, 150.0, 150.0),
            (2501, 3000, 150.0, 112.5),
            (3001, 3125, 112.5, 112.5),
            (3126, 3625, 112.5, 75.0),
            (3626, 3750, 75.0, 75.0),
            (3751, 4250, 75.0, 37.5),
            (4251, 4375, 37.5, 37.5),
            (4376, 4875, 37.5, 0.0),
            (4876, 5000, 0.0, 0.0),
        ]:
            md.cooling.input.varying_conditions.add(
                type='TEMP0', istep1=istep1, istep2=istep2,
                value1=value1, value2=value2
            )

    def run(self, md: RefinementProtocol):
        input_st = self.remove_water_and_alt_conformations()
        self.check_occupancy(input_st)
        self.protonate(input_st)
        self.build_unit_cell()
        self.prepare_draft_simulation()
        self.prepare_b_factors()
        self.prepare_structure_factors()
        self.prepare_structure(input_st)
        self.prepare_xray_prmtop()
        self.prepare_files_for_next_stages(md)


class RefinementProtocol(MdProtocol):
    def __init__(self):
        wd = Path("protocol_wd")
        self.mkdir(wd)
        MdProtocol.__init__(self, name="B0", wd=wd)

        self.sander = SanderCommand()

        self.sander = PmemdCommand()
        self.sander.executable = ["pmemd.cuda"]
        self.sander.allow_small_box = True

        config_dir = Path(__file__).parent / "data" / "config"
        global_config = Config.from_config_dir(config_dir)

        self.prepare = Prepare(
            name="prepare",
            input_pdb_path=(Path.cwd() / "2msi.pdb").absolute().__str__(),
            input_mtz_path=(Path.cwd() / "2msi.mtz").absolute().__str__(),
            config=global_config
        )

        self.minimize = SingleSanderCall("minimize")
        self.heating = SingleSanderCall("heating")
        self.evolution = SingleSanderCall("evolution")
        self.cooling = SingleSanderCall("cooling")


def main():
    md = RefinementProtocol()
    with ChangeDirectory(md.wd):
        md.run()


if __name__ == '__main__':
    main()
