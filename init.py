from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import List, Union

import gemmi
from amber_runner.MD import MdProtocol, PmemdCommand, SingleSanderCall, Step
from remote_runner import Task

from arx.utils import check_call


def count_polymer_residues(
    st: gemmi.Structure,
) -> int:
    non_polymer_residue_names = ["WAT", "Cl-", "Na+"]
    count = 0
    for model in st:
        for chain in model:
            for residue in chain:  # type: gemmi.Residue
                if residue.name in non_polymer_residue_names:
                    return count
                else:
                    count += 1
    return count


class Prepare(Step):
    @classmethod
    def write_sf_dat_file(
        cls, mtz: "gemmi.Mtz", output_filename: Union[os.PathLike, str]
    ):
        """Write .tab file for pmemd.arx
        :param mtz: mtz file with P1 symmetry
        :param output_filename: output .dat filename
        """
        import re

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
            raise RuntimeError(
                f"MTZ file missing R-FREE-FLAG column "
                f"(pattern: `{r_free_pattern_string}`)"
                f"\nPresent columns: {[column.label for column in mtz.columns]}"
            )

        H, K, L, FOBS, SIGMA_FOBS = [
            mtz.column_with_label(label) or missing_column(label)
            for label in ("H", "K", "L", "FOBS", "SIGFOBS")
        ]

        n_positive_r_flags = sum(R_FREE_FLAG)
        flag_is_one = n_positive_r_flags > len(R_FREE_FLAG) / 2

        with open(output_filename, "w") as output:
            output.write(f"{len(H)} 0\n")

            for h, k, l, fobs, sigma, r_flag in zip(
                H, K, L, FOBS, SIGMA_FOBS, R_FREE_FLAG
            ):
                r = r_flag if flag_is_one else 1 - r_flag
                output.write(
                    f"{h:3.0f} {k:3.0f} {l:3.0f} {fobs:15.8e} {sigma:15.8e} {r:1.0f}\n"
                )

    def __init__(
        self,
        name: str,
        parm7: Path,
        rst7: Path,
        pdb: Path,
        mtz: Path,
    ):
        Step.__init__(self, name)

        self.charge = None
        self.input_pdb_path: Path = pdb
        self.input_mtz_path: Path = mtz
        self.wbox_prmtop = parm7
        self.wbox_inpcrd_path = rst7
        self.wbox_pdb = pdb

    @property
    def wbox_xray_prmtop_path(self):
        return self.step_dir / "wbox.xray.parm7"

    @property
    def structure_factors_dat(self):
        return str(self.step_dir / "sf.dat")

    def prepare_structure_factors(self):
        import gemmi

        mtz = gemmi.read_mtz_file(str(self.input_mtz_path))
        self.write_sf_dat_file(mtz=mtz, output_filename=self.structure_factors_dat)

    def prepare_xray_prmtop(self):
        tmp_parm = self.step_dir / "tmp.parm7"
        parmed_add_xray_parameters_in = self.step_dir / "parmed.add-xray-parameters.in"
        parmed_in = f"""
addPdb {self.wbox_pdb} elem strict allicodes
lmod
parmout {tmp_parm}
go
"""

        with open(parmed_add_xray_parameters_in, "w") as f:
            f.write(parmed_in)
        check_call(
            [
                "parmed",
                "--overwrite",
                str(self.wbox_prmtop),
                str(parmed_add_xray_parameters_in),
            ]
        )

        check_call(
            [
                "add_xray",
                "-i",
                str(tmp_parm),
                "-o",
                str(self.wbox_xray_prmtop_path),
                "-scattering",
                "xray",
            ]
        )

    def prepare_files_for_next_stages(self, md: RefinementProtocol):
        from amber_runner.inputs import AmberInput

        # Set global attributes
        md.sander.prmtop = str(self.wbox_xray_prmtop_path)
        md.sander.inpcrd = str(self.wbox_inpcrd_path)
        md.sander.refc = str(self.wbox_inpcrd_path)

        n_polymer_residues = count_polymer_residues(gemmi.read_pdb(str(self.wbox_pdb)))

        # Configure minimize
        md.minimize.input.cntrl(
            imin=1,
            maxcyc=500,
            ncyc=150,
            ntb=1,
            ntr=0,
            cut=8.0,
        )
        # Configure heating
        md.heating.input.cntrl(
            cut=8.0,
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
                residue_id_ranges=[(1, n_polymer_residues)],
            )
        )

        # Configure evolution
        md.evolution.input.cntrl(
            imin=0,
            irest=1,
            ntx=5,
            iwrap=1,
            ntb=1,
            ntt=3,
            gamma_ln=3.0,
            ig=-1,
            tempi=298.0,
            temp0=298.0,
            ntp=0,
            pres0=1.0,
            taup=2.0,
            cut=8.0,
            ntr=0,
            ntc=2,
            ntf=2,
            nstlim=5000,
            nscm=100,
            dt=0.002,
            ntpr=100,
            ntwx=100,
            ntwr=5000,
            ioutfm=1,
        )

        md.evolution.input._get("xray")(
            spacegroup_name="P1",
            pdb_infile=str(self.wbox_pdb),
            pdb_read_coordinates=False,
            reflection_infile=self.structure_factors_dat,
            atom_selection_mask=f":1-{n_polymer_residues}",
            xray_weight_initial=0.0,
            xray_weight_final=1.0,
            target="ml",
            bulk_solvent_model="afonine-2013",
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
            spacegroup_name="P1",
            pdb_infile=str(self.wbox_pdb),
            pdb_read_coordinates=False,
            reflection_infile=self.structure_factors_dat,
            atom_selection_mask=f":1-{n_polymer_residues}",
            xray_weight_initial=1.0,
            xray_weight_final=1.0,
            target="ml",
            bulk_solvent_model="afonine-2013",
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
                type="TEMP0", istep1=istep1, istep2=istep2, value1=value1, value2=value2
            )

    def run(self, md: RefinementProtocol):
        self.prepare_structure_factors()
        self.prepare_xray_prmtop()
        self.prepare_files_for_next_stages(md)


class ConvertToPdb(Step):
    def run(self, md: "RefinementProtocol"):
        import tempfile

        from arx.prepare import copy_coordinates, read_pdb, write_pdb

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_pdb = Path(tmp_dir) / "tmp.pdb"
            # Don't change directory to avoid problems
            # with relative paths in md.sander.*
            with open(tmp_pdb, "wb") as f:
                check_call(
                    [
                        "ambpdb",
                        # ambpdb doesn't work well with xray-modified prmtop
                        # use original parm7
                        "-p",
                        str(md.prepare.wbox_prmtop),
                        "-c",
                        str(md.sander.inpcrd),
                    ],
                    stdout=f,
                )
            final = read_pdb(tmp_pdb)

        initial = read_pdb(md.prepare.wbox_pdb)
        result = copy_coordinates(initial, reference=final)
        write_pdb(result, self.step_dir / "final.pdb")


class RefinementProtocol(MdProtocol):
    def __init__(self, pdb: Path, mtz: Path, parm7: Path, rst7: Path, output_dir: Path):
        wd = output_dir
        wd.mkdir(mode=0o755, exist_ok=True, parents=True)
        MdProtocol.__init__(self, name="B0", wd=wd)

        self.sander = PmemdCommand()
        self.sander.executable = ["pmemd.cuda"]
        self.sander.allow_small_box = True

        self.prepare = Prepare(
            name="prepare",
            pdb=pdb,
            mtz=mtz,
            parm7=parm7,
            rst7=rst7,
        )

        self.minimize = SingleSanderCall("minimize")
        self.heating = SingleSanderCall("heating")
        self.evolution = SingleSanderCall("evolution")
        self.cooling = SingleSanderCall("cooling")
        self.convert_to_pdb = ConvertToPdb("convert_to_pdb")


def create_tasks(subset: str):
    tasks = []
    for pdb_code in ["2msi"]:
        input_dir = Path.cwd() / "data" / "input" / pdb_code / subset
        output_dir = Path.cwd() / "data" / "output" / pdb_code / subset

        # Copy input files to trajectory folder
        input_copy = output_dir / "inputs"
        input_copy.mkdir(exist_ok=True, parents=True)
        pdb = input_copy / "wbox.pdb"
        mtz = input_copy / f"{pdb_code}.mtz"
        rst7 = input_copy / "wbox.rst7"
        parm7 = input_copy / "wbox.parm7"

        shutil.copy(input_dir / "wbox.pdb", pdb)
        shutil.copy(input_dir / f"{pdb_code}.mtz", mtz)
        shutil.copy(input_dir / "wbox.rst7", rst7)
        shutil.copy(input_dir / "wbox.parm7", parm7)

        md = RefinementProtocol(
            pdb=pdb.relative_to(output_dir),
            mtz=mtz.relative_to(output_dir),
            rst7=rst7.relative_to(output_dir),
            parm7=parm7.relative_to(output_dir),
            output_dir=output_dir,
        )
        md.save(md.wd / "state.dill")
        tasks.append(md)
    return tasks


def run_all_locally(tasks: List[Task]):
    import remote_runner
    from remote_runner import LocalWorker, Pool

    remote_runner.log_to(".remote-runner.log", level="DEBUG")
    workers = [LocalWorker()]
    Pool(workers).run(tasks)


def run_sequentially_inplace(tasks: List[Task]):
    from remote_runner.utility import ChangeDirectory

    for md in tasks:
        with ChangeDirectory(md.wd):
            md.run()


def main():
    tasks = create_tasks("prepared")
    assert tasks
    run_sequentially_inplace(tasks)


if __name__ == "__main__":
    main()
