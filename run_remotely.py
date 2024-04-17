from pathlib import Path
import remote_runner
from remote_runner import Pool, LocalSlurmWorker, Task
from glob import glob
import pandas as pd


def main():
    remote_runner.log_to(".remote-runner.log")
    
    # path_to_csv = '/home/polina/xray_refinment/1_filter_pdb/biounit_polym_filtered_with_space_groups_monomers.csv'
    # pdb_df = pd.read_csv(path_to_csv, sep=',')
    # run_dirnames = pdb_df['pdb_id']
    run_dirnames = ['1nko_rapper']    
    paths_to_run_dirs = []

    for dirname in run_dirnames:
        paths_to_run_dirs += sorted(Path("./").glob(f"4_protocol_run/output/{dirname}/*/state.dill"))

    workers = [LocalSlurmWorker(
        remote_user_rc="""
export PYTHONPATH=$PYTHONPATH:$(pwd)
source /opt/amber22_with_patches/amber.sh
source /home/olebedenko/xray-refinment/venv_amber22_with_patches/bin/activate
source /home/polina/xray-refinement/venv/bin/activate
""",
        sbatch_args=[
            "--cpus-per-task=2",  # cpu here means threads
            "--partition=full",
            "--gpus=1",
            # "--nodelis=bionmr-mom-004",
        ]
    ) for _ in range(len(run_dirnames))]

    tasks = [
        Task.load(path)
        for path in paths_to_run_dirs
    ]
    print(len(tasks))
    Pool(workers).run(tasks)


if __name__ == "__main__":
    main()
