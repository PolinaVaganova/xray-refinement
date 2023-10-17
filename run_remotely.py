from pathlib import Path

import remote_runner
from remote_runner import LocalSlurmWorker, Pool, Task


def main():
    remote_runner.log_to(".remote-runner.log")

    workers = [
        LocalSlurmWorker(
            remote_user_rc="""
unset PYTHONPATH
<see run.sh for an example>
""",
            sbatch_args=[
                "--cpus-per-task=1",  # cpu here means threads
                "--gpus=1",
            ],
        )
        for _ in range(4)
    ]

    tasks = [
        Task.load(state)
        for state in sorted(Path("./").glob("data/output/*/*/state.dill"))
    ]
    print(len(tasks))
    Pool(workers).run(tasks)


if __name__ == "__main__":
    main()
