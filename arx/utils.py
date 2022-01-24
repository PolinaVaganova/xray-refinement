import os
import subprocess
import sys
from contextlib import contextmanager
from pathlib import Path


@contextmanager
def chdir(path: Path):
    """Sets the cwd within the context

    Args:
        path (Path): The path to the cwd

    Yields:
        None
    """

    origin = Path().absolute()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(origin)


def check_call(*args, **kwargs) -> None:
    """
    A substitute for subprocess.check_call that suppresses output
    unless an error occurs

    :param args: Same as in subprocess.Popen
    :param kwargs: Same as in subprocess.Popen
    :return: None
    """
    if "stderr" not in kwargs:
        stderr_capture = subprocess.PIPE
        kwargs["stderr"] = stderr_capture
    else:
        stderr_capture = None

    if "stdout" not in kwargs:
        kwargs["stdout"] = subprocess.DEVNULL

    with subprocess.Popen(*args, **kwargs) as p:
        try:
            return_code = p.wait()
        except:  # noqa: E722
            p.kill()
            if stderr_capture:
                if p.stderr:
                    sys.stderr.write(p.stderr.read())
            raise
    cmd = kwargs.get("args") or args[0]
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)
