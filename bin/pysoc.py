#!/usr/bin/env python3.6

"""
Driver script for pysoc.
"""

import subprocess

def run(command):
    """
    Run an external program, watching for certain types of errors.

    :param command: An iterable of the command to be run (eg, ["ls", "./", "-l"])
    """
    try:
        subprocess.run(
            command,
            check = True
            )
    except FileNotFoundError as e:
        # The specified command could not be found.
        raise Exception("Could not find external command '{}'; is PySOC setup correctly?".format(command[0])) from e

# First we run the python script which performs setup (?).
run(("soc.py",))

# Then the fortran program.
run(("soc_td",))


