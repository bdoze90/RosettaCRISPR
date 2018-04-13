"""File for overloading the subprocess function to run whatever Rosetta binary you wish to run."""

import subprocess, os


class RosettaSubprocess:

    def __init__(self, program):
        self.rosetta_directory = "/Users/brianmendoza/rosetta_bin_mac_2018.12.60119_bundle/main/source/bin/"
        self.program = program
        self.inputs = list()

    def run_process(self):
        total_bin_dir = self.rosetta_directory + self.program
        args_array = [total_bin_dir] + self.inputs
        subprocess.run(args_array)  # use run for single processor machines, Popen and a pool for multi

    def set_inputs(self,in_put):
        self.inputs = in_put

