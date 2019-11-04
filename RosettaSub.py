"""File for overloading the subprocess function to run whatever Rosetta binary you wish to run."""

import subprocess, os, time


class RosettaSubprocess:

    def __init__(self, program, max_p, process_list):
        self.rosetta_directory = "/home/trinhlab/Documents/rosetta_bin_linux_2018.12.60119_bundle/main/source/bin/"
        self.program = program
        self.inputs = list()

        # Variables for the batch process
        self.process_list = process_list
        self.max_processes = max_p
        self.next_no = 0
        self.Processes = list()

        self.processes_complete = 0

    def start_new(self):
        self.inputs[1] = self.process_list[self.next_no]  # sets the "blank" variable to the appropriate file
        total_bin_dir = self.rosetta_directory + self.program
        args_array = [total_bin_dir] + self.inputs
        proc = subprocess.Popen(args_array)
        self.next_no += 1
        self.Processes.append(proc)

    def set_inputs(self, in_put):
        self.inputs = in_put

    def check_running(self):
        for p in range(len(self.Processes)-1, 0, -1):
            if self.Processes[p].poll() is not None:
                print("Process has finished.  Adding to completed queue.")
                self.processes_complete += 1
                del self.Processes[p]

        while (len(self.Processes) < self.max_processes) and (self.next_no < len(self.process_list)):
            self.start_new()

    # Calling this function initializes the calculation
    def run_batch(self,sleeptime=240):
        self.check_running()
        while len(self.Processes) > 0:
            time.sleep(sleeptime)
            self.check_running()
            if self.processes_complete + 1 == len(self.process_list):  # The last process is not deleted in self.Processes
                break


class RosettaSingleProcess:

    def __init__(self, program):
        self.rosetta_directory = "/home/trinhlab/Documents/rosetta_bin_linux_2018.12.60119_bundle/main/source/bin/"
        self.program = program
        self.inputs = list()

    def run_process(self):
        total_bin_dir = self.rosetta_directory + self.program
        args_array = [total_bin_dir] + self.inputs
        subprocess.run(args_array)  # use run for single processor machines, Popen and a pool for multi

    def set_inputs(self,in_put):
        self.inputs = in_put