import multiprocessing
import subprocess
import os
import sys

if os.path.exists("/net/home/pwojcik/.local/lib/python2.7/site-packages"):
    sys.path.insert(0, "/net/home/pwojcik/.local/lib/python2.7/site-packages")
import time
from runnerConfigClass import *


class Runner(RunnerConfig):
    def __init__(self):
        RunnerConfig.__init__(self)

    def run_slurm_param_value(self, paramValuePairs, isAres: bool = False):
        """
        Runs jobs on ARES
        """
        if isAres:
            pathToAppend = f"/net/ascratch/people/plgjczarnecki/LAO-STO-QD/RUN"
        else:
            pathToAppend = f"RUN"

        for pair in paramValuePairs:
            if pair[0] != "calculation_parameters":
                pathToAppend = pathToAppend + f"_{pair[1]}_{pair[2]}"

        runner_cwd = os.getcwd()
        if isAres:
            path = pathToAppend
        else:
            path = os.path.join(runner_cwd, pathToAppend)

        output_dir = f"OutputData"
        if not os.path.exists(path):
            os.mkdir(path)
            os.mkdir(os.path.join(path, output_dir))
        os.chdir(path)

        nml = self.LAO_STO_QD_default_nml()  # creating default namelist
        for pair in paramValuePairs:
            nml[pair[0]][pair[1]] = pair[2]  # editing all key-value pairs

        with open("./OutputData/quantum_dot.nml", "w") as nml_file:
            f90nml.write(nml, nml_file, sort=False)
        # setting up slurm script
        with open("job.sh", "w") as job_file:
            if isAres:
                print(self.job_header_ares, file=job_file)
            else:
                print(self.job_header, file=job_file)

            print("cd " + path, file=job_file)
            print(os.path.join(runner_cwd, "..", "bin", "lao_sto_qd.x"), file=job_file)

        # queue slurm job
        simulate = subprocess.run(["sbatch", "job.sh"])
        os.chdir(runner_cwd)
        return path  # for sequential runner