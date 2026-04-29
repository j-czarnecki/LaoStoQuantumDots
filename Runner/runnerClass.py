import multiprocessing
import subprocess
import os
import sys

if os.path.exists("/net/home/pwojcik/.local/lib/python2.7/site-packages"):
  sys.path.insert(0, "/net/home/pwojcik/.local/lib/python2.7/site-packages")
import time
from runnerConfigClass import *

SCRATCH_PATH = os.getenv("SCRATCH")
HOME_PATH = os.getenv("HOME")


class Runner(RunnerConfig):
  def __init__(self):
    RunnerConfig.__init__(self)

  def run_slurm_param_value(self, paramValuePairs, runsDir: str, machine: str = "default"):
    """
    Runs jobs on ARES
    """
    newRunPath = self.__createRunDirStructure(runsDir, paramValuePairs)
    runnerCwd = os.getcwd()
    os.chdir(newRunPath)

    # Creating input.nml in given directory
    self.__createAndWriteInputNml(paramValuePairs)

    # setting up slurm script
    with open("job.sh", "w") as jobFile:
      print(self.jobHeader[machine], file=jobFile)
      print("cd " + newRunPath, file=jobFile)
      print(os.path.join(runnerCwd, "..", "bin", "lao_sto_qd.x"), file=jobFile)

    # queue slurm job
    simulate = subprocess.run(["sbatch", "job.sh"])
    os.chdir(runnerCwd)
    return newRunPath  # for sequential runner

  def __createRunDirStructure(self, runsDir: str, paramValuePairs: list[tuple[str, str, float | list[float]]]) -> str:
    pathToAppend = os.path.join(SCRATCH_PATH, runsDir)
    os.makedirs(pathToAppend, exist_ok=True)
    pathToAppend = os.path.join(pathToAppend, "RUN")

    for pair in paramValuePairs:
      if pair[0] != "self_consistency" and not isinstance(pair[2], list):
        pathToAppend = pathToAppend + f"_{pair[1]}_{pair[2]}"

    path = pathToAppend

    outputDir = f"OutputData"
    if not os.path.exists(path):
      os.mkdir(path)
      os.mkdir(os.path.join(path, outputDir))

    return path

  def __createAndWriteInputNml(self, paramValuePairs: list) -> None:
    # Getting namelist with parameters
    nml = self.LAO_STO_QD_default_nml()

    for pair in paramValuePairs:
      nml[pair[0]][pair[1]] = pair[2]  # editing all key-value pairs

    with open(os.path.join("OutputData", "quantum_dot.nml"), "w") as nmlFile:
      f90nml.write(nml, nmlFile, sort=False)
