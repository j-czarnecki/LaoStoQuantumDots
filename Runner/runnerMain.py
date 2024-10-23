from runnerClass import *
import re

def main():
  runner = Runner()
  nml_name = "external_parameters"
  param_name = "Bz"
  Bz_min = 0
  Bz_max = 40.0
  Bz_steps = 20
  dBz = abs(Bz_max - Bz_min) / Bz_steps
  Bz_table = [(nml_name, param_name, Bz_min + i * dBz) for i in range(Bz_steps + 1)]

  for Bz in Bz_table:
     runner.run_slurm_param_value([Bz], isAres=True)

if __name__ == "__main__":
    main()