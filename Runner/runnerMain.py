from runnerClass import *
import re


def main():
    runner = Runner()
    nml_name = "external_parameters"
    param_name = "Bz"
    Bz_min = 0
    Bz_max = 12.0
    Bz_steps = 12
    dBz = abs(Bz_max - Bz_min) / Bz_steps
    Bz_table = [(nml_name, param_name, Bz_min + i * dBz) for i in range(Bz_steps + 1)]
    #Bz_table = [(nml_name, param_name, B) for B in [0.0, 10.0, 12.0]]

    # for Bz in Bz_table:
    #     runner.run_slurm_param_value([Bz], isAres=True)

    runner.run_slurm_param_value([(nml_name, param_name, 1.0)], isAres=True)
    runner.run_slurm_param_value([(nml_name, param_name, 10.0)], isAres=True)
    runner.run_slurm_param_value([(nml_name, param_name, 12.0)], isAres=True)



if __name__ == "__main__":
    main()
