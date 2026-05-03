import f90nml
import textwrap


class RunnerConfig:
  def __init__(self):
    # Using dedent to remove indentation and align every row with the first column of file
    self.jobHeader = {
      "default": textwrap.dedent(
        """\
                #!/bin/bash
                ##### Amount of cores per task
                #SBATCH --cpus-per-task=16
                ##### Partition name
                #SBATCH -p cpu
                ##### Name of job in queuing system
                #SBATCH --job-name=DOS
                #SBATCH --output=\"output.out\"    # Path to the standard output and error files relative to the working directory
                """
      ),
      "ares": textwrap.dedent(
        """\
                #!/bin/bash -l
                ## Job name
                #SBATCH -J v_tuning_NNN
                ## Number of allocated nodes
                #SBATCH -N 1
                ## Number of tasks per node (by default this corresponds to the number of cores allocated per node)
                #SBATCH --ntasks-per-node=48
                ## Memory allocated per core (default is 5GB), comment if mem for whole job should be taken
                ##SBATCH --mem-per-cpu=3800MB
                ## Memory allocated for whole job, comment if mem-per-cpu should be taken
                #SBATCH --mem=16GB
                ## Max task execution time (format is HH:MM:SS)
                #SBATCH --time=72:00:00
                ## Name of grant to which resource usage will be charged
                #SBATCH -A plgktosto111-cpu
                ## Name of partition
                #SBATCH -p plgrid
                ## Name of file to which standard output will be redirected
                #SBATCH --output="output.out"
                ## Name of file to which the standard error stream will be redirected
                #SBATCH --error="error.err"
                """
      ),
      "helios": textwrap.dedent(
        """\
                #!/bin/bash -l
                ## Job name
                #SBATCH -J v_tuning_NNN
                ## Number of allocated nodes
                #SBATCH -N 1
                ## Number of tasks per node (by default this corresponds to the number of cores allocated per node)
                #SBATCH --ntasks-per-node=48
                ## Memory allocated per core (default is 5GB), comment if mem for whole job should be taken
                ##SBATCH --mem-per-cpu=3800MB
                ## Memory allocated for whole job, comment if mem-per-cpu should be taken
                #SBATCH --mem=16GB
                ## Max task execution time (format is HH:MM:SS)
                #SBATCH --time=72:00:00
                ## Name of grant to which resource usage will be charged
                #SBATCH -A plgktosto111-cpu
                ## Name of partition
                #SBATCH -p plgrid
                ## Name of file to which standard output will be redirected
                #SBATCH --output="output.out"
                ## Name of file to which the standard error stream will be redirected
                #SBATCH --error="error.err"
                module load GCC/13.2.0 OpenMPI/5.0.3 FlexiBLAS/3.3.1 ScaLAPACK/2.2.0-fb gimkl/2023b
                """
      ),
    }

  def LAO_STO_QD_default_nml(self):
    parser = f90nml.Parser()
    params_nml = parser.reads(
      f"&calculation_parameters \
                Nx=80, \
                Ny=80, \
                dx=.39, \
                norbs=6, \
                nstate_1=22, \
                nstate_2 = 20, \
                k_electrons = 2, \
                dt = 5e-6, \
                t_max = 5./ \
              &physical_parameters \
                th=0.04, \
                tl=0.875, \
                td=0.04, \
                dso=0.01, \
                drso=0.02, \
                dE=0.047, \
                g=3.0, \
                eps_r = 100.0/ \
              &external_parameters \
                omega = 18.689e-3, \
                Bx=0.0, \
                By=0.0, \
                Bz=0.0, \
                domega_ac = 2e-6, \
                omega_ac_max = 30e-3, \
                f_ac = 1e6, \
                Vb = 0.0, \
                V0 = 0.0, \
                d_image = 5.0/ \
              &self_consistency \
                max_sc_iter = 30, \
                eps_potential = 1e-6, \
                sc_alpha = 0.2/"
    )
    return params_nml
