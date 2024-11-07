import f90nml
import textwrap


class RunnerConfig:
    def __init__(self):

        # Using dedent to remove indentation and align every row with the first column of file
        self.job_header = textwrap.dedent(
            """\
        #!/bin/bash
        #SBATCH --job-name=LAO_STO_QD              # Job name
        #SBATCH --partition tera-cpu       # we specify to run the process on gpu nodes
        #SBATCH --ntasks-per-node=1        # Maximum number of tasks on each node
        #SBATCH --time=72:00:00            # Wall time limit (days-hrs:min:sec)
        #SBATCH --mem-per-cpu=3800MB         # Memory (i.e. RAM) per processor
        #SBATCH --output=\"output.out\"    # Path to the standard output and error files relative to the working directory
        """
        )

        self.job_header_ares = textwrap.dedent(
            """\
        #!/bin/bash -l
        ## Job name
        #SBATCH -J LAO-STO-QD
        ## Number of allocated nodes
        #SBATCH -N 1
        ## Number of tasks per node (by default this corresponds to the number of cores allocated per node)
        #SBATCH --ntasks-per-node=1
        ## Memory allocated per core (default is 5GB)
        #SBATCH --mem-per-cpu=3800MB
        ## Max task execution time (format is HH:MM:SS)
        #SBATCH --time=168:00:00
        ## Name of grant to which resource usage will be charged
        #SBATCH -A plglaosto111-cpu
        ## Name of partition
        #SBATCH -p plgrid-long
        ## Name of file to which standard output will be redirected
        #SBATCH --output="output.out"
        ## Name of file to which the standard error stream will be redirected
        #SBATCH --error="error.err"
        """
        )

    def LAO_STO_QD_default_nml(self):
        parser = f90nml.Parser()
        params_nml = parser.reads(
            f"&calculation_parameters \
                Nx=80, \
                Ny=80, \
                dx=.39, \
                norbs=6, \
                nstate_1=50, \
                nstate_2 = 20, \
                k_electrons = 2, \
                dt = 1e-5, \
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
                domega_ac = 1e-4, \
                omega_ac_max = 30e-3, \
                f_ac = 1e6/"
        )
        return params_nml
