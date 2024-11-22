from dataReaderClass import *
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.collections import LineCollection
from matplotlib import cm
from scipy.interpolate import interp1d
from matplotlib.colors import LinearSegmentedColormap
import numpy as np


class Plotter(DataReader):
    def __init__(self, runsPath: str, matchPattern: str, physicalParams: tuple[str]):
        DataReader.__init__(self, runsPath, matchPattern, physicalParams)
        self.__initRcParams()

    def __initRcParams(self):
        plt.rcParams["text.usetex"] = True
        plt.rcParams["font.family"] = "serif"
        plt.rcParams["font.serif"] = "Computer Modern Roman"
        plt.rcParams["font.sans-serif"] = "Computer Modern Sans serif"
        plt.rcParams["font.monospace"] = "Computer Modern Typewriter"
        plt.rcParams["axes.titlesize"] = 24
        plt.rcParams["axes.labelsize"] = 24
        plt.rcParams["xtick.labelsize"] = 20
        plt.rcParams["ytick.labelsize"] = 20
        # Optionally, add custom LaTeX preamble
        plt.rcParams["text.latex.preamble"] = (
            r"\usepackage{amsmath} \usepackage{amsfonts} \usepackage{amssymb}"
        )

        # Choose a seaborn palette
        palette = sns.color_palette("hsv", 7)  # has to specify number of lines

        # Set the color cycle
        plt.rcParams["axes.prop_cycle"] = plt.cycler(color=palette)

        # Set rcParams for tighter layout
        plt.rcParams["figure.autolayout"] = True
        plt.rcParams["figure.constrained_layout.use"] = True
        plt.rcParams["axes.linewidth"] = 1.2

        # Set rcParams to show ticks on both left and right sides
        plt.rcParams["xtick.direction"] = "in"
        plt.rcParams["ytick.direction"] = "in"
        plt.rcParams["xtick.bottom"] = True
        plt.rcParams["ytick.left"] = True
        plt.rcParams["xtick.top"] = True
        plt.rcParams["ytick.right"] = True

        plt.rcParams["legend.fontsize"] = 12
        plt.rcParams["legend.title_fontsize"] = 14

        plt.rcParams["axes.xmargin"] = 0.01

    def PlotSingleElectronPsi(self):
        plt.figure()
        plt.scatter(self.psi1[5]['kx'],
                    self.psi1[5]['ky'],
                    c = self.psi1[5]['rePsi_xy_up']**2 + self.psi1[5]['imPsi_xy_up']**2,
                    cmap = 'inferno')
        plt.savefig('../Plots/Psi.png')
        plt.close()

    def PlotSingleElectronEnergies(self, colorParam: str):
        plt.figure()
        fig, ax = plt.subplots(figsize = (6,10))

        Bz, *_ = zip(*self.params)
        color = [df[colorParam].tolist() for df in self.expectations1]
        for state in range(0,50):
            energy = [sublist[state] for sublist in self.energies1]
            colorState = [sublist[state] for sublist in color]

            #This smoothens the colorlines
            num_points = 500
            interp_Bz = np.linspace(min(Bz), max(Bz), num_points)
            interp_energy = interp1d(Bz, energy, kind='linear')(interp_Bz)
            interp_color = interp1d(Bz, colorState, kind='linear')(interp_Bz)

            #Create segments for colored lines
            points = np.array([interp_Bz, interp_energy]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)

            # Create LineCollection with color mapping
            lc = LineCollection(segments, array=interp_color, cmap='coolwarm', norm=plt.Normalize(-1, 1))
            lc.set_linewidth(2)
            ax.add_collection(lc)

        ax.autoscale()
        plt.colorbar(lc, ax=ax, label=r"$\langle s_z \rangle$")
        plt.xlabel(r'$B_z$ (T)')
        plt.ylabel(r'$E$ (meV)')
        plt.locator_params(nbins=3, axis='x')
        plt.locator_params(nbins=10, axis='y')
        fig.tight_layout(pad = 3.0)
        plt.savefig('../Plots/Energies1.png', dpi = 300)
        plt.close()

    def PlotMultiElectronEnergies(self, colorParam: str):
        spins_cmap = LinearSegmentedColormap.from_list("spins_cmap", ["deepskyblue", "black", "deeppink"])
        plt.figure()
        fig, ax = plt.subplots(figsize = (6,10))

        Bz, *_ = zip(*self.params)
        color = [df[colorParam].tolist() for df in self.expectations2]
        for state in range(0,20):
            energy = [sublist[state] for sublist in self.energies2]
            colorState = [sublist[state] for sublist in color]

            #This smoothens the colorlines
            num_points = 500
            interp_Bz = np.linspace(min(Bz), max(Bz), num_points)
            interp_energy = interp1d(Bz, energy, kind='linear')(interp_Bz)
            interp_color = interp1d(Bz, colorState, kind='linear')(interp_Bz)

            #Create segments for colored lines
            points = np.array([interp_Bz, interp_energy]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)

            # Create LineCollection with color mapping
            lc = LineCollection(segments, array=interp_color, cmap=spins_cmap, norm=plt.Normalize(-2, 2))
            lc.set_linewidth(2)
            ax.add_collection(lc)


        ax.autoscale()
        plt.colorbar(lc, ax=ax, label=r"$\langle S_z \rangle$")
        plt.xlabel(r'$B_z$ (T)')
        plt.ylabel(r'$E$ (meV)')
        plt.locator_params(nbins=3, axis='x')
        plt.locator_params(nbins=10, axis='y')
        fig.tight_layout(pad = 3.0)
        plt.savefig('../Plots/Energies2.png', dpi = 300)
        plt.close()


    def PlotTimeDependence(self):
        df = self.LoadTimeDependence('RUN_Bz_0.0')
        print(df.iloc[0])
        print(np.sum(df.iloc[0]))
        # Choose a seaborn palette
        palette = sns.color_palette("hsv", 6)  # has to specify number of lines

        # Set the color cycle
        plt.rcParams["axes.prop_cycle"] = plt.cycler(color=palette)

        for i in range(1,6):
            plt.plot(df['t'], df[f'c_{i}'], label = fr'$|c_{i}|^2$')
        plt.legend()
        plt.xlabel(r'$t$ (ns)')
        plt.ylabel(r'$|c_n|^2$')
        plt.savefig('../Plots/TimeDependence.png', dpi = 300)

    def PlotTimeMaxCoeffs(self):
        df = self.LoadMaxCoeffs('RUN_Bz_12')
        print(df.iloc[0])
        print(np.sum(df.iloc[0]))
        # Choose a seaborn palette
        palette = sns.color_palette("hsv", 9)  # has to specify number of lines

        # Set the color cycle
        plt.rcParams["axes.prop_cycle"] = plt.cycler(color=palette)

        for i in range(2, 9):
            plt.plot(df['omega_ac'], df[f'c_{i}'], label = fr'$|c_{i}|^2$')
        plt.title(r'$B_z = 10.0$ (T)')
        plt.legend(loc = 'upper right')
        plt.xlabel(r'$\hbar \omega_{AC}$ (meV)')
        plt.ylabel(r'$max(|c_n|^2(t))$')
        plt.savefig('../Plots/CMax.png', dpi = 300)

    def PlotTimeSingleMaxCoeffs(self):
        df = self.LoadSingleMaxCoeffs('RUN_Bz_10')
        print(df.iloc[0])
        print(np.sum(df.iloc[0]))
        # Choose a seaborn palette
        palette = sns.color_palette("hsv", 9)  # has to specify number of lines

        # Set the color cycle
        plt.rcParams["axes.prop_cycle"] = plt.cycler(color=palette)

        for i in range(2, 9):
            plt.plot(df['omega_ac'], df[f'c_{i}'], label = fr'$|c_{i}|^2$')
        plt.title(r'$B_z = 10.0$ (T)')
        plt.legend(loc = 'upper right')
        plt.xlabel(r'$\hbar \omega_{AC}$ (meV)')
        plt.ylabel(r'$max(|c_n|^2(t))$')
        plt.savefig('../Plots/CMax_single.png', dpi = 300)

