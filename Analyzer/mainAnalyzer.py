from plotterClass import *


def main():

    # plotter = Plotter(runsPath = "/net/ascratch/people/plgjczarnecki/LAO-STO-QD-comparison",
    #                         matchPattern = "RUN_.*", physicalParams = ("bz", ))
    plotter = Plotter(
        runsPath="/home/jczarnecki/LaoStoQdResults/LAO-STO-QD-dense-B",
        matchPattern="RUN_.*",
        physicalParams=("bz",),
    )
    # plotter.LoadSingleElectronPsi()
    plotter.LoadSingleElectronEnergies()
    plotter.LoadSingleElectronExpectations()
    plotter.LoadMultiElectronEnergies()
    plotter.LoadMultiElectronExpectations()
    plotter.SortData()
    # plotter.PlotSingleElectronPsi()
    plotter.PlotSingleElectronEnergies("parity")
    plotter.PlotMultiElectronEnergies("parity")
    plotter.PlotTimeMaxCoeffs()
    plotter.PlotTimeSingleMaxCoeffs()
    return


if __name__ == "__main__":
    main()
