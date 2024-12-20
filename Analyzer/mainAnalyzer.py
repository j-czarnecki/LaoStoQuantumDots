from plotterClass import *


def main():

    plotter = Plotter(runsPath = "/net/ascratch/people/plgjczarnecki/LAO-STO-QD-edsr-2e",
                            matchPattern = "RUN_.*", physicalParams = ("bz", ))
    # plotter = Plotter(
    #     runsPath="/home/jczarnecki/LaoStoQdResults/LAO-STO-QD-testing-2e",
    #     matchPattern="RUN_.*",
    #     physicalParams=("bz",),
    # )
    # plotter.LoadSingleElectronPsi()
    # plotter.LoadSingleElectronEnergies()
    # plotter.LoadSingleElectronExpectations()
    # plotter.LoadMultiElectronEnergies()
    # plotter.LoadMultiElectronExpectations()
    # plotter.SortData()
    # plotter.PlotSingleElectronPsi()
    # plotter.PlotSingleElectronEnergies("s_z")
    # plotter.PlotMultiElectronEnergies("parity")
    plotter.PlotTimeMaxCoeffs()
    plotter.PlotTimeSingleMaxCoeffs()
    plotter.PrintNxmElems()
    plotter.PrintSingleNxmElems()
    return


if __name__ == "__main__":
    main()
