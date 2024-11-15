from plotterClass import *
def main():

    plotter = Plotter(runsPath = "/net/ascratch/people/plgjczarnecki/LAO-STO-QD",
                            matchPattern = "RUN_.*", physicalParams = ("bz", ))
    #plotter.LoadSingleElectronPsi()
    plotter.LoadSingleElectronEnergies()
    plotter.LoadSingleElectronExpectations()
    plotter.LoadMultiElectronEnergies()
    plotter.LoadMultiElectronExpectations()
    plotter.SortData()
    #plotter.PlotSingleElectronPsi()
    plotter.PlotSingleElectronEnergies('s_z')
    plotter.PlotMultiElectronEnergies('s_z')
    plotter.PlotTimeMaxCoeffs()
    return


if __name__ == "__main__":
    main()