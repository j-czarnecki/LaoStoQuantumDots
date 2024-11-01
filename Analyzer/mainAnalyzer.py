from plotterClass import *
def main():

    plotter = Plotter(runsPath = "/net/ascratch/people/plgjczarnecki/LAO-STO-QD",
                            matchPattern = "RUN_.*", physicalParams = ("bz", ))
    #plotter.LoadSingleElectronPsi()
    plotter.LoadSingleElectronEnergies()
    plotter.LoadSingleElectronExpectations()
    plotter.SortData()
    #plotter.PlotSingleElectronPsi()
    plotter.PlotSingleElectronEnergies('s_z')
    return


if __name__ == "__main__":
    main()