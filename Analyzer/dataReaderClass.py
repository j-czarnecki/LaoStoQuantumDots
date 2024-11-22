import pandas as pd
import numpy as np
import os
import sys
import re

if os.path.exists("/net/home/pwojcik/.local/lib/python2.7/site-packages"):
    sys.path.insert(0, "/net/home/pwojcik/.local/lib/python2.7/site-packages")
import f90nml

class DataReader:
    def __init__(self, runsPath: str, matchPattern: str, physicalParams: tuple[str]):
        self.runsPath: str = runsPath
        self.matchPattern: str = matchPattern

        self.psi1: list = []
        self.psi2: list = []
        self.energies1: list = []
        self.energies2: list = []
        self.expectations1: list = []
        self.expectations2: list = []

        self.directories: list[str] = [dir for dir in os.listdir(self.runsPath) if re.match(self.matchPattern, dir)]
        self.params: list[list] = []

        self.__initParams(physicalParams)

    def __initParams(self, physicalParams):
        for dir in self.directories:
            nmlPath = os.path.join(self.runsPath, dir, 'OutputData', 'quantum_dot.nml')
            with open(nmlPath) as nmlFile:
                nml = f90nml.read(nmlFile)
                paramsValueList = []
                for param in physicalParams:
                    paramsValueList.append(nml['external_parameters'][param])
                self.params.append(tuple(paramsValueList))
        print(self.params)
    def SortData(self):
        sortedIndeces = sorted(range(len(self.params)), key = lambda x: self.params[x])
        #self.psi1 = [self.psi1[i] for i in sortedIndeces]
        #self.psi2 = [self.psi2[i] for i in sortedIndeces]
        self.energies1 = [self.energies1[i] for i in sortedIndeces]
        self.energies2 = [self.energies2[i] for i in sortedIndeces]
        self.expectations1 = [self.expectations1[i] for i in sortedIndeces]
        self.expectations2 = [self.expectations2[i] for i in sortedIndeces]
        self.params = sorted(self.params)

    def LoadSingleElectronPsi(self):
        for dir in self.directories:
            #print(dir)
            for n in range(1, 51): #TODO: This 50 should be deduced from .nml file!
                psiPath = os.path.join(self.runsPath, dir, 'OutputData', f'Psi_1_n{n}.dat')
                if os.path.exists(psiPath):
                    psi1 = pd.read_fwf(psiPath, skiprows = 1, infer_nrows = 100,
                                       names = ['kx', 'ky',
                                                'rePsi_xy_up', 'imPsi_xy_up', 'rePsi_xy_down', 'imPsi_xy_down',
                                                'rePsi_xz_up', 'imPsi_xz_up', 'rePsi_xz_down', 'imPsi_xz_down',
                                                'rePsi_yz_up', 'imPsi_yz_up', 'rePsi_yz_down', 'imPsi_yz_down'])
                    self.psi1.append(psi1)
                else:
                    print("File does not exists, skipping: ", psiPath)
                    continue
        return

    def LoadSingleElectronEnergies(self):
        for dir in self.directories:
            #print(dir)
            energiesPath = os.path.join(self.runsPath, dir, 'OutputData', f'Energies1.dat')
            if os.path.exists(energiesPath):
                energies = pd.read_fwf(energiesPath, skiprows = 1, infer_nrows = 100, names = ['state', 'E'])
                self.energies1.append(list(energies['E']))
            else:
                print("File does not exists, skipping: ", energiesPath)
                continue
        return
    def LoadSingleElectronExpectations(self):
        for dir in self.directories:
            #print(dir)
            expectationsPath = os.path.join(self.runsPath, dir, 'OutputData', f'Expectations_1.dat')
            if os.path.exists(expectationsPath):
                expectations = pd.read_fwf(expectationsPath, skiprows = 1, infer_nrows = 100, names = ['state', 's_x', 's_y', 's_z', 'd_xy', 'd_xz', 'd_yz'])
                self.expectations1.append(expectations)
            else:
                print("File does not exists, skipping: ", expectationsPath)
                continue
        return

    def LoadMultiElectronEnergies(self):
        for dir in self.directories:
            #print(dir)
            energiesPath = os.path.join(self.runsPath, dir, 'OutputData', f'Energies2.dat')
            if os.path.exists(energiesPath):
                energies = pd.read_fwf(energiesPath, skiprows = 1, infer_nrows = 100, names = ['state', 'E'])
                self.energies2.append(list(energies['E']))
            else:
                print("File does not exists, skipping: ", energiesPath)
                continue
        return

    def LoadMultiElectronExpectations(self):
        for dir in self.directories:
            print(dir)
            expectationsPath = os.path.join(self.runsPath, dir, 'OutputData', f'Expectations_2.dat')
            if os.path.exists(expectationsPath):
                expectations = pd.read_fwf(expectationsPath, skiprows = 1, infer_nrows = 100, names = ['state', 'x', 's_x', 's_y', 's_z'])
                self.expectations2.append(expectations)
            else:
                print("File does not exists, skipping: ", expectationsPath)
                continue
        return

    def LoadTimeDependence(self, dir):
        header = ['t']
        for i in range(1, 21):
            header.append(f'c_{i}')

        timeDependencePath = os.path.join(self.runsPath, dir, 'OutputData', f'Time_dependent.dat')
        if os.path.exists(timeDependencePath):
            timeDependence = pd.read_fwf(timeDependencePath, skiprows = 0, infer_nrows = 100, names = header)
            return timeDependence
        else:
            print("File does not exists, skipping: ", timeDependencePath)
            return

    def LoadMaxCoeffs(self, dir):
        header = ['omega_ac']
        for i in range(1, 21):
            header.append(f'c_{i}')

        cMaxPath = os.path.join(self.runsPath, dir, 'OutputData', f'C_max_time.dat')
        if os.path.exists(cMaxPath):
            cMax = pd.read_fwf(cMaxPath, skiprows = 1, infer_nrows = 100, names = header)
            cMax.sort_values('omega_ac', inplace = True)
            return cMax
        else:
            print("File does not exists, skipping: ", cMaxPath)
            return

    def LoadSingleMaxCoeffs(self, dir):
        header = ['omega_ac']
        for i in range(1, 51):
            header.append(f'c_{i}')

        cMaxPath = os.path.join(self.runsPath, dir, 'OutputData', f'C_single_max_time.dat')
        if os.path.exists(cMaxPath):
            cMax = pd.read_fwf(cMaxPath, skiprows = 1, infer_nrows = 100, names = header)
            cMax.sort_values('omega_ac', inplace = True)
            return cMax
        else:
            print("File does not exists, skipping: ", cMaxPath)
            return
