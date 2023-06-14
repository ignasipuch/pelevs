from openbabel import pybel
import os
import shutil
import pandas as pd

class DockingAnalyzer:

    def __init__(self):
        
        self.receptor = os.listdir('1_input_files/receptor')[0]
        self.ligands = os.listdir('1_input_files/ligands')[0]
        self.docking_tool = None
        self.experimental_data = None

    def _glideDockingResultsChecker(self):
        
        path_docking = '3_docking_job/job'
        path_results = os.path.join(path_docking,[x for x in os.listdir(path_docking) if x.endswith('.csv')][0])

        if not os.path.isfile(path_results):
            raise Exception('ResultsMissingError: Before initializing the object the results must be downloaded and at {}'.format(path_docking))
        
    def _glideDataFrameRetriever(self):
        # From the whole dataframe keep the best conformations with their score.
        pass

    def _rdockDockingResultsChecker(self):

        path_docking = '3_docking_job/job/results'
        path_results = [x for x in os.listdir(path_docking) if x.endswith('.sdf')][0]

        if len(path_results) == 0 and path_results[0].endswith('sd'):
            raise Exception('ResultsMissingError: Before initializing the object the results must be downloaded and located at {}'.format(path_docking))
        
    def _rdockDataFrameGenerator(self):
        # Generate a dataframe with the name of the molecle and the score
        pass

    def glideAnalysis(self):
        
        self._glideDockingResultsChecker()

    def rdockAnalysis(self):

        self._rdockDockingResultsChecker()

