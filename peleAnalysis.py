import os
import pandas as pd
import shutil

class PELEAnalyzer:

    def __init__(self):
        """
        Initialize object and assign atributes.
        """

        self.sampling = None
        self.protein_portion = None
        self.forcefield = None
        self.perturbation_protocol = None
        self.docking_tool = None
        self.experimental_data = None
        self.calculated_data = None

    def _PELEFoldersHierarchy(self):
        pass

    # First of all we need to define how are we going to download the data and with 
    # what format. Take a look at multistate library by Martin.