from pelevs.inputPrepare import InputPreparation
from pelevs.dockingJob import DockingJob
from pelevs.dockingAnalysis import DockingAnalyzer
from pelevs.peleJob import PELEJob
from pelevs.peleAnalysis import PELEAnalyzer

class pelevs(InputPreparation, DockingJob, DockingAnalyzer, PELEJob, PELEAnalyzer):
    def __init__(self):
        super().__init__()