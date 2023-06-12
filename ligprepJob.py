import os
import shutil
import pandas as pd

class LigPrepJob:
    """
    Attributes
    ==========
    receptor : str
        Name of the receptor's file (sdf, mol2 or pdb).
    ligands : str
        Name of the file were ligands in a csv with SMILES are lcoated.
    receptor_format : str
        Receptor's format.

    Methods
    =======
    ligPrepJob(self, ligands_input, ligands_output, pH, pH_tolerance, conformations)
        Prepare a ligprep job with the inputted ligands.

    Hidden Methods
    ==============
    _filesChecker(self, receptor, ligands)
        Check if the formats inputted are correct.
    _prepareFolders(self, ligands)
        Prepare paths were inputs and outputs are going to be stored.
    """
    
    def __init__(self, receptor, ligands):
        """
        Check inputted formats and prepare folder system for further
        processing.

        Parameters
        ==========
        receptor : str
            Name of the file with the receptor
        ligands : str
            Name of the csv file with SMILES and id.
        """

        self.receptor = receptor
        self.ligands = ligands
        self.receptor_format = None
            
        self._filesChecker(receptor, ligands)
        self._prepareFolders(receptor, ligands)

    def _filesChecker(self, receptor, ligands):
        """
        Check validity of inputted files.

        Parameters
        ==========
        receptor : str
            Name of the file with the receptor
        ligands : str
            Name of the csv file with SMILES and id.
        """
        
        ligands_format = ligands.split('.')[-1]
            
        if ligands_format == 'csv':
            if pd.read_csv(ligands).shape[1] == 2:
                print(' -     The ligands were passed in a csv with SMILES and id.')
            else:
                raise Exception('FormatLigandsError: The format of the ligand file should be a csv with smiles and id.')
        else:
            raise Exception('FormatLigandsError: The format of the ligand file should be a csv with smiles and id.')
            
        receptor_file_format = receptor.split('.')[-1]
        
        if (receptor_file_format == 'pdb') or (receptor_file_format == 'sd') or (receptor_file_format == 'sdf') or (receptor_file_format == 'mol2'):
                self.receptor_format = receptor_file_format
                print(' -     The receptor is a(n) {} file.'.format(self.receptor_format))
        else:
            raise Exception('FormatLigandsError: The format of the ligand files is not supported')
          
    def _prepareFolders(self, receptor, ligands):
        """
        Generate folders and move files in their corresponding paths.

        Parameters
        ==========
        receptor : str
            Name of the file with the receptor
        ligands : str
            Name of the csv file with SMILES and id.
        """
        
        if not os.path.isdir('1_input_files'):
            os.mkdir('1_input_files')
            
        if not os.path.isdir('1_input_files/receptor'):
            os.mkdir('1_input_files/receptor')
        if not os.path.isdir('1_input_files/ligands'):
            os.mkdir('1_input_files/ligands')
            
        if not os.path.isdir('2_ligprep_job'):
            os.mkdir('2_ligprep_job')

        shutil.move(receptor, os.path.join('1_input_files/receptor', receptor))
        shutil.move(ligands, os.path.join('1_input_files/ligands', ligands))
             
    def ligPrepJob(self, ligands_input, ligands_output, pH=7., pH_tolerance=2., conformations=4):
        """
        Generate a job folder with necessary files to launch a ligprep job.

        Parameters
        ==========
        ligands_input : str
            Name of the csv file with ligands.
        ligands_output : str
            Name of the ligprep output.
        pH : float
            pH at which the protonation will take place.
        pH_tolerance : float
            pH tolerance at which the protonation will take place.
        conformations : int
            Maximum number of conformations per inputted ligand.
        """
        
        if not os.path.isdir('2_ligprep_job/job'):
            os.mkdir('2_ligprep_job/job')
            
        with open ('2_ligprep_job/job/ligprep.sh', 'w') as filein:
            filein.writelines('$SCHRODINGER/ligprep -retain_i -ph {pH_val} -pht {pHt_val} -bff 16 -g -s {conf} -epik -icsv {ligands_in} -osd {ligands_out}.sdf'.format(pH_val = pH, pHt_val = pH_tolerance, conf = conformations, ligands_in = ligands_input, ligands_out = ligands_output))
            
        shutil.copy('1_input_files/{}'.format(ligands_input), '2_ligprep_job/job')

