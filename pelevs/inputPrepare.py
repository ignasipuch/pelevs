import os
import shutil
import pandas as pd


class InputPreparation:
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

    def __init__(self, ligands, receptor=None):
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
        self.ligands_file = os.path.basename(ligands)
        self.receptor_format = None

        self._filesChecker(self.receptor, self.ligands)
        self._prepareFolders(self.receptor, self.ligands)

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

        if os.path.isdir('1_input_files/'):
            pass

        else:
            ligands_format = ligands.split('.')[-1]

            if ligands_format == 'csv':
                if pd.read_csv(ligands).shape[1] == 2:
                    print(' -     The ligands were passed in a csv with SMILES and id.')
                else:
                    raise Exception(
                        'FormatLigandsError: The format of the ligand file should be a csv with smiles and id.')
            elif ligands_format == 'pdb':
                print(' -     The ligand was passed in a pdb.')  
            else:
                raise Exception(
                    'FormatLigandsError: The format of the ligand file should be a csv with smiles and id.')

            if receptor == None:
                print(' -     No receptor was passed. Only ligand-related jobs will be possible.')

            else:
                receptor_file_format = receptor.split('.')[-1]

                if (receptor_file_format == 'pdb') or (receptor_file_format == 'sd') or (receptor_file_format == 'sdf') or (receptor_file_format == 'mol2'):
                    self.receptor_format = receptor_file_format
                    print(
                        ' -     The receptor is a(n) {} file.'.format(self.receptor_format))
                else:
                    raise Exception(
                        'FormatLigandsError: The format of the ligand files is not supported')

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

        # Creating main directories
        if not os.path.isdir('1_input_files'):
            os.mkdir('1_input_files')

        # Sorting ligands file
        if not os.path.isdir('1_input_files/ligands'):
            os.mkdir('1_input_files/ligands')

        if os.path.isfile(os.path.join('1_input_files/ligands', self.ligands_file)):
            print(' -     Ligand file is already in 1_input_files/ligands.')
        else:
            shutil.copy(ligands, os.path.join(
                '1_input_files/ligands', self.ligands_file))

        # Sorting receptor file
        if receptor is not None:
            if not os.path.isdir('1_input_files/receptor'):
                os.mkdir('1_input_files/receptor')

            if len(os.listdir('1_input_files/receptor')) != 0:
                pass
            else:
                shutil.copy(receptor, os.path.join(
                    '1_input_files/receptor',  os.path.basename(receptor)))
        
    def setUpLigPrepJob(self, ligands_output='ligands_out', pH=7., pH_tolerance=2., conformations=4):
        """
        Generate a job folder with necessary files to launch a ligprep job.

        Parameters
        ==========
        ligands_output : str
            Name of the ligprep output.
        pH : float
            pH at which the protonation will take place.
        pH_tolerance : float
            pH tolerance at which the protonation will take place.
        conformations : int
            Maximum number of conformations per inputted ligand.
        """

        ligands_input = self.ligands_file

        if not os.path.isdir('2_ligprep_job'):
            os.mkdir('2_ligprep_job')

        if not os.path.isdir('2_ligprep_job/job'):
            os.mkdir('2_ligprep_job/job')

        with open('2_ligprep_job/job/ligprep.sh', 'w') as filein:
            filein.writelines('$SCHRODINGER/ligprep -retain_i -ph {pH_val} -pht {pHt_val} -bff 16 -g -s {conf} -epik -icsv {ligands_in} -osd {ligands_out}.sdf'.format(
                pH_val=pH, pHt_val=pH_tolerance, conf=conformations, ligands_in=ligands_input, ligands_out=ligands_output))

        shutil.copy('1_input_files/ligands/{}'.format(ligands_input),
                    '2_ligprep_job/job')

    def setUpQMParametrization(self):
        """
        Generate a job folder with necessary files to launch a qm parametrization job.
        """

        def _ligandFormatChecker():
            """
            Check format of the ligand to begin preparation of the job.
            """

            possible_formats = ['pdb','mae']
            if self.ligands_file.split('.')[1] in possible_formats:
                pass
            else:
                raise Exception('LifandFormatError: This method only accepts either pdb or mae format for the ligand. Please check the format is correct.')

        def _generateDirectory():
            """
            Generate directories to store all the files.
            """

            ligand = self.ligands_file.split('.')[0]
            qm_path = '1_input_files/qm'
            qm_job_path = os.path.join(qm_path,ligand,'qm_job')
            ligands_folder_path = '1_input_files/ligands/'
            script_path = 'dockprotocol/scripts/qm.py'

            if not os.path.isdir(qm_path):
                os.mkdir(qm_path)            
            if not os.path.isdir(qm_job_path):
                os.makedirs(qm_job_path)

            print(' -     Setting up a qm parametrization job at: {}.'.format(qm_job_path))
            print(' -     Ligand to parametrize: {}.'.format(ligand))

            ligand_path = os.path.join(ligands_folder_path,self.ligands_file)

            shutil.copy(ligand_path, qm_job_path)
            shutil.copy(script_path, qm_job_path)

            return ligand

        def _generateRun(ligand):
            """
            Generate a run file to send the job in a local machine.
            """

            run_path = '1_input_files/qm/{}/qm_job/run.sh'.format(ligand)

            with open(run_path, 'w') as filein:
                filein.writelines('$SCHRODINGER/run python3 qm.py -f {}\n'.format(self.ligands_file))

            print(' -     Send job to local machine (cactus, bubbles, blossom) to send job.')

        _ligandFormatChecker()
        ligand = _generateDirectory()
        _generateRun(ligand)
