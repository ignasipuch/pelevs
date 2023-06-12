import os
import shutil
import pandas as pd

class LigPrepJob:
    
    def __init__(self, receptor, ligands):
        self.receptor = receptor
        self.ligands = ligands
        self.receptor_format = None
            
        self._filesChecker(receptor, ligands)
        self._prepareFolders(receptor, ligands)

    def _filesChecker(self, receptor, ligands):
        
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
        
        if not os.path.isdir('2_ligprep_job/job'):
            os.mkdir('2_ligprep_job/job')
            
        with open ('2_ligprep_job/job/ligprep.sh', 'w') as filein:
            filein.writelines('$SCHRODINGER/ligprep -retain_i -ph {pH_val} -pht {pHt_val} -bff 16 -g -s {conf} -epik -icsv {ligands_in} -osd {ligands_out}.sdf'.format(pH_val = pH, pHt_val = pH_tolerance, conf = conformations, ligands_in = ligands_input, ligands_out = ligands_output))
            
        shutil.copy('1_input_files/{}'.format(ligands_input), '2_ligprep_job/job')

