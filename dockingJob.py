from openbabel import pybel
import os
import shutil

class DockingJob:
   
    def __init__(self, receptor='1_input_files/receptor/' + os.listdir('1_input_files/receptor/')[0], ligands='2_ligprep_job/job/' + [x for x in os.listdir('2_ligprep_job/job/') if x.endswith('.sdf')][0]):
      
        self.receptor = receptor.split('/')[-1]
        self.ligands = ligands.split('/')[-1]
        self.docking_tool = None
        self.grid = None
        self.reference_ligand = None
        
        self._ligandsChecker(ligands)
        self._folderPreparation()
        
    def _ligandsChecker(self, ligands):
        
        if (ligands.split('.')[-1] != 'sdf') or (ligands.split('.')[-1] != 'sd'):
            pass
        else: 
            raise Exception('LigandsFileError: The ligands file shoud be the output from ligprep in sdf format.')
    
    def _folderPreparation(self):
        
        if not os.path.isdir('3_docking_job'):
            os.mkdir('3_docking_job')
            
        if not os.path.isdir('3_docking_job/job'):
            os.mkdir('3_docking_job/job')

    def _rdockReceptorFormatChecker(self, receptor):
        
        if len(receptor.split('.')) == 2:
            receptor_name, receptor_format = receptor.split('.')
            
        else:
            raise Exception('ReceptorNameError: Receptor name should only have one dot separating the name and the format.')

        if receptor_format != 'mol2':
            receptor_generator = pybel.readfile(receptor_format, '1_input_files/receptor/' + receptor)
            receptor_molecule = next(receptor_generator)  # Convert generator to molecule object
            receptor_molecule.write("mol2", "3_docking_job/job/{}.mol2".format(receptor_name), overwrite=True)
            
        else: pass
        
    def _paramFilesWriter(self, receptor, reference_ligand):
            
        parameter_file = os.path.join('3_docking_job/job', 'parameter_file.prm')
            
        if not os.path.isfile(parameter_file):
            with open(parameter_file, 'w') as fileout:
                fileout.writelines(
                'RBT_PARAMETER_FILE_V1.00\n'
                'TITLE rdock\n'
                '\n'
                'RECEPTOR_FILE ' + receptor + '\n'
                'RECEPTOR_FLEX 3.0\n'
                '\n'
                '##################################################################\n'
                '### CAVITY DEFINITION: REFERENCE LIGAND METHOD\n'
                '##################################################################\n'
                'SECTION MAPPER\n'
                '    SITE_MAPPER RbtLigandSiteMapper\n'
                '    REF_MOL ' + reference_ligand + '\n'
                '    RADIUS 6.0\n'
                '    SMALL_SPHERE 1.0\n'
                '    MIN_VOLUME 100\n'
                '    MAX_CAVITIES 1\n'
                '    VOL_INCR 0.0\n'
                '   GRIDSTEP 0.5\n'
                'END_SECTION\n'
                '\n'
                '#################################\n'
                '#CAVITY RESTRAINT PENALTY\n'
                '#################################\n'
                'SECTION CAVITY\n'
                '    SCORING_FUNCTION RbtCavityGridSF\n'
                '    WEIGHT 1.0\n'
                'END_SECTION\n')
                
    def _rdockGridGenerator(self):

        with open('3_docking_job/job/run_grid.sh', 'w') as fileout:
            fileout.writelines(
                'rbcavity -was -d -r parameter_file.prm > parameter_file.log\n'
            )   
            
    def _rdockJobSplitter(self):
        pass
           
    def setGlideDocking(self, grid_file, forcefield='OPLS_2005'):
        
        self.grid_file = grid_file
        self.docking_tool = 'glide'
        
        shutil.copy(grid_file, '3_docking_job/job')
        shutil.copy('2_ligprep_job/job/' + self.ligands, '3_docking_job/job')
            
        with open('3_docking_job/job/glide_job.in', 'w') as filein:
            filein.writelines([
                'FORCEFIELD   {}\n'.format(forcefield),
                'GRIDFILE   {}\n'.format(grid_file),
                'LIGANDFILE   {}\n'.format(self.ligands),
                'PRECISION   SP\n'
            ])
            
        print(' - Glide job generated succesfully.')
                             
    def setRdockDocking(self, reference_ligand):
                   
        self.reference_ligand = reference_ligand
        self.docking_tool = 'rdock'
        
        self._rdockReceptorFormatChecker(self.receptor)
        
        shutil.copy(reference_ligand, '3_docking_job/job')
        shutil.copy('2_ligprep_job/job/' + self.ligands, '3_docking_job/job')
        shutil.copy('1_input_files/receptor/' + self.receptor, '3_docking_job/job')
        
        self._paramFilesWriter(self.receptor, self.reference_ligand)
        self._rdockGridGenerator()