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
        
        if not os.path.isdir('3_docking_job/job/runs'):
            os.mkdir('3_docking_job/job/runs')       

        with open('3_docking_job/job/runs/grid.sh', 'w') as fileout:
            fileout.writelines(
                'rbcavity -was -d -r parameter_file.prm > parameter_file.log\n'
            )   
            
    def _rdockJobSplitter(self, ligands, cpus_docking):
        
        # Generating split file 
        if not os.path.isfile('3_docking_job/job/splitMols.sh'):
            with open ('3_docking_job/job/splitMols.sh', 'w') as filein:
                filein.writelines(
                '#!/bin/bash\n'
                '#Usage: splitMols.sh <input> #Nfiles <outputRoot>\n'
                'module load rdock\n'
                'fname=$1\n'
                'nfiles=$2\n'
                'output=$3\n'
                'molnum=$(grep -c \'$$$$\' $fname)\n'
                'echo " - $molnum molecules found"\n'
                'echo " - Dividing \'$fname\' into $nfiles files"\n'
                'rawstep=`echo $molnum/$nfiles | bc -l`\n'
                'let step=$molnum/$nfiles\n'
                'if [ ! `echo $rawstep%1 | bc` == 0 ]; then\n'
                '        let step=$step+1;\n'
                'fi;\n'
                'sdsplit -$step -o$output $1\n'    
                )

                
        # Generating splitted ligand files    
        with open('3_docking_job/job/runs/split.sh', 'w') as fileout:
            fileout.writelines(
                '../splitMols.sh ../{ligands_file} {cpus} ../ligands\n'.format(ligands_file=ligands, cpus=cpus_docking)
            ) 
            
    def _rdockRunFilesGenerator(self, cpus_docking):
        
        # Generating folders 
        if not os.path.isdir('3_docking_job/job/ligands'):
            os.mkdir('3_docking_job/job/ligands')
            
        if not os.path.isdir('3_docking_job/job/results'):
            os.mkdir('3_docking_job/job/results')
    
        # Generating run files
        for i in range(1,cpus_docking+1):
            with open('3_docking_job/job/runs/run{}'.format(i),'w') as fileout:
                fileout.writelines([
                '#!/bin/sh\n',
                '#SBATCH --job-name=rdock' + str(i) + ' \n',
                '#SBATCH --time=01:00:00\n',
                '#SBATCH --ntasks=1\n',
                '#SBATCH --output=rdock.out\n',
                '#SBATCH --error=rdock.err\n',
                '\n',
                'module load rdock\n',
                'module load ANACONDA/2019.10\n',
                'module load intel\n',
                'module load mkl\n',
                'module load impi\n',
                'module load gcc\n',
                'module load boost/1.64.0\n',
                '\n',
                '\n',
                'rbdock -i ../ligands/split{val}.sd -o ../results/split{val}_out -r ../parameter_file.prm -p dock.prm -n 50\n'.format(val=i),
                'done\n'
                ])
                
        with open('3_docking_job/job/run_prepare_rDock.sh','w') as fileout:
            fileout.writelines(
                '#!/bin/bash\n'
                '# Run runs/run_grid.sh\n'
                'source runs/run_grid.sh\n'
                '\n'
                '# Run runs/run_split.sh\n'
                'source runs/run_split.sh\n'       
            )
            
        with open('3_docking_job/job/run_rDock.sh','w') as fileout:
            fileout.writelines(
                '#!/bin/bash\n'
                'for d in run*; do echo ${d}; sbatch ${d}; done'   
            )
            
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
                          
    def setRdockDocking(self, reference_ligand, ligands, cpus_docking):
                   
        self.reference_ligand = reference_ligand
        self.docking_tool = 'rdock'
        
        self._rdockReceptorFormatChecker(self.receptor)
        
        shutil.copy(reference_ligand, '3_docking_job/job')
        shutil.copy('2_ligprep_job/job/' + self.ligands, '3_docking_job/job')
        shutil.copy('1_input_files/receptor/' + self.receptor, '3_docking_job/job')
        
        self._paramFilesWriter(self.receptor, self.reference_ligand)
        self._rdockGridGenerator()
        self._rdockJobSplitter(ligands,cpus_docking)
        self._rdockRunFilesGenerator(cpus_docking)