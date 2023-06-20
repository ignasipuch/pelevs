import os 
import shutil
from openbabel import openbabel as ob
from Bio.PDB import PDBParser, PDBIO
import re
from openbabel import pybel

class PELE:
    """
    Attributes
    ==========
    receptor : str
        Name of the receptor's file (sdf, mol2 or pdb).

    Methods
    =======
    glideAnalysis(self, experimental_data, column_name)
        Calculate energetic correlation between Glide's predictions and
        experimental data. Also it plots the distribution of time spent
        per ligand.

    Hidden Methods
    ==============
    _correlationPlotter(self, x, y, docking_method)
        Plotts x and y in a z_score format and stores 
        the image.
    """

    def __init__(self):
        """
        Initialize object and assign atributes.
        """

        self._folderPreparation()

    def _folderPreparation(self):
        
        if not os.path.isdir('4_pele_simulation'):
            os.mkdir('4_pele_simulation')

        if not os.path.isdir('4_pele_simulation/docking_input'):
            os.mkdir('4_pele_simulation/docking_input')

        if not os.path.isdir('4_pele_simulation/pele_simulation'):
            os.mkdir('4_pele_simulation/pele_simulation')

    def _folderHierarchy(self, force_field, truncated, perturbation_protocol, rescoring_method):

        if force_field is None:
            force_field = 'opls'
            forcefield_list = [False, force_field]
        else: 
            print(' - Forcefield chosen : {}.'.format(force_field))
            forcefield_list = [True, force_field]

        if truncated is None:
            truncated = 'truncated'
            truncated_list = [False, truncated]
        else: 
            print(' - Protein section chosen : {}.'.format(truncated))
            truncated_list = [True, truncated]

        if perturbation_protocol is None:
            perturbation_protocol = 'refinement'
            perturbation_list = [False, perturbation_protocol]
        else: 
            print(' - Perturbation protocol chosen : {}.'.format(perturbation_protocol))
            perturbation_list = [True, perturbation_protocol]

        print(' - Sampling chosen : {}.'.format(rescoring_method))
        rescoring_method_list = [True, rescoring_method]

        return forcefield_list, truncated_list, perturbation_list, rescoring_method_list

    def _PELEJobManager(self, forcefield_list, truncated_list, perturbation_protocol_list, rescoring_method_list):
        
        path = '4_pele_simulation/'
        pele_simulation_path = '4_pele_simulation/pele_simulation'

        lists = [forcefield_list, truncated_list, perturbation_protocol_list, rescoring_method_list]
        directory = ''
        cont = 1

        for lst in lists:
            if lst[0] is True:
                if cont == 1:
                    directory_path = os.path.join(path, lst[1])
                    pele_directory_path = os.path.join(pele_simulation_path, lst[1])

                directory_name = lst[1]
                new_directory = os.path.join(path, directory_name)

                if not os.path.isdir(new_directory):
                    os.makedirs(new_directory, exist_ok=True)

                directory = os.path.join(directory,directory_name)
                path = new_directory

                cont += 1
        
        # Move all to pele_simulation_path
        if not os.path.isdir(directory_path):
            shutil.move(directory_path, pele_simulation_path)
        elif not os.path.isdir(os.path.join(pele_simulation_path, directory)):
            shutil.move(path, pele_directory_path)
            shutil.rmtree(directory_path)

        return os.path.join(pele_simulation_path, directory)

    def _PELEJobChecker(self, forcefield_list, truncated_list, perturbation_protocol_list, rescoring_method_list):

        path = '4_pele_simulation/pele_simulation'

        directories = []
        directories.append(path)

        lists = [forcefield_list, truncated_list, perturbation_protocol_list, rescoring_method_list]
        options = [inner_list[1] for inner_list in lists if inner_list[0]]

        path_to_create = '/'.join(options)
        path_simulation = os.path.join(path, path_to_create)

        previous_simulation_bool = False

        # Create path for new job
        if not os.path.isdir(path_simulation):
            os.mkdir(path_simulation)

        # Generate list of directories to check
        for lst in lists:
            if lst[0] is True:
                path = os.path.join(path,lst[1])
                directories.append(path)

        # Checks if there are other created jobs obtains path_to_retrieve
        for directory in directories[:-1]:
            content_in_directory = [x for x in os.listdir(directory) if not x.startswith('.')]
            if len(content_in_directory) > 1:
                previous_simulation_bool = True
                if len(content_in_directory) < 5:
                    previously_generated_simulation = [x for x in content_in_directory if x not in options][0]
                    path_to_retrieve = os.path.join(directory, previously_generated_simulation)
                else: 
                    path_to_retrieve = directory
            else: path_to_retrieve = None

        # Copies files if there are files to copy
        if path_to_retrieve is not None:
            for ligand in os.listdir(path_to_retrieve):
                file = os.path.join(path_to_retrieve, ligand, ligand + '.pdb')
                path_to_ligand = os.path.join(path_simulation,ligand)
                if not os.path.isdir(path_to_ligand):
                    os.mkdir(path_to_ligand)

                shutil.copy(file, path_to_ligand)

        return previous_simulation_bool

    def _PDBConversor(self, file_in, path_out):
        """
        Converts the input file to PDB format using Open Babel.
        """

        file_name, file_format = file_in.split('.')
        file_name = os.path.basename(file_name)  

        if not file_format == 'pdb':

            conv = ob.OBConversion()
            conv.SetInAndOutFormats(file_in.split('.')[-1], 'pdb')
            mol = ob.OBMol()
            conv.ReadFile(mol, file_in)
            conv.WriteFile(mol, os.path.join(path_out,'{}.pdb'.format(file_name)))

    def _PDBMerger(self, receptor, ligand):
        
        def _receptorModifier(receptor):

            file_mod_prot = receptor.split('.pdb')[0] + '_mod.pdb'

            # Assigning a random residue number.
            parser = PDBParser()
            io = PDBIO()

            structure_prot = parser.get_structure('prot', receptor)

            index = 1000

            for res in structure_prot.get_residues():
                res.id = (' ',index,' ')
                index += 1

            io.set_structure(structure_prot)
            io.save(file_mod_prot)

            # Counting number of atoms of the protein
            ligand_cont_num = 0

            with open (file_mod_prot) as filein:
                for line in filein:
                    sline = line.split()
                    if sline[0] == 'ATOM':
                        ligand_cont_num += 1

            return file_mod_prot, ligand_cont_num

        def _ligandAtomChainNumberModifier(ligand):

            file_mod1_lig = ligand.split('.pdb')[0] + '_m.pdb'
            file_mod_lig = ligand.split('.pdb')[0] + '_mod.pdb'

            cont = 1
            with open(ligand) as filein:
                with open(file_mod1_lig, 'w') as fileout:
                    for line in filein:
                        if line.startswith('HETATM'):
                            sline = line.split()
                            beginning_of_line = line[0:13]
                            end_of_line = line[17:]
                            atom = ''.join([c for c in sline[2] if c.isalpha()])
                            new_line = beginning_of_line + (atom + str(cont)).ljust(4) + end_of_line
                            fileout.writelines(new_line)
                            cont += 1   

            # Change the chain name and number
            with open(file_mod1_lig, 'r') as filein:
                with open(file_mod_lig, 'w') as fileout:
                    for line in filein:                      
                        if re.search('UN.......', line):
                            line = re.sub('UN.......','LIG L 900', line)
                            fileout.writelines(line)
                        else:
                            fileout.writelines(line)

            return file_mod_lig  

        def _receptorLigandMerger(receptor, ligand):

            output_file = 'intermediate.pdb'
            path_files = os.path.dirname(ligand)
            writing_path = os.path.join(path_files, output_file)

            os.system('obabel {receptor} {ligand} -O {output}'.format(receptor=receptor,ligand=ligand,output=writing_path))

            return writing_path

        def _PELEInputAdapter(merged, output_file, ligand_cont_num):

            # Changing the atom numeration of the ligand.
            with open (merged, 'r') as filein:
                with open(output_file, 'w') as fileout:
                    for line in filein:
                        sline = line.split()
                        if sline[0] == 'ATOM' or sline[0] == 'HETATM':
                            if re.search('HETATM.....', line):
                                # Considering length of the protein
                                if len(str(ligand_cont_num + 1)) < 5:
                                    line = re.sub('HETATM.....','HETATM ' + str(ligand_cont_num + 1), line) 
                                elif len(str(ligand_cont_num + 1)) == 5:
                                    line = re.sub('HETATM.....','HETATM' + str(ligand_cont_num + 1), line) 

                                fileout.writelines(line)
                                ligand_cont_num += 1
                            else:
                                fileout.writelines(line)

        def _intermediateFilesRemover(output_file):

            working_directory = os.path.dirname(output_file)
            input_PELE_file = os.path.basename(output_file)

            for file in os.listdir(working_directory):
                if file != input_PELE_file:
                    os.remove(os.path.join(working_directory,file))

        ligand_name = os.path.basename(ligand)
        working_directory = os.path.dirname(ligand)
        output_file = os.path.join(working_directory,ligand_name)

        file_mod_prot, ligand_cont_num = _receptorModifier(receptor)
        file_mod_lig = _ligandAtomChainNumberModifier(ligand)
        intermediate = _receptorLigandMerger(file_mod_prot,file_mod_lig)
        _PELEInputAdapter(intermediate, output_file, ligand_cont_num)
        _intermediateFilesRemover(output_file)

    def _PELESimulationFiles(self, path, pdb_file, force_field, truncated, perturbation_protocol, rescoring_method):
  
        with open(os.path.join(path, 'input.yaml'), 'w') as fileout:
            fileout.writelines(
                'complex_data: "' + pdb_file + '"\n'
                'complex_ligand_selection: "L:900"\n'
                'static_name: false\n'
                'name: "' + pdb_file.split('.')[0] + '"\n'
                'verbosity: info\n')
            
        
            if truncated == 'truncated':
                fileout.writelines(
                    'flexible_region_radius: 8\n'
                    'frozen_region_radius: 13\n')

            elif truncated == 'full':
                pass

            else:
                fileout.writelines(
                    'flexible_region_radius: 8\n'
                    'frozen_region_radius: 13\n')

            if perturbation_protocol == 'if':
                fileout.writelines('pele_perturbation_level: 3\n')
            elif perturbation_protocol == 'refinement':
                fileout.writelines('pele_perturbation_level: 2\n')
            elif perturbation_protocol == 'minimization':
                fileout.writelines('pele_perturbation_level: 0\n')
            else:
                raise Exception('InvalidProtocol: You have to choose either if or refinement.')

            if force_field == 'opls':
                fileout.writelines('pele_forcefield: opls2005\n')
            elif force_field == 'openff':
                fileout.writelines('pele_forcefield: openff-2.0.0\n')
            else:
                raise Exception('InvalidForceField: You have to choose either opls or openff.')

            if rescoring_method == 'xlong':
                fileout.writelines(
                    'pele_steps: 25\n'
                    'adaptive_epochs: 5\n'
                    'cpus: 48\n')
            elif rescoring_method == 'short':
                fileout.writelines(
                    'pele_steps: 10\n'
                    'adaptive_epochs: 1\n'
                    'cpus: 16\n')
            elif rescoring_method == 'long':
                fileout.writelines(
                    'pele_steps: 10\n'
                    'adaptive_epochs: 5\n'
                    'cpus: 32\n')
            elif rescoring_method == 'xshort':
                fileout.writelines(
                    'pele_steps: 5\n'
                    'adaptive_epochs: 1\n'
                    'cpus: 8\n')
            else:
                raise Exception('RescoringMethodError: The rescoring method entered is not valid.\
                     Try: xshort, short, long, or xlong.')

            fileout.writelines(
                'pipeline:\n'
                '   - flow: induced_fit_refinement\n'
                )

        with open(os.path.join(path, 'run_plat'), 'w') as fileout:
            fileout.writelines(
                '#!/bin/bash\n'
                '#SBATCH -J ' + pdb_file.split('.pdb')[0] + '\n'
                '#SBATCH --output=PELE.out\n'
                '#SBATCH --error=PELE.err\n')

            if rescoring_method == 'xlong':
                fileout.writelines(
                    '#SBATCH --ntasks=48\n'
                )
            elif rescoring_method == 'short':
                fileout.writelines(
                    '#SBATCH --ntasks=16\n'
                )
            elif rescoring_method == 'long':
                fileout.writelines(
                    '#SBATCH --ntasks=32\n'
                )
            elif rescoring_method == 'xshort':
                fileout.writelines(
                    '#SBATCH --ntasks=8\n'
                )
            else:
                raise Exception('RescoringMethodError: The rescoring method entered is not valid.\
                     Try: rescoring, short or long.')

            fileout.writelines(
                '#SBATCH --time=00-06:00:00\n'
                '\n'
                'module load ANACONDA/2019.10\n'
                'module load intel mkl impi gcc # 2> /dev/null\n'
                'module load impi\n'
                'module load boost/1.64.0\n'
                '\n'
                'source activate /gpfs/projects/bsc72/conda_envs/nbdsuite/0.0.1b3\n'
                '\n'
                'python -m nbdsuite.main input.yaml\n'
            )

    def _equibindDockingPoseRetriever(self, simulation_path):

        # Generating paths
        docked_jobs_origin = '3_docking_job/job/equibind_results' 
        docked_jobs_destination = '4_pele_simulation/docking_input/ligands'
        receptor_origin = '3_docking_job/job/equibind_calculations' 
        receptor_destination = '4_pele_simulation/docking_input/receptor'
        pele_simulation_path = simulation_path

        # Copying docked ligands and generating folders
        if not os.path.isdir(docked_jobs_destination):
            os.mkdir(docked_jobs_destination)

        for conformer in [x for x in os.listdir(docked_jobs_origin) if os.path.isdir(os.path.join(docked_jobs_origin,x))]:
            ligand_docked = os.listdir(os.path.join(docked_jobs_origin,conformer))[0]
            shutil.copy(os.path.join(docked_jobs_origin,conformer,ligand_docked), os.path.join(docked_jobs_destination,'{}.sdf'.format(conformer)))

            if not os.path.isdir(os.path.join(pele_simulation_path,conformer)):
                os.mkdir(os.path.join(pele_simulation_path,conformer))

        # Copying receptor
        receptor_folder = [x for x in os.listdir(receptor_origin) if os.path.isdir(os.path.join(receptor_origin,x))][0]
        receptor = [x for x in os.listdir(os.path.join(receptor_origin,receptor_folder)) if x.endswith('_protein.pdb')][0]

        if not os.path.isdir(receptor_destination):
            os.mkdir(receptor_destination)
            shutil.copy(os.path.join(receptor_origin,receptor_folder,receptor), os.path.join(receptor_destination,'{}.pdb'.format('receptor')))

    def _equibindPELEInputGenerator(self, previous_simulation_bool, simulation_path, force_field, truncated, perturbation_protocol, rescoring_method):

        docked_ligands_path = '4_pele_simulation/docking_input/ligands'
        receptor_path = '4_pele_simulation/docking_input/receptor'
        pele_simulation_path = simulation_path

        if not previous_simulation_bool:
        
            receptor = os.listdir(receptor_path)[0]

            print(' - Converting {} ligands to pdb.'.format(len(os.listdir(docked_ligands_path))))

            # Convert all the ligands to pdb and store them
            for ligand in os.listdir(docked_ligands_path):
                ligand_name = ligand.split('.')[0]
                self._PDBConversor(os.path.join(docked_ligands_path, ligand), os.path.join(pele_simulation_path,ligand_name))
                shutil.copy(os.path.join(receptor_path,receptor), os.path.join(pele_simulation_path,ligand_name))

            print(' - Merging the {} ligands to the receptor.'.format(len(os.listdir(docked_ligands_path))))

            # Merging all the inputs
            for ligand_folder in [x for x in os.listdir(pele_simulation_path) if x != '.ipynb_checkpoint']:
                ligand_folder_path = os.path.join(pele_simulation_path,ligand_folder)
                receptor = [x for x in os.listdir(ligand_folder_path) if x.startswith('receptor')][0]
                ligand = [x for x in os.listdir(ligand_folder_path) if x != receptor][0]
                self._PDBMerger(os.path.join(ligand_folder_path,receptor),os.path.join(ligand_folder_path,ligand))

        print(' - Generating files for the {} simulations.'.format(len(os.listdir(docked_ligands_path))))

        for ligand_folder in [x for x in os.listdir(pele_simulation_path) if x != '.ipynb_checkpoint']:
            working_path = os.path.join(pele_simulation_path, ligand_folder)
            input_simulation_file = os.listdir(working_path)[0]  
            self._PELESimulationFiles(working_path, input_simulation_file, force_field, truncated, perturbation_protocol, rescoring_method)   

    def setRdockToPELESimulation(self):
        # Rdock is going to be difficoult due to deprotonation.
        pass

    def setGlideToPELESimulation(self):
        # Glide also has the handicap that the separation of files has to be done
        # with a machine with a SHCRODINGER's license.
        pass
    
    def setEquibindToPELESimulation(self, rescoring_method, force_field=None, truncated=None, perturbation_protocol=None):

        forcefield_list, truncated_list, perturbation_list, rescoring_method_list = self._folderHierarchy(force_field, truncated, perturbation_protocol, rescoring_method)
        simulation_path = self._PELEJobManager(forcefield_list, truncated_list, perturbation_list, rescoring_method_list)
        previous_simulation_bool = self._PELEJobChecker(forcefield_list, truncated_list, perturbation_list, rescoring_method_list)
        self._equibindDockingPoseRetriever(simulation_path)
        self._equibindPELEInputGenerator(previous_simulation_bool, simulation_path, forcefield_list[1], truncated_list[1], perturbation_list[1], rescoring_method_list[1])
        
        
