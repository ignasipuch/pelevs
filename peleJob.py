import os 
import shutil
from openbabel import openbabel as ob
from Bio.PDB import PDBParser, PDBIO
import re
import pandas as pd

class PELE:
    """
    Attributes
    ==========
    sampling : str
        Sampling method chosen by the user.
    protein_portion : str 
        Selection of the protein with which PELE is going to
        perform the simulation.
    forcefield : str
        Force-field used to perform the PELE simulation.
    perturbation_protocol : str
        Type of perturbation used in the PELE simulation.
    docking_tool : str
        Docking tool used to obtain the ligand poses.

    Methods
    =======
    setGlideToPELESimulation(self, rescoring_method, force_field, truncated, perturbation_protocol)
        With the results of the Glide docking, prepares a PELE simulation to be run at MN4.
    setRdockToPELESimulation(self, rescoring_method, force_field, truncated, perturbation_protocol)
        With the results of the Equibind docking, prepares a PELE simulation to be run at MN4.
    setEquibindToPELESimulation(self, rescoring_method, force_field, truncated, perturbation_protocol)
        With the results of the rDock docking, prepares a PELE simulation to be run at MN4.

    Hidden Methods
    ==============
    _folderPreparation(self)
        Prepare necessary folders to store data.
    _folderHierarchy(self, force_field, truncated, perturbation_protocol, rescoring_method)
        Stores which variables have been set by the user.
    _PELEJobManager(self, forcefield_list, truncated_list, perturbation_protocol_list, rescoring_method_list)
        Creates hierarchical directories to store the inputs for the PELE simulations.
    _PELEJobChecker(self, forcefield_list, truncated_list, perturbation_protocol_list, rescoring_method_list)
        Checks if another pele job has been created previously to not repeat processes.
    _PDBConversor(self, file_in, path_out)
        Converts files from whatever format to pdb.
    _PDBMerger(self, receptor, ligand)
        Merges ligand and receptor into a single file.
    _PELESimulationFiles(self, path, pdb_file, force_field, truncated, perturbation_protocol, rescoring_method)
        Generates the yaml and the runner for each individual PELE simulation.
    _PELERunner(self, simulation_path)
        Generates a general runner.
    """

    def __init__(self):
        """
        Initialize object and assign atributes.
        """

        self.sampling = None
        self.protein_portion = None
        self.forcefield = None
        self.perturbation_protocol = None
        self.docking_tool = None

    def _folderPreparation(self):
        """
        Generate folders where data from the docking is going to be stored.
        """
        
        if not os.path.isdir('4_pele_simulation'):
            os.mkdir('4_pele_simulation')

        if not os.path.isdir('4_pele_simulation/docking_input'):
            os.mkdir('4_pele_simulation/docking_input')

        if not os.path.isdir('4_pele_simulation/pele_simulation'):
            os.mkdir('4_pele_simulation/pele_simulation')

    def _folderHierarchy(self, force_field, truncated, perturbation_protocol, rescoring_method):
        """
        Store which variables have been set by the user and to what have been set.

        Parameters
        ==========
        force_field : str 
            Forcefield of the simulation: opls or openff.
        truncated : str 
            Type of simulation for the PELE simulation: truncated or full.
        perturbation_protocol : str 
            Type of perturbation for the PELE simulation: if (induced fit), 
            refinement or minimization.
        rescoring_method : str 
            Amount of sampling for the PELE simulation: xshort, short, long, xlong.

        Returns
        =======
        forcefield_list : list
            List with [bool, str] where the string is the forcefield chosen (or default)
            and the bool indicates whether the forcefield option was chosen or default.
        truncated_list : list
            List with [bool, str] where the string is the protein portion chosen (or default)
            and the bool indicates whether the protein portion option was chosen or default.
        perturbation_list : list
            List with [bool, str] where the string is the perturbation chosen (or default)
            and the bool indicates whether the perturbation option was chosen or default.
        rescoring_method_list : list
            List with [bool, str] where the string is the rescoring method chosen (or default)
            and the bool indicates whether the rescoring method  option was chosen or default.
        """

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

        print(' - Sampling chosen: {}.'.format(rescoring_method))
        rescoring_method_list = [True, rescoring_method]

        self.sampling = rescoring_method
        self.protein_portion = truncated
        self.forcefield = force_field
        self.perturbation_protocol = perturbation_protocol

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
            if directory_path != path:
                shutil.rmtree(directory_path)

        print(' - Jobs will be stored at: {}'.format(os.path.join(pele_simulation_path, directory)))

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
            print(' - Previous simulation job found at: {}'.format(path_to_retrieve))
            print(' - Copying pdb files from: {}'.format(path_to_retrieve))

            for ligand in [x for x in os.listdir(path_to_retrieve) if x != 'general_runner.sh']:
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
 
    def _PELERunner(self, simulation_path):

        with open(os.path.join(simulation_path, 'general_runner.sh'), 'w') as fileout:
            fileout.writelines(
                'for d in *; if [ "$d" != "general_runner.sh" ]; then ; cd $d ; sbatch run_plat ; cd .. ; done\n'
            )

    def setGlideToPELESimulation(self, rescoring_method, force_field=None, truncated=None, perturbation_protocol=None):

        def _glideMaegztoPDB():

            maegz_to_pdb_path = '3_docking_job/job/output_pdb_files'

            if os.path.isdir(maegz_to_pdb_path):
                print(' - Splitting of the maegz file already done.')

            else:
                print(' - Splitting the outputted maegz file into individual pdbs.')
                print(' - Only the best tautomer/stereoisomer in saved.')

                os.system('$SCHRODINGER/run python3 dockprotocol/scripts/glide_to_pdb.py -jn {}'.format('glide_job'))

        def _glideDockingPoseRetriever(simulation_path):

            # Generating paths
            docked_jobs_origin = '3_docking_job/job/output_pdb_files' 
            docked_jobs_destination = '4_pele_simulation/docking_input/ligands'
            receptor_origin = '1_input_files/receptor/' 
            receptor_destination = '4_pele_simulation/docking_input/receptor'
            pele_simulation_path = simulation_path

            # Copying docked ligands and generating folders
            if not os.path.isdir(docked_jobs_destination):
                os.mkdir(docked_jobs_destination)

            for ligand_docked in [x for x in os.listdir(docked_jobs_origin) if os.path.isfile(os.path.join(docked_jobs_origin,x))]:
                shutil.copy(os.path.join(docked_jobs_origin, ligand_docked), os.path.join(docked_jobs_destination))

                if not os.path.isdir(os.path.join(pele_simulation_path,ligand_docked.split('.')[0])):
                    os.mkdir(os.path.join(pele_simulation_path,ligand_docked.split('.')[0]))

            # Copying receptor
            receptor = os.path.join(receptor_origin,os.listdir(receptor_origin)[0])

            if not os.path.isdir(receptor_destination):
                os.mkdir(receptor_destination)
                shutil.copy(receptor, os.path.join(receptor_destination,'{}.pdb'.format('receptor')))

        def _glidePELEInputGenerator(previous_simulation_bool, simulation_path, force_field, truncated, perturbation_protocol, rescoring_method):

            docked_ligands_path = '4_pele_simulation/docking_input/ligands'
            receptor_path = '4_pele_simulation/docking_input/receptor'
            pele_simulation_path = simulation_path

            if not previous_simulation_bool:
            
                receptor = os.listdir(receptor_path)[0]

                # Store all the ligands and receptor
                for ligand in os.listdir(docked_ligands_path):
                    ligand_name = ligand.split('.')[0]
                    shutil.copy(os.path.join(receptor_path,receptor), os.path.join(pele_simulation_path,ligand_name))
                    shutil.copy(os.path.join(docked_ligands_path,ligand), os.path.join(pele_simulation_path,ligand_name))

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

            print(' - Job created to run at MN4.')
            print(' - Send pele_simulation_folder to perform the simulations and run:\n bash glide_runner.sh.')

            self.docking_tool = 'glide'
  
        self._folderPreparation()
        forcefield_list, truncated_list, perturbation_list, rescoring_method_list = self._folderHierarchy(force_field, truncated, perturbation_protocol, rescoring_method)
        simulation_path = self._PELEJobManager(forcefield_list, truncated_list, perturbation_list, rescoring_method_list)
        previous_simulation_bool = self._PELEJobChecker(forcefield_list, truncated_list, perturbation_list, rescoring_method_list)
        _glideMaegztoPDB()
        _glideDockingPoseRetriever(simulation_path)
        _glidePELEInputGenerator(previous_simulation_bool, simulation_path, forcefield_list[1], truncated_list[1], perturbation_list[1], rescoring_method_list[1])
        self._PELERunner(simulation_path)

    def setRdockToPELESimulation(self, rescoring_method, force_field=None, truncated=None, perturbation_protocol=None):
        
        def _sdfSplitterAndSelector():

            path_docking = '3_docking_job/job/results'
            path_docked_ligands = '4_pele_simulation/docking_input/ligands'

            # Generating storage folders
            if not os.path.isdir(path_docked_ligands):
                os.mkdir(path_docked_ligands)
            
            print(' - Splitting all the outputted sdf files.')

            # sdf splitting
            for file in [x for x in os.listdir(path_docking) if x != '.ipynb_checkpoints']:
                input_file = os.path.join(path_docking, file)
                    
                with open(input_file, 'r') as file:
                    sdf_content = file.read()
                    
                previous_entry_id = None
                entries = sdf_content.strip().split("$$$$")
                auxiliar_counter = 0
                
                for entry in entries[:-1]:
                    auxiliar_counter += 1
                    count = auxiliar_counter - 50*((auxiliar_counter - 1)//50)
                    entry = entry.strip()

                    if entry:
                        entry_lines = entry.split("\n")
                        entry_id = None
                        for i in range(len(entry_lines)):
                            if entry_lines[i].strip().startswith(">  <s_lp_Variant>"):
                                entry_id = entry_lines[i+1].strip()
                                break

                        if entry_id == previous_entry_id:
                            previous_entry_id = entry_id
                            output_file = os.path.join(path_docked_ligands,'{entry_id}_{c}.sdf'.format(entry_id=entry_id, c=count))

                            with open(output_file, 'w') as outfile:
                                outfile.write(entry)
                        else:
                            previous_entry_id = entry_id
                            output_file = os.path.join(path_docked_ligands,'{entry_id}_{c}.sdf'.format(entry_id=entry_id, c=count))                 
                            with open(output_file, 'w') as outfile:
                                outfile.write(entry)

        def _rdockDockingPoseRetriever(simulation_path):

            # Generating paths
            docked_jobs_origin = '4_pele_simulation/docking_input/ligands'
            docking_job_path = '3_docking_job/job/'
            csv_path = '3_docking_job/rDock_best_poses.csv'

            receptor = [x for x in os.listdir(docking_job_path) if x.endswith('.mol2')][0]
            receptor_origin_path = os.path.join(docking_job_path, receptor)
            receptor_destination = '4_pele_simulation/docking_input/receptor'
            pele_simulation_path = simulation_path

            # Selecting which ligands we copy
            df = pd.read_csv(csv_path)
            df['file_name'] = df.apply(lambda row: f"{row['ligand']}-{row['conformer']}_{row['docking_conformation']}.sdf", axis=1)

            file_list = df['file_name'].tolist()

            # Copying docked ligands and generating folders
            if not os.path.isdir(pele_simulation_path):
                os.mkdir(pele_simulation_path)

            for conformer in file_list:
                ligand = conformer.split('-')[0]
                if not os.path.isdir(os.path.join(pele_simulation_path,ligand)):
                    os.mkdir(os.path.join(pele_simulation_path,ligand))

                shutil.copy(os.path.join(docked_jobs_origin,conformer), os.path.join(pele_simulation_path,ligand,'{}.sdf'.format(ligand)))

            # Copying receptor
            if not os.path.isdir(receptor_destination):
                os.mkdir(receptor_destination)
                shutil.copy(receptor_origin_path, os.path.join(receptor_destination,'{}.mol2'.format('receptor')))

        self._folderPreparation()
        forcefield_list, truncated_list, perturbation_list, rescoring_method_list = self._folderHierarchy(force_field, truncated, perturbation_protocol, rescoring_method)
        simulation_path = self._PELEJobManager(forcefield_list, truncated_list, perturbation_list, rescoring_method_list)
        previous_simulation_bool = self._PELEJobChecker(forcefield_list, truncated_list, perturbation_list, rescoring_method_list)
        _sdfSplitterAndSelector()
        _rdockDockingPoseRetriever(simulation_path)

 
    def setEquibindToPELESimulation(self, rescoring_method, force_field=None, truncated=None, perturbation_protocol=None):

        def _equibindDockingPoseRetriever(simulation_path):

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

        def _equibindPELEInputGenerator(previous_simulation_bool, simulation_path, force_field, truncated, perturbation_protocol, rescoring_method):

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

            self.docking_tool = 'equibind'

        self._folderPreparation()
        forcefield_list, truncated_list, perturbation_list, rescoring_method_list = self._folderHierarchy(force_field, truncated, perturbation_protocol, rescoring_method)
        simulation_path = self._PELEJobManager(forcefield_list, truncated_list, perturbation_list, rescoring_method_list)
        previous_simulation_bool = self._PELEJobChecker(forcefield_list, truncated_list, perturbation_list, rescoring_method_list)
        _equibindDockingPoseRetriever(simulation_path)
        _equibindPELEInputGenerator(previous_simulation_bool, simulation_path, forcefield_list[1], truncated_list[1], perturbation_list[1], rescoring_method_list[1])
        self._PELERunner(simulation_path)
        
