import os
import shutil
from openbabel import openbabel as ob
import re
import pandas as pd
import importlib.resources
import pelevs


class PELEJob:
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
    PELEDownloader(self)
        With all the simulations generated, creates a bash file to run all of them together.

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
        """
        Organize and create the hierarchy of folders required for the necessities of the 
        simulations wanted. If specified, the hierarchy will be:
        {force field}/{receptor portion}/{perturbation protocol}/{rescoring sampling}.

        Parameters
        ==========
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

        Returns
        =======
        simulation_path : str
            Path where all the PELE simulation files are going to be stored.
        """

        path = '4_pele_simulation/'
        pele_simulation_path = '4_pele_simulation/pele_simulation'

        lists = [forcefield_list, truncated_list,
                 perturbation_protocol_list, rescoring_method_list]
        directory = ''
        cont = 1

        for lst in lists:
            if lst[0] is True:
                if cont == 1:
                    directory_path = os.path.join(path, lst[1])
                    pele_directory_path = os.path.join(
                        pele_simulation_path, lst[1])

                directory_name = lst[1]
                new_directory = os.path.join(path, directory_name)

                if not os.path.isdir(new_directory):
                    os.makedirs(new_directory, exist_ok=True)

                directory = os.path.join(directory, directory_name)
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

        simulation_path = os.path.join(pele_simulation_path, directory)

        return simulation_path

    def _PELEJobChecker(self, forcefield_list, truncated_list, perturbation_protocol_list, rescoring_method_list):
        """
        Check if a previous PELE simulation job has been created to avoid loss of time. 
        If it exists, it copies the pdb file to new wanted PELE job directory.

        Parameters
        ==========
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

        Returns
        =======
        previous_simulation_bool : bool
            Boolean to determine whether or not a previous PELE job exists.
        """

        path = '4_pele_simulation/pele_simulation'

        directories = []
        directories.append(path)

        lists = [forcefield_list, truncated_list,
                 perturbation_protocol_list, rescoring_method_list]
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
                path = os.path.join(path, lst[1])
                directories.append(path)

        # Checks if there are other created jobs obtains path_to_retrieve
        for directory in directories[:-1]:
            content_in_directory = [x for x in os.listdir(
                directory) if not x.startswith('.')]
            if len(content_in_directory) > 1:
                previous_simulation_bool = True
                if len(content_in_directory) < 5:
                    previously_generated_simulation = [
                        x for x in content_in_directory if x not in options][0]
                    path_to_retrieve = os.path.join(
                        directory, previously_generated_simulation)
                else:
                    path_to_retrieve = directory
            else:
                path_to_retrieve = None

        # Copies files if there are files to copy
        if path_to_retrieve is not None:
            print(' - Previous simulation job found at: {}'.format(path_to_retrieve))
            print(' - Copying pdb files from: {}'.format(path_to_retrieve))

            for ligand in [x for x in os.listdir(path_to_retrieve) if x != 'general_runner.sh']:
                file = os.path.join(path_to_retrieve, ligand, ligand + '.pdb')
                path_to_ligand = os.path.join(path_simulation, ligand)
                if not os.path.isdir(path_to_ligand):
                    os.mkdir(path_to_ligand)

                shutil.copy(file, path_to_ligand)

        return previous_simulation_bool

    def _PDBConversor(self, file_in, path_out):
        """
        Converts the input file to PDB format using Open Babel.

        Parameters
        ==========
        file_in : str
            Name of the file to be converted.
        path_out : str
            Path to store the output file.
        """

        file_name, file_format = file_in.split('.')
        file_name = os.path.basename(file_name)

        if not file_format == 'pdb':

            conv = ob.OBConversion()
            conv.SetInAndOutFormats(file_in.split('.')[-1], 'pdb')
            mol = ob.OBMol()
            conv.ReadFile(mol, file_in)
            conv.WriteFile(mol, os.path.join(
                path_out, '{}.pdb'.format(file_name)))

    def _PDBMerger(self, receptor, ligand):
        """
        Merges receptor and ligand into a single file with the characteristics required by
        PELE to work properly.

        Parameters
        ==========
        receptor : str
            Path to the receptor file we want to merge.
        ligand : str
            Path to the ligand file we want to merge.
        """

        def _receptorModifier(receptor):
            """
            Changes residues numbers and counts number of atoms in the receptor.

            Parameters
            ==========
            receptor : str
                Path to the receptor file we want to merge.

            Returns
            =======
            file_mod_prot : str
                Path to the modified version of the receptor.
            ligand_cont_num : int
                Number of atoms in the receptor.
            """

            file_mod_prot = receptor.split('.pdb')[0] + '_mod.pdb'

            with open(receptor, 'r') as pdb_file:
                lines = pdb_file.readlines()

            modified_lines = []
            residue_index = 1000
            previous_residue = None
            previous_chain_residue = None

            # Changeing residue number
            for line in lines:

                # Change lines of the receptor
                if line.startswith('ATOM'):
                    residue_letters = line[17:20].strip()
                    chain_residue = line[21:26].strip()

                    if (residue_letters != previous_residue) or (chain_residue != previous_chain_residue):
                        residue_index += 1

                    line = line[:22] + str(residue_index) + line[26:]
                    modified_lines.append(line)
                    previous_residue = residue_letters
                    previous_chain_residue = chain_residue

                elif line.startswith('TER'):
                    pass

                # Change lines of waters/metals
                elif line.startswith('HETATM'):
                    residue_letters = line[17:20].strip()
                    chain_residue = line[21:26].strip()

                    if (residue_letters != previous_residue) or (chain_residue != previous_chain_residue):
                        residue_index += 1

                    line = line[:22] + str(residue_index) + line[26:]
                    modified_lines.append(line)
                    previous_residue = residue_letters
                    previous_chain_residue = chain_residue

                else:
                    pass

            # Write the modified lines back to the PDB file
            with open(file_mod_prot, 'w') as pdb_file:
                pdb_file.writelines(modified_lines)

            # Counting number of atoms of the protein
            ligand_cont_num = 0

            with open(file_mod_prot) as filein:
                for line in filein:
                    sline = line.split()
                    if sline[0] == 'ATOM':
                        ligand_cont_num += 1

            return file_mod_prot, ligand_cont_num

        def _ligandAtomChainNumberModifier(ligand):
            """
            Modifies the ligand file to have the characteristics needed for the PELE simulation:
            like the ligand chain or the chain name.

            Parameters
            ==========
            ligand : str
                Path to the ligand file we want to merge.

            Returns
            =======
            file_mod_lig : str
                Path to the modified version of the ligand.
            """

            file_mod1_lig = ligand.split('.pdb')[0] + '_m.pdb'
            file_mod_lig = ligand.split('.pdb')[0] + '_mod.pdb'

            # Modifying the atom names
            cont = 1
            with open(ligand) as filein:
                with open(file_mod1_lig, 'w') as fileout:
                    for line in filein:
                        if line.startswith('HETATM'):
                            sline = line.split()
                            beginning_of_line = line[0:13]
                            end_of_line = line[17:]
                            atom = ''.join(
                                [c for c in sline[2] if c.isalpha()])
                            new_line = beginning_of_line + \
                                (atom + str(cont)).ljust(4) + end_of_line
                            fileout.writelines(new_line)
                            cont += 1

            # Change the chain name and number
            with open(file_mod1_lig, 'r') as filein:
                with open(file_mod_lig, 'w') as fileout:
                    for line in filein:
                        if re.search('UN.......', line):
                            line = re.sub('UN.......', 'LIG L 900', line)
                            fileout.writelines(line)
                        else:
                            fileout.writelines(line)

            return file_mod_lig

        def _receptorLigandMerger(receptor, ligand):
            """
            Merges the receptor and the ligand into a single file.

            Parameters
            ==========
            receptor : str
                Path to the receptor file we want to merge.
            ligand : str
                Path to the ligand file we want to merge.

            Returns
            =======
            writing_path : str
                Path to the newly created pdb with both molecules.
            """

            output_file = 'intermediate.pdb'
            path_files = os.path.dirname(ligand)
            writing_path = os.path.join(path_files, output_file)

            # Joining ligand and receptor and directing warnings to out.txt
            os.system('obabel {receptor} {ligand} -O {output} 2> out.txt'.format(
                receptor=receptor, ligand=ligand, output=writing_path))

            return writing_path

        def _inputAdapter(merged, output_file, ligand_cont_num):
            """
            Changes the atom numeration of the ligand according to the number
            of atoms of the receptor.

            Parameters
            ==========
            merged : str
                Path to the merged file.
            output_file : str
                Path of the definitive merged file.
            ligand_cont_num : int
                Number of atoms of the receptor.
            """

            # Changing the atom numeration of the ligand.
            with open(merged, 'r') as filein:
                with open(output_file, 'w') as fileout:
                    for line in filein:
                        sline = line.split()
                        if sline[0] == 'ATOM' or sline[0] == 'HETATM':
                            if re.search('HETATM.....', line):
                                # Considering length of the protein
                                if len(str(ligand_cont_num + 1)) < 5:
                                    line = re.sub(
                                        'HETATM.....', 'HETATM ' + str(ligand_cont_num + 1), line)
                                elif len(str(ligand_cont_num + 1)) == 5:
                                    line = re.sub(
                                        'HETATM.....', 'HETATM' + str(ligand_cont_num + 1), line)

                                fileout.writelines(line)
                                ligand_cont_num += 1
                            else:
                                fileout.writelines(line)

        def _intermediateFilesRemover(output_file):
            """
            Removes all the intermediate files that have been 
            generated in the process.

            Parameters
            ==========
            output_file : str
                Path of the definitive merged file.
            """

            working_directory = os.path.dirname(output_file)
            input_PELE_file = os.path.basename(output_file)

            for file in os.listdir(working_directory):
                if file != input_PELE_file:
                    os.remove(os.path.join(working_directory, file))

        ligand_name = os.path.basename(ligand)
        working_directory = os.path.dirname(ligand)
        output_file = os.path.join(working_directory, ligand_name)

        file_mod_prot, ligand_cont_num = _receptorModifier(receptor)
        file_mod_lig = _ligandAtomChainNumberModifier(ligand)
        intermediate = _receptorLigandMerger(file_mod_prot, file_mod_lig)
        _inputAdapter(intermediate, output_file, ligand_cont_num)
        _intermediateFilesRemover(output_file)

    def _PELESimulationFiles(self, path, pdb_file, force_field, truncated, perturbation_protocol, rescoring_method):
        """
        Generates the yaml and the run files for the specific conditions 
        inputted by the user.

        Parameters
        ==========
        path : str
            Path where the simulation is stored.
        pdb_file : str
            Name of the pdb file for which the yaml and run are needed.
        force_field : str
            Name of the forcefield chosen (or default)
        truncated : str
            Name of the protein portion chosen (or default)
        perturbation_protocol : str
            Name of the perturbation chosen (or default)
        rescoring_method : str
            Name of the rescoring method chosen (or default)
        """

        # YAML file generation
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
                raise Exception(
                    'InvalidProtocol: You have to choose either if or refinement.')

            if force_field == 'opls':
                fileout.writelines('pele_forcefield: opls2005\n')
            elif force_field == 'openff':
                fileout.writelines('pele_forcefield: openff-2.0.0\n')
            else:
                raise Exception(
                    'InvalidForceField: You have to choose either opls or openff.')

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

        # Run file generation
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
                'source activate /gpfs/projects/bsc72/conda_envs/nbdsuite/0.0.1b4\n'
                '\n'
                'python -m nbdsuite.main input.yaml\n'
            )

    def _PELERunner(self, simulation_path):
        """
        Generates a general runner that iterates over the individual runners.

        Parameters
        ==========
        simulation_path : str
            Path were the simulation is being stored.

        """
        with open(os.path.join(simulation_path, 'general_runner.sh'), 'w') as fileout:
            fileout.writelines(
                'for d in *; do if [ "$d" != "general_runner.sh" ]; then cd "$d"; sbatch run_plat; cd ..; fi; done\n'
            )

    def setGlideToPELESimulation(self, rescoring_method, force_field=None, truncated=None, perturbation_protocol=None):
        """
        From a Glide docking, obtain the best scores per ligand and the conformations associated
        to then prepare a PELE simulation with certain variables to be chosen.

        Parameters
        ==========
        rescoring_method : str
            Length of the sampling in the PELE simulation: 
                1. xshort (cpus = 8; epochs = 1; steps = 5)
                2. short (cpus = 16; epochs = 1; steps = 10)
                3. long (cpus = 32; epochs = 5; steps = 10)
                4. xlong (cpus = 48; epochs = 5; steps = 25)
        force_field : str
            Force field to be used in the PELE simulation: 
                1. opls (opls2005) -> Default
                2. openff (openff-2.0.0)
        truncated : str
            Portion of the receptor to be used in the simulation:
                1. truncated (flexible_region_radius = 8; frozen_region_radius = 13) -> Default
                2. full (whole protein)
        perturbation_protocol : str
            Strength of the perturbations tried by PELE in each Monte Carlo step:
                1. minimization (pele_perturbation_level = 0)
                2. refinement (pele_perturbation_level = 2) -> Default
                3. if (induced_fit: pele_perturbation_level = 3)
        """

        def _glideMaegztoPDB():
            """
            If not done previously, split the Glide's output maegz file into its components 
            in pdb format by using the script in: scripts/glide_to_pdb.py.
            """

            maegz_to_pdb_path = '3_docking_job/job/output_pdb_files'

            if os.path.isdir(maegz_to_pdb_path):
                print(' - Splitting of the maegz file already done.')

            else:
                print(' - Splitting the outputted maegz file into individual pdbs.')
                print(' - Only the best tautomer/stereoisomer is saved.')

                # Generating folder
                os.mkdir(maegz_to_pdb_path)

                os.system(
                    '$SCHRODINGER/run python3 scripts/glide_to_pdb.py -jn {}'.format('glide_job'))

        def _glideDockingPoseRetriever(simulation_path):
            """
            Creation of necessary folders and copying important files.

            Parameters
            ==========
            simulation_path : str
                Path where the PELE simulation is being stored.
            """

            # Generating paths
            docked_jobs_origin = '3_docking_job/job/output_pdb_files'
            docked_jobs_destination = '4_pele_simulation/docking_input/ligands'
            receptor_origin = '1_input_files/receptor/'
            receptor_destination = '4_pele_simulation/docking_input/receptor'
            pele_simulation_path = simulation_path

            # Copying docked ligands and generating folders
            if not os.path.isdir(docked_jobs_destination):
                os.mkdir(docked_jobs_destination)

            for ligand_docked in [x for x in os.listdir(docked_jobs_origin) if os.path.isfile(os.path.join(docked_jobs_origin, x))]:
                shutil.copy(os.path.join(docked_jobs_origin, ligand_docked), os.path.join(
                    docked_jobs_destination))

                if not os.path.isdir(os.path.join(pele_simulation_path, ligand_docked.split('.')[0])):
                    os.mkdir(os.path.join(pele_simulation_path,
                                          ligand_docked.split('.')[0]))

            # Copying receptor
            receptor = os.path.join(
                receptor_origin, os.listdir(receptor_origin)[0])

            if not os.path.isdir(receptor_destination):
                os.mkdir(receptor_destination)
                shutil.copy(receptor, os.path.join(
                    receptor_destination, '{}.pdb'.format('receptor')))

        def _glidePELEInputGenerator(previous_simulation_bool, simulation_path, force_field, truncated, perturbation_protocol, rescoring_method):
            """
            Merging receptors and ligands and generating the necessary yamls and runs.

            Parameters
            ==========
            previous_simulation_bool : bool
                Boolean to determine whether or not a previous PELE job exists.
            simulation_path : str
                Path where the PELE simulation is being stored.
            rescoring_method : str
                Length of the sampling in the PELE simulation.
            force_field : str
                Force field to be used in the PELE simulation.
            truncated : str
                Portion of the receptor to be used in the simulation.
            perturbation_protocol : str
                Strength of the perturbations tried by PELE in each Monte Carlo step.
            """

            docked_ligands_path = '4_pele_simulation/docking_input/ligands'
            receptor_path = '4_pele_simulation/docking_input/receptor'
            pele_simulation_path = simulation_path

            if not previous_simulation_bool:

                receptor = os.listdir(receptor_path)[0]

                # Store all the ligands and receptor
                for ligand in os.listdir(docked_ligands_path):
                    ligand_name = ligand.split('.')[0]
                    shutil.copy(os.path.join(receptor_path, receptor),
                                os.path.join(pele_simulation_path, ligand_name))
                    shutil.copy(os.path.join(docked_ligands_path, ligand), os.path.join(
                        pele_simulation_path, ligand_name))

                print(
                    ' - Merging the {} ligands to the receptor...'.format(len(os.listdir(docked_ligands_path))))

                # Merging all the inputs
                for ligand_folder in [x for x in os.listdir(pele_simulation_path) if (x != '.ipynb_checkpoint') and (x != 'general_runner.sh')]:
                    ligand_folder_path = os.path.join(
                        pele_simulation_path, ligand_folder)
                    receptor = [x for x in os.listdir(
                        ligand_folder_path) if x.startswith('receptor')][0]
                    ligand = [x for x in os.listdir(
                        ligand_folder_path) if x != receptor][0]
                    self._PDBMerger(os.path.join(ligand_folder_path, receptor), os.path.join(
                        ligand_folder_path, ligand))

            # Deleting openbabel merging information.
            if os.path.isfile('out.txt'):
                os.remove('out.txt')

            print(
                ' - Generating yaml and run files.')

            # Generating yamls and runs
            for ligand_folder in [x for x in os.listdir(pele_simulation_path) if (x != '.ipynb_checkpoint') and (x != 'general_runner.sh')]:
                working_path = os.path.join(
                    pele_simulation_path, ligand_folder)
                input_simulation_file = os.listdir(working_path)[0]
                self._PELESimulationFiles(working_path, input_simulation_file,
                                          force_field, truncated, perturbation_protocol, rescoring_method)

            print(' - Job created to run at MN4.')
            print(
                ' - Send pele_simulation_folder to perform the simulations and run:\n   bash general_runner.sh.')

            self.docking_tool = 'glide'

        self._folderPreparation()
        forcefield_list, truncated_list, perturbation_list, rescoring_method_list = self._folderHierarchy(
            force_field, truncated, perturbation_protocol, rescoring_method)
        simulation_path = self._PELEJobManager(
            forcefield_list, truncated_list, perturbation_list, rescoring_method_list)
        previous_simulation_bool = self._PELEJobChecker(
            forcefield_list, truncated_list, perturbation_list, rescoring_method_list)
        _glideMaegztoPDB()
        _glideDockingPoseRetriever(simulation_path)
        _glidePELEInputGenerator(previous_simulation_bool, simulation_path,
                                 forcefield_list[1], truncated_list[1], perturbation_list[1], rescoring_method_list[1])
        self._PELERunner(simulation_path)

    def setRdockToPELESimulation(self, rescoring_method, force_field=None, truncated=None, perturbation_protocol=None):
        """
        From an rDock docking, obtain the best scores per ligand and the conformations associated
        to then prepare a PELE simulation with certain variables to be chosen.

        Parameters
        ==========
        rescoring_method : str
            Length of the sampling in the PELE simulation: 
                1. xshort (cpus = 8; epochs = 1; steps = 5)
                2. short (cpus = 16; epochs = 1; steps = 10)
                3. long (cpus = 32; epochs = 5; steps = 10)
                4. xlong (cpus = 48; epochs = 5; steps = 25)
        force_field : str
            Force field to be used in the PELE simulation: 
                1. opls (opls2005) -> Default
                2. openff (openff-2.0.0)
        truncated : str
            Portion of the receptor to be used in the simulation:
                1. truncated (flexible_region_radius = 8; frozen_region_radius = 13) -> Default
                2. full (whole protein)
        perturbation_protocol : str
            Strength of the perturbations tried by PELE in each Monte Carlo step:
                1. minimization (pele_perturbation_level = 0)
                2. refinement (pele_perturbation_level = 2) -> Default
                3. if (induced_fit: pele_perturbation_level = 3)
        """

        def _sdfSplitterAndSelector():
            """
            If not done previously, split the rDock's output sdf files into its components.
            """

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
                            output_file = os.path.join(
                                path_docked_ligands, '{entry_id}_{c}.sdf'.format(entry_id=entry_id, c=count))

                            with open(output_file, 'w') as outfile:
                                outfile.write(entry)
                        else:
                            previous_entry_id = entry_id
                            output_file = os.path.join(
                                path_docked_ligands, '{entry_id}_{c}.sdf'.format(entry_id=entry_id, c=count))
                            with open(output_file, 'w') as outfile:
                                outfile.write(entry)

        def _rdockDockingPoseRetriever(simulation_path):
            """
            Creation of necessary folders, copying the receptor, selecting ligands with 
            best scores, copying their structures and deleting the rest.

            Parameters
            ==========
            simulation_path : str
                Path where the PELE simulation is being stored.
            """

            # Generating paths
            docked_jobs_origin = '4_pele_simulation/docking_input/ligands'
            docking_job_path = '3_docking_job/job/'
            csv_path = '3_docking_job/rDock_best_poses.csv'

            receptor = [x for x in os.listdir(
                docking_job_path) if x.endswith('.mol2')][0]
            receptor_origin_path = os.path.join(docking_job_path, receptor)
            receptor_destination = '4_pele_simulation/docking_input/receptor'
            pele_simulation_path = simulation_path

            # Selecting which ligands we copy
            df = pd.read_csv(csv_path)
            df['file_name'] = df.apply(
                lambda row: f"{row['ligand']}-{row['conformer']}_{row['docking_conformation']}.sdf", axis=1)

            file_list = df['file_name'].tolist()

            # Deleting ligands
            for file_name in os.listdir(docked_jobs_origin):
                file_path = os.path.join(docked_jobs_origin, file_name)
                ligand_name = file_name.split('-')[0] + '.sdf'
                new_ligand_path = os.path.join(docked_jobs_origin, ligand_name)

                if file_name not in file_list:
                    if os.path.isfile(file_path):
                        os.remove(file_path)
                else:
                    os.rename(file_path, new_ligand_path)

            # Copying docked ligands and generating folders
            if not os.path.isdir(pele_simulation_path):
                os.mkdir(pele_simulation_path)

            for conformer in os.listdir(docked_jobs_origin):
                ligand = conformer.split('.')[0]
                if not os.path.isdir(os.path.join(pele_simulation_path, ligand)):
                    os.mkdir(os.path.join(pele_simulation_path, ligand))

            # Copying receptor
            if not os.path.isdir(receptor_destination):
                os.mkdir(receptor_destination)
                shutil.copy(receptor_origin_path, os.path.join(
                    receptor_destination, '{}.mol2'.format('receptor')))

        def _rdockPELEInputGenerator(previous_simulation_bool, simulation_path, force_field, truncated, perturbation_protocol, rescoring_method):
            """
            Converting ligands and receptor to then merge them and generate the necessary yamls and runs.

            Parameters
            ==========
            previous_simulation_bool : bool
                Boolean to determine whether or not a previous PELE job exists.
            simulation_path : str
                Path where the PELE simulation is being stored.
            rescoring_method : str
                Length of the sampling in the PELE simulation.
            force_field : str
                Force field to be used in the PELE simulation.
            truncated : str
                Portion of the receptor to be used in the simulation.
            perturbation_protocol : str
                Strength of the perturbations tried by PELE in each Monte Carlo step.
            """

            docked_ligands_path = '4_pele_simulation/docking_input/ligands'
            receptor_path = '4_pele_simulation/docking_input/receptor'
            pele_simulation_path = simulation_path

            if not previous_simulation_bool:

                receptor = os.listdir(receptor_path)[0]
                converted_receptor = 'receptor.pdb'

                print(' - Converting receptor from mol2 to pdb.')

                # Converting receptor from mol2 to pdb.
                self._PDBConversor(os.path.join(
                    receptor_path, receptor), receptor_path)

                print(
                    ' - Converting {} ligands to pdb.'.format(len(os.listdir(docked_ligands_path))))

                # Convert all the ligands to pdb and store them
                for ligand in os.listdir(docked_ligands_path):
                    ligand_name = ligand.split('.')[0]
                    self._PDBConversor(os.path.join(docked_ligands_path, ligand), os.path.join(
                        pele_simulation_path, ligand_name))
                    shutil.copy(os.path.join(receptor_path, converted_receptor), os.path.join(
                        pele_simulation_path, ligand_name))

                print(
                    ' - Merging the {} ligands to the receptor.'.format(len(os.listdir(docked_ligands_path))))

                # Merging all the inputs
                for ligand_folder in [x for x in os.listdir(pele_simulation_path) if x != '.ipynb_checkpoint']:
                    ligand_folder_path = os.path.join(
                        pele_simulation_path, ligand_folder)
                    receptor = [x for x in os.listdir(
                        ligand_folder_path) if x.startswith('receptor')][0]
                    ligand = [x for x in os.listdir(
                        ligand_folder_path) if x != receptor][0]
                    self._PDBMerger(os.path.join(ligand_folder_path, receptor), os.path.join(
                        ligand_folder_path, ligand))

            print(
                ' - Generating yaml and run files.')

            # Generating yamls and runs
            for ligand_folder in [x for x in os.listdir(pele_simulation_path) if (x != '.ipynb_checkpoint') and (x != 'general_runner.sh')]:
                working_path = os.path.join(
                    pele_simulation_path, ligand_folder)
                input_simulation_file = [x for x in os.listdir(
                    working_path) if x.endswith('.pdb')][0]
                self._PELESimulationFiles(working_path, input_simulation_file,
                                          force_field, truncated, perturbation_protocol, rescoring_method)

            print(' - Job created to run at MN4.')
            print(
                ' - Send pele_simulation_folder to perform the simulations and run:\n   bash general_runner.sh.')

            self.docking_tool = 'rdock'

        self._folderPreparation()
        forcefield_list, truncated_list, perturbation_list, rescoring_method_list = self._folderHierarchy(
            force_field, truncated, perturbation_protocol, rescoring_method)
        simulation_path = self._PELEJobManager(
            forcefield_list, truncated_list, perturbation_list, rescoring_method_list)
        previous_simulation_bool = self._PELEJobChecker(
            forcefield_list, truncated_list, perturbation_list, rescoring_method_list)
        _sdfSplitterAndSelector()
        _rdockDockingPoseRetriever(simulation_path)
        _rdockPELEInputGenerator(previous_simulation_bool, simulation_path,
                                 forcefield_list[1], truncated_list[1], perturbation_list[1], rescoring_method_list[1])
        self._PELERunner(simulation_path)

    def setEquibindToPELESimulation(self, rescoring_method, force_field=None, truncated=None, perturbation_protocol=None):
        """
        From an Equibind docking prepare a PELE simulation with certain variables to be chosen. 
        Since Equibind does not yield any score, all the outputs from the docking are simulated with
        PELE. The nomenclature for each folder/pdb is: {ligand}_{ligprep_conformer}.

        Parameters
        ==========
        rescoring_method : str
            Length of the sampling in the PELE simulation: 
                1. xshort (cpus = 8; epochs = 1; steps = 5)
                2. short (cpus = 16; epochs = 1; steps = 10)
                3. long (cpus = 32; epochs = 5; steps = 10)
                4. xlong (cpus = 48; epochs = 5; steps = 25)
        force_field : str
            Force field to be used in the PELE simulation: 
                1. opls (opls2005) -> Default
                2. openff (openff-2.0.0)
        truncated : str
            Portion of the receptor to be used in the simulation:
                1. truncated (flexible_region_radius = 8; frozen_region_radius = 13) -> Default
                2. full (whole protein)
        perturbation_protocol : str
            Strength of the perturbations tried by PELE in each Monte Carlo step:
                1. minimization (pele_perturbation_level = 0)
                2. refinement (pele_perturbation_level = 2) -> Default
                3. if (induced_fit: pele_perturbation_level = 3)
        """

        def _equibindDockingPoseRetriever(simulation_path):
            """
            Creation of necessary folders, copying the receptor and all the 
            ligands.

            Parameters
            ==========
            simulation_path : str
                Path where the PELE simulation is being stored.
            """

            # Generating paths
            docked_jobs_origin = '3_docking_job/job/equibind_results'
            docked_jobs_destination = '4_pele_simulation/docking_input/ligands'
            receptor_origin = '3_docking_job'
            receptor_destination = '4_pele_simulation/docking_input/receptor'
            pele_simulation_path = simulation_path

            # Copying docked ligands and generating folders
            if not os.path.isdir(docked_jobs_destination):
                os.mkdir(docked_jobs_destination)

            for conformer in [x for x in os.listdir(docked_jobs_origin) if os.path.isdir(os.path.join(docked_jobs_origin, x))]:
                ligand_docked = os.listdir(os.path.join(
                    docked_jobs_origin, conformer))[0]
                shutil.copy(os.path.join(docked_jobs_origin, conformer, ligand_docked), os.path.join(
                    docked_jobs_destination, '{}.sdf'.format(conformer)))

                if not os.path.isdir(os.path.join(pele_simulation_path, conformer)):
                    os.mkdir(os.path.join(pele_simulation_path, conformer))

            # Copying receptor
            receptor = [x for x in os.listdir(
                receptor_origin) if x.endswith('.pdb')][0]
            print(
                ' - Reference receptor: {}'.format(os.path.join(receptor_origin, receptor)))

            if not os.path.isdir(receptor_destination):
                os.mkdir(receptor_destination)
                shutil.copy(os.path.join(receptor_origin, receptor), os.path.join(
                    receptor_destination, '{}.pdb'.format('receptor')))

        def _equibindPELEInputGenerator(previous_simulation_bool, simulation_path, force_field, truncated, perturbation_protocol, rescoring_method):
            """
            If not done previously, converting ligands to then merge them to the receptor
            and generate the necessary yamls and runs.

            Parameters
            ==========
            previous_simulation_bool : bool
                Boolean to determine whether or not a previous PELE job exists.
            simulation_path : str
                Path where the PELE simulation is being stored.
            rescoring_method : str
                Length of the sampling in the PELE simulation.
            force_field : str
                Force field to be used in the PELE simulation.
            truncated : str
                Portion of the receptor to be used in the simulation.
            perturbation_protocol : str
                Strength of the perturbations tried by PELE in each Monte Carlo step.
            """

            docked_ligands_path = '4_pele_simulation/docking_input/ligands'
            receptor_path = '4_pele_simulation/docking_input/receptor'
            pele_simulation_path = simulation_path

            if not previous_simulation_bool:

                receptor = os.listdir(receptor_path)[0]

                print(
                    ' - Converting {} ligands to pdb.'.format(len(os.listdir(docked_ligands_path))))

                # Convert all the ligands to pdb and store them
                for ligand in os.listdir(docked_ligands_path):
                    ligand_name = ligand.split('.')[0]
                    self._PDBConversor(os.path.join(docked_ligands_path, ligand), os.path.join(
                        pele_simulation_path, ligand_name))
                    shutil.copy(os.path.join(receptor_path, receptor),
                                os.path.join(pele_simulation_path, ligand_name))

                print(
                    ' - Merging the {} ligands to the receptor.'.format(len(os.listdir(docked_ligands_path))))

                # Merging all the inputs
                for ligand_folder in [x for x in os.listdir(pele_simulation_path) if x != '.ipynb_checkpoint']:
                    ligand_folder_path = os.path.join(
                        pele_simulation_path, ligand_folder)
                    receptor = [x for x in os.listdir(
                        ligand_folder_path) if x.startswith('receptor')][0]
                    ligand = [x for x in os.listdir(
                        ligand_folder_path) if x != receptor][0]
                    self._PDBMerger(os.path.join(ligand_folder_path, receptor), os.path.join(
                        ligand_folder_path, ligand))

            print(
                ' - Generating yaml and run files.')

            # Generating yamls and runs
            for ligand_folder in [x for x in os.listdir(pele_simulation_path) if x != '.ipynb_checkpoint']:
                working_path = os.path.join(
                    pele_simulation_path, ligand_folder)
                input_simulation_file = [x for x in os.listdir(
                    working_path) if x.endswith('.pdb')][0]
                self._PELESimulationFiles(working_path, input_simulation_file,
                                          force_field, truncated, perturbation_protocol, rescoring_method)

            print(' - Job created to run at MN4.')
            print(
                ' - Send pele_simulation_folder to perform the simulations and run:\n   bash general_runner.sh.')

            self.docking_tool = 'equibind'

        self._folderPreparation()
        forcefield_list, truncated_list, perturbation_list, rescoring_method_list = self._folderHierarchy(
            force_field, truncated, perturbation_protocol, rescoring_method)
        simulation_path = self._PELEJobManager(
            forcefield_list, truncated_list, perturbation_list, rescoring_method_list)
        previous_simulation_bool = self._PELEJobChecker(
            forcefield_list, truncated_list, perturbation_list, rescoring_method_list)
        _equibindDockingPoseRetriever(simulation_path)
        _equibindPELEInputGenerator(previous_simulation_bool, simulation_path,
                                    forcefield_list[1], truncated_list[1], perturbation_list[1], rescoring_method_list[1])
        self._PELERunner(simulation_path)

    def PELEDownloader(self):
        """
        Generates and copies files in order to make the volume of data downloaded
        as small as possible. This method is thought to be used after all the jobs wanted 
        have been generated.
        """

        # Generating paths
        path_destination = self.docking_tool
        path_script = 'dockprotocol/scripts/pele_downloader.py'
        path_simulation = '4_pele_simulation'
        file_name = 'download_files.sh'
        download_file_path = os.path.join(path_simulation, file_name)

        # Copying script
        shutil.copy(path_script, path_simulation)

        # Generating runner
        with open(download_file_path, 'w') as filein:
            filein.writelines(
                'python pele_downloader.py -o {}'.format(path_destination)
            )

        print(' - Run:\n   bash {}\n   After PELE simulations have been performed.'.format(file_name))
