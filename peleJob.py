import os 
import shutil
from openbabel import openbabel as ob
from Bio.PDB import PDBParser, PDBIO, Select

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

    def _PDBConversor(self, file_in, path_out):
        """
        Converts the input file to PDB format using Open Babel.
        """

        file_name, file_format = file_in.split('.')

        if not file_format == 'pdb':

            conv = ob.OBConversion()
            conv.SetInAndOutFormats(file_in.split('.')[-1], 'pdb')
            mol = ob.OBMol()
            conv.ReadFile(mol, file_in)
            conv.WriteFile(mol, os.path.join(path_out,'{}.pdb'.format(file_name)))

    def _PDBMerger(self, receptor, ligand):
        """
        Merges the receptor and ligand into a single PDB file with a single entry.
        """

        file_name = ligand.split('.')[0]

        merged_pdb = ligand
        storage_path = os.path.join('4_pele_simulation/pele_simulation',file_name,merged_pdb)

        parser = PDBParser()

        # Parse the receptor and ligand PDB files
        receptor_structure = parser.get_structure('receptor', receptor)
        ligand_structure = parser.get_structure('ligand', ligand)

        # Create a new structure to hold the merged receptor and ligand
        merged_structure = parser.get_structure('merged', receptor)

        # Get the model and chain objects of the merged structure
        merged_model = merged_structure[0]
        merged_chain_receptor = merged_model['A']
        merged_chain_ligand = merged_model['L']

        # Get the last residue number of the receptor
        last_residue_number = max(
            residue.id[1] for chain in receptor_structure.get_chains() for residue in chain
        )

        # Add the receptor atoms to the merged structure
        receptor_model = receptor_structure[0]
        for chain in receptor_model:
            for residue in chain:
                residue.id = (' ', residue.id[1] + 1, ' ')
                for atom in residue:
                    atom_copy = atom.copy()
                    atom_copy.parent = merged_chain_receptor
                    merged_chain_receptor.add(atom_copy)

        # Add the ligand atoms to the merged structure
        ligand_model = ligand_structure[0]
        for chain in ligand_model:
            for residue in chain:
                residue.id = (' ', last_residue_number + 1, ' ')
                residue_copy = residue.copy()
                for atom in residue:
                    atom_copy = atom.copy()
                    atom_copy.parent = residue_copy
                    residue_copy.add(atom_copy)
                merged_chain_ligand.add(residue_copy)
                last_residue_number += 1

        # Save the merged structure to a PDB file
        io = PDBIO()
        io.set_structure(storage_path)
        io.save(storage_path, Select())

    def _equibindDockingPoseRetriever(self):

        # Generating paths
        docked_jobs_origin = '3_docking_job/job/equibind_results' 
        docked_jobs_destination = '4_pele_simulation/docking_input/ligands'
        receptor_origin = '3_docking_job/job/equibind_calculations' 
        receptor_destination = '4_pele_simulation/docking_input/receptor'
        pele_simulation_path = '4_pele_simulation/pele_simulation'

        # Copying docked ligands and generating folders
        if not os.path.isdir(docked_jobs_destination):
            os.mkdir(docked_jobs_destination)

        for conformer in [x for x in os.listdir(docked_jobs_origin) if os.path.isdir(x)]:
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

    def _equibindPELEInputGenerator(self):

        docked_ligands_path = '4_pele_simulation/docking_input/ligands'
        receptor_path = '4_pele_simulation/docking_input/receptor'
        pele_simulation_path = '4_pele_simulation/pele_simulation'

        receptor = os.listdir(receptor_path)[0]

        # Convert all the ligands to pdb and store them
        for ligand in os.listdir(docked_ligands_path):
            ligand_name = ligand.split('.')[0]
            self._PDBConversor(os.path.join(docked_ligands_path, ligand), os.path.join(pele_simulation_path,ligand_name))
            shutil.copy(os.path.join(receptor_path,receptor), os.path.join(pele_simulation_path,ligand_name))

        # Merging all the inputs
        for ligand_folder in os.listdir(pele_simulation_path):
            ligand_folder_path = os.path.join(pele_simulation_path,ligand_folder)
            receptor = [x for x in os.listdir(ligand_folder_path) if x.startswith('receptor')][0]
            ligand = [x for x in os.listdir(ligand_folder_path) if x != receptor][0]

            self._PDBMerger(receptor,ligand)

    def setRdockToPELESimulation(self):
        # Rdock is going to be difficoult due to deprotonation.
        pass

    def setGlideToPELESimulation(self):
        # Glide also has the handicap that the separation of files has to be done
        # with a machine with a SHCRODINGER's license.
        pass
    
    def setEquibindToPELESimulation(self):

        self._equibindDockingPoseRetriever()
        self._equibindPELEInputGenerator()
     