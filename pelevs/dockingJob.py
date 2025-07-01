from openbabel import pybel
from Bio.PDB import PDBParser
import os
import shutil
import MDAnalysis
from openbabel import openbabel as ob
import warnings

warnings.filterwarnings("ignore")


class DockingJob:
    """
    Attributes
    ==========
    receptor : str
        Name of the receptor's file (sdf, mol2 or pdb).
    ligands : str
        Name of the file were ligands in a csv with SMILES are located.
    docking_tool : str
        Docking software wanted.
    grid_file : str
        Name of Glide's grid file (.zip)
    reference_ligand : str
        Name of the bound ligand that is going to be used to generate the cavity
        and the grid with rDock.
    ligand_score : str
        Ligand we want to obtain a rescore with.
    output_models : int
        Number of conformations generated per docking.

    Methods
    =======
    setGlideDocking(self, grid_file, forcefield)
        Prepare a Glide docking folder with necessary files to send to a
        machine with SCHRÖDINGER's license and launch the job.
    setRdockDocking(self, reference_ligand, ligands, cpus_docking)
        Prepare an rDock docking folder with necessary files to send to MN4
        and launch the job.
    setEquibindDocking(self, ligands, receptor)
        Prepare an Equibind docking folder with necessary files to send to
        CTE-POWER and launch the job.
    rDockRescore(self, complete_structure)
        Prepare a rescoring job with rDock.
    GlideRescore(self, grid_file, ligand)
        Prepare a rescoring job with Glide.

    Hidden Methods
    ==============
    _ligandsChecker(self, ligands)
        Checks if the ligands inputted come from a
        ligprep job with sdf format
    _folderPreparation(self)
        Prepares folder paths for the new docking job.
    _glidePrepareJob(self, grid_file, forcefield)
        Prepare files for a docking with Glide
    _rdockRescorePreparation(self, complete_structure)
        Prepare the files to launch an rdock rescore
        simulation.
    _rdockReceptorFormatChecker(self, receptor)
        Check receptor's file format and change it to mol2.
    _rdockFileCopier(self, reference_ligand)
        Copy important files to their assignesd places.
    _rdockParamFilesWriter(self, receptor, reference_ligand)
        Generate param file for the cavity generation and docking
        with rDock.
    _rdockGridGenerator(self)
        Generate files to send cavity and grid files generator job.
    _rdockJobSplitter(self, ligands, cpus_docking)
        Split the total number of molecules into number of cpus available.
    _rdockRunFilesGenerator(self, cpus_docking)
        Generate all the necessary run files.
    _equibindReceptorFormatChecker(self, receptor)
        Check if formats for equibind docking are correct.
    _equibindSplitLigands(self, ligands)
        Splitting ligprep's output sdf file into individual sdfs.
    _equibindFolderPreparation(self, receptor)
        Preparing folder system for equibind docking simulations.
    _equibindFilesPreparation(self)
        Writing necessary inference.yml and run files.
    """

    def __init__(self):
        """
        Prepare files and folders coming from the liprep job in a
        directory to perform  acertain kind of docking.
        """

        self.prepared_ligands = False

        if os.path.isdir("1_input_files/receptor/"):
            receptor = (
                "1_input_files/receptor/" + os.listdir("1_input_files/receptor/")[0]
            )
        else:
            raise Exception(
                "MissingReceptorFile: Receptor file should be located at '1_input_files/receptor/'"
            )

        if os.path.isdir("2_ligprep_job/job/"):
            ligands = (
                "2_ligprep_job/job/"
                + [x for x in os.listdir("2_ligprep_job/job/") if x.endswith(".sdf")][0]
            )

        elif os.path.isdir("1_input_files/ligands"):
            files = os.listdir("1_input_files/ligands")
            if len(files) == 1 and files[0].split(".")[1] == "sdf":
                ligands = "1_input_files/ligands/" + files[0]
                self.prepared_ligands = True
                print(" - The prepared ligands are already prepared and in sdf format.")

        else:
            raise Exception(
                "MissingLigandsFile: Ligands file should be located at '2_ligprep_job/job/'"
            )

        self.receptor = receptor.split("/")[-1]
        self.ligands = ligands.split("/")[-1]
        self.docking_tool = None
        self.grid_file = None
        self.reference_ligand = None
        self.sphere = None
        self.ligand_score = None
        self.output_models = None

        self._ligandsChecker(ligands)
        self._folderPreparation()

    def _ligandsChecker(self, ligands):
        """
        Check validity of inputted ligands (sdf).

        Parameters
        ==========
        ligands : str
            Name of the csv file with SMILES and id.
        """

        if (ligands.split(".")[-1] != "sdf") or (ligands.split(".")[-1] != "sd"):
            pass
        else:
            raise Exception(
                "LigandsFileError: The ligands file shoud be the output from ligprep in sdf format."
            )

    def _folderPreparation(self):
        """
        Prepare docking folder and job folder to be sent
        to calculate.

        Parameters
        ==========
        ligands : str
            Name of the csv file with SMILES and id.
        """

        if not os.path.isdir("3_docking_job"):
            os.mkdir("3_docking_job")

        if not os.path.isdir("3_docking_job/job"):
            os.mkdir("3_docking_job/job")

    def _glidePrepareJob(self, grid_file, forcefield, protocol, output_model):
        """
        Copy files to job folder and generate necessary .in file
        to perform Glide simulation.

        Parameters
        ==========
        grid_file : str
            Name of the grid file (.zip) to dock ligands.
        forcefield : str
            Name of the forcefield to be used in the Glide docking.
        protocol : str
            Name of the protocol used to score the docking.
        output_model : int
            Number of output conformations in the docking.
        """

        input_ligand_path = "1_input_files/ligands/"
        docking_job_path = "3_docking_job/job"
        ligprep_path = "2_ligprep_job/job/"
        glide_score_path = "3_docking_job/glide_score"
        grid_file_name = os.path.basename(grid_file)

        if protocol == "dock":
            shutil.copy(grid_file, docking_job_path)

            if not self.prepared_ligands:
                shutil.copy(ligprep_path + self.ligands, docking_job_path)
            else:
                shutil.copy(input_ligand_path + self.ligands, docking_job_path)

            grid_path = os.path.join(docking_job_path, "glide_job.sh")
            job_path = os.path.join(docking_job_path, "glide_job.in")

        elif protocol == "score":

            ligand_name = os.path.basename(self.ligand_score).split(".")[0]
            ligand_score_path = os.path.join(glide_score_path, ligand_name)

            if not os.path.isdir(glide_score_path):
                os.mkdir(glide_score_path)

            if not os.path.isdir(ligand_score_path):
                os.mkdir(ligand_score_path)

            shutil.copy(grid_file, ligand_score_path)

            grid_path = os.path.join(ligand_score_path, "glide_score.sh")
            job_path = os.path.join(ligand_score_path, "glide_score.in")

        with open(grid_path, "w") as filein:
            if protocol == "dock":
                filein.writelines(
                    '"${SCHRODINGER}/glide" glide_job.in -OVERWRITE -adjust -HOST localhost:1 -TMPLAUNCHDIR'
                )
            elif protocol == "score":
                filein.writelines(
                    '"${SCHRODINGER}/glide" glide_score.in -OVERWRITE -adjust -HOST localhost:1 -TMPLAUNCHDIR'
                )

        with open(job_path, "w") as filein:
            filein.writelines(
                ["GRIDFILE   {}\n".format(grid_file_name), "PRECISION   SP\n"]
            )

            if protocol == "dock":
                filein.writelines(
                    [
                        "LIGANDFILE   {}\n".format(self.ligands),
                        "FORCEFIELD   {}\n".format(forcefield),
                        "POSES_PER_LIG   {}\n".format(output_model),
                        "POSTDOCK_NPOSE   {}\n".format(output_model),
                    ]
                )

            if protocol == "score":
                filein.writelines(
                    [
                        "LIGANDFILE   {}\n".format(os.path.basename(self.ligand_score)),
                        "DOCKING_METHOD   inplace\n",
                        "POSTDOCK   False\n",
                    ]
                )

        if protocol == "score":
            shutil.copy(self.ligand_score, ligand_score_path)

        print(
            " - Glide job generated successfully with grid {grid} and forcefield {ff}.".format(
                grid=grid_file, ff=forcefield
            )
        )

    def _rdockRescorePreparation(self, complete_structure):
        """
        Prepare the files for rDock Rescoring by splitting the complete structure into receptor and ligand

        Parameters
        ==========
        complete_structure : str
            File name of the complexed structure
        """

        def split_pdb(complete_structure, receptor_file, ligand_file, reference_ligand):
            """
            Splits the complex file into ligand, reference ligand, and receptor.

            Parameters
            ==========
            complete_structure : str
                File name of the complexed structure.
            receptor_file : str
                File name of the receptor structure.
            ligand_file : str
                File name of the ligand structure.
            reference_ligand : str
                File name of the reference ligand structure.
            """

            # Retrieve chains
            pdb = PDBParser().get_structure("complex", complete_structure)

            for chain in pdb.get_chains():
                if chain.id == "L":
                    ligand_chain = chain.id
                else:
                    receptor_chain = chain.id

            # Import the pdb file as universe
            u = MDAnalysis.Universe(complete_structure)

            # Select with chains
            receptor = u.select_atoms("chainID {}".format(receptor_chain))
            ligand = u.select_atoms("chainID {}".format(ligand_chain))

            # Save separated files as pdb
            receptor.write(receptor_file)
            ligand.write(ligand_file)
            ligand.write(reference_ligand)

            self.ligands = os.path.basename(ligand_file.split(".")[0] + ".sdf")
            self.reference_ligand = os.path.basename(ligand_file.split(".")[0] + ".sdf")
            self.receptor = os.path.basename(receptor_file.split(".")[0] + ".mol2")

        def conversor(file_in, format_out, path_out):
            """
            Converts the input file to PDB format using Open Babel.

            Parameters
            ==========
            file_in : str
                Name of the file we want to convert.
            format_out : str
                Format we want to convert to.
            path_out : str
                Path to the directory we want the converted file to keep.
            """

            file = os.path.basename(file_in)
            file_name, file_format = file.split(".")

            conv = ob.OBConversion()
            conv.SetInAndOutFormats(file_format, format_out)
            mol = ob.OBMol()
            conv.ReadFile(mol, file_in)
            conv.WriteFile(
                mol,
                os.path.join(
                    path_out,
                    "{name}.{format}".format(name=file_name, format=format_out),
                ),
            )

        rdock_path = "3_docking_job/rdock_score"

        if not os.path.isdir(rdock_path):
            os.mkdir(rdock_path)

        file_structure = os.path.basename(complete_structure).split(".")[0]
        path_ligand = os.path.join(rdock_path, file_structure)

        if not os.path.isdir(path_ligand):
            os.mkdir(path_ligand)

        shutil.copy(complete_structure, path_ligand)

        # Splitting file name
        input_pdb_file = os.path.join(path_ligand, "{}.pdb".format(file_structure))
        receptor_pdb_file = "{}/receptor.pdb".format(path_ligand)
        ligand_pdb_file = "{}/ligand.pdb".format(path_ligand)
        reference_ligand_pdb_file = "{}/reference_ligand.pdb".format(path_ligand)

        split_pdb(
            input_pdb_file,
            receptor_pdb_file,
            ligand_pdb_file,
            reference_ligand_pdb_file,
        )

        # Converting
        conversor(receptor_pdb_file, "mol2", path_ligand)
        conversor(ligand_pdb_file, "sdf", path_ligand)
        conversor(reference_ligand_pdb_file, "sdf", path_ligand)

        # Removing
        os.remove(receptor_pdb_file)
        os.remove(ligand_pdb_file)
        os.remove(reference_ligand_pdb_file)

        self.ligand_score = file_structure

    def _rdockReceptorFormatChecker(self, receptor):
        """
        Check receptor's format and transform it if necessary to mol2.

        Parameters
        ==========
        receptor : str
            File name of the receptor.
        """

        if len(receptor.split(".")) == 2:
            receptor_name, receptor_format = receptor.split(".")

        else:
            raise Exception(
                "ReceptorNameError: Receptor name should only have one dot separating the name and the format."
            )

        if receptor_format != "mol2":

            print(" - Changing receptor's format to mol2.")
            receptor_generator = pybel.readfile(
                receptor_format, "1_input_files/receptor/" + receptor
            )
            receptor_molecule = next(receptor_generator)
            receptor_molecule.write(
                "mol2",
                "3_docking_job/job/{}.mol2".format(receptor_name),
                overwrite=True,
            )

        else:
            pass

    def _rdockFileCopier(self, reference_ligand):
        """
        Copy files to job folder.

        Parameters
        ==========
        reference_ligand : str
            File name of the ligand used to generate cavity
            and grid for rDock.
        """

        os.makedirs("3_docking_job/job/ligands", exist_ok=True)

        if bool(reference_ligand):

            reference_ligand_name = os.path.basename(reference_ligand)

            if os.path.isfile(reference_ligand):
                shutil.copy(
                    reference_ligand,
                    os.path.join("1_input_files/ligands", reference_ligand_name),
                )
                shutil.move(
                    reference_ligand,
                    os.path.join("3_docking_job/job", reference_ligand_name),
                )
            elif not os.path.isfile(reference_ligand):
                if os.path.isfile(
                    os.path.join("3_docking_job/job", reference_ligand_name)
                ):
                    pass
                elif os.path.isfile(
                    os.path.join("1_input_files/ligands", reference_ligand_name)
                ):
                    shutil.copy(
                        os.path.join("1_input_files/ligands", reference_ligand_name),
                        os.path.join("3_docking_job/job", reference_ligand_name),
                    )
                else:
                    raise Exception(
                        "ReferenceLigandError: Tha path is not correct and it has not been found in 3_docking_job/job."
                    )

        if self.prepared_ligands:
            shutil.copy(
                "1_input_files/ligands/" + self.ligands_to_dock,
                "3_docking_job/job/" + self.ligands_to_dock,
            )
        else:
            shutil.copy("2_ligprep_job/job/" + self.ligands, "3_docking_job/job")

        shutil.copy(
            "1_input_files/receptor/" + self.receptor,
            "3_docking_job/job/" + self.receptor,
        )

    def _rdockParamFilesWriter(self, receptor, reference_ligand, sphere, protocol):
        """
        Writing parameter file to run rDock simulations and cavity
        generation.

        Parameters
        ==========
        receptor : str
            File name of the receptor.
        reference_ligand : str
            File name of the ligand used to generate cavity
            and grid for rDock.
        protocol : str
            Name of the protocol used to score the ligand(s).
        """

        receptor_mol2 = "{}.mol2".format(receptor.split(".")[0])

        if protocol == "dock":
            parameter_file = os.path.join("3_docking_job/job", "parameter_file.prm")

        elif protocol == "score":
            parameter_file = os.path.join(
                "3_docking_job/rdock_score/{}".format(self.ligand_score),
                "parameter_file.prm",
            )

        if not os.path.isfile(parameter_file):
            with open(parameter_file, "w") as fileout:
                if bool(reference_ligand):
                    fileout.writelines(
                        "RBT_PARAMETER_FILE_V1.00\n"
                        "TITLE rdock\n"
                        "\n"
                        "RECEPTOR_FILE " + receptor_mol2 + "\n"
                        "RECEPTOR_FLEX 3.0\n"
                        "\n"
                        "##################################################################\n"
                        "### CAVITY DEFINITION:\n"
                        "##################################################################\n"
                        "SECTION MAPPER\n"
                        "    SITE_MAPPER RbtLigandSiteMapper\n"
                        "    REF_MOL " + reference_ligand + "\n"
                        "    RADIUS 6.0\n"
                        "    SMALL_SPHERE 1.0\n"
                        "    MIN_VOLUME 100\n"
                        "    MAX_CAVITIES 1\n"
                        "    VOL_INCR 0.0\n"
                        "   GRIDSTEP 0.5\n"
                        "END_SECTION\n"
                        "\n"
                        "#################################\n"
                        "#CAVITY RESTRAINT PENALTY\n"
                        "#################################\n"
                        "SECTION CAVITY\n"
                        "    SCORING_FUNCTION RbtCavityGridSF\n"
                        "    WEIGHT 1.0\n"
                        "END_SECTION\n"
                    )
                elif bool(sphere):
                    sphere_center, sphere_radius = sphere
                    fileout.writelines(
                        "RBT_PARAMETER_FILE_V1.00\n"
                        "TITLE rdock\n"
                        "\n"
                        "RECEPTOR_FILE " + receptor_mol2 + "\n"
                        "RECEPTOR_FLEX 3.0\n"
                        "\n"
                        "##################################################################\n"
                        "### CAVITY DEFINITION:\n"
                        "##################################################################\n"
                        "SECTION MAPPER\n"
                        "    SITE_MAPPER RbtSphereSiteMapper\n"
                        "    CENTER " + sphere_center + "\n"
                        "    RADIUS " + str(sphere_radius) + "\n"
                        "    SMALL_SPHERE 1.0\n"
                        "    MIN_VOLUME 100\n"
                        "    MAX_CAVITIES 1\n"
                        "    VOL_INCR 0.0\n"
                        "   GRIDSTEP 0.5\n"
                        "END_SECTION\n"
                        "\n"
                        "#################################\n"
                        "#CAVITY RESTRAINT PENALTY\n"
                        "#################################\n"
                        "SECTION CAVITY\n"
                        "    SCORING_FUNCTION RbtCavityGridSF\n"
                        "    WEIGHT 1.0\n"
                        "END_SECTION\n"
                    )

    def _rdockGridGenerator(self, protocol):
        """
        Generate run file to generate cavity and grid for rDock.

        Parameters
        ==========
        protocol : str
            Name of the protocol used to score the ligand(s).
        """
        if protocol == "dock":
            with open("3_docking_job/job/grid.sh", "w") as fileout:
                fileout.writelines(
                    "module load rdock\n"
                    "rbcavity -was -d -r parameter_file.prm > parameter_file.log\n"
                )

        elif protocol == "score":
            with open(
                "3_docking_job/rdock_score/{}/grid.sh".format(self.ligand_score), "w"
            ) as fileout:
                fileout.writelines(
                    "module load rdock\n"
                    "rbcavity -was -d -r parameter_file.prm > parameter_file.log\n"
                )

    def _rdockJobSplitter(self, ligands, cpus_docking, protocol):
        """
        Generate script and run files to split the sdf with N
        molecules between M cpus.

        Parameters
        ==========
        ligands : str
            Name of the outputted ligands by ligprep.
        cpus_docking : int
            Number of cpus to use for the rDock docking.

        """

        print(
            " - Splitting {ligands_file}'s molecules into {cpus} (different) file(s).".format(
                ligands_file=ligands, cpus=cpus_docking
            )
        )

        if protocol == "dock":

            # Generating split file
            if not os.path.isfile("3_docking_job/job/splitMols.sh"):
                with open("3_docking_job/job/splitMols.sh", "w") as filein:
                    filein.writelines(
                        "#!/bin/bash\n"
                        "#Usage: splitMols.sh <input> #Nfiles <outputRoot>\n"
                        "fname=$1\n"
                        "nfiles=$2\n"
                        "output=$3\n"
                        "molnum=$(grep -c '$$$$' $fname)\n"
                        'echo " - $molnum molecules found"\n'
                        "echo \" - Dividing '$fname' into $nfiles files\"\n"
                        'echo " "\n'
                        "rawstep=`echo $molnum/$nfiles | bc -l`\n"
                        "let step=$molnum/$nfiles\n"
                        "if [ ! `echo $rawstep%1 | bc` == 0 ]; then\n"
                        "        let step=$step+1;\n"
                        "fi;\n"
                        "sdsplit -$step -o$output $1\n"
                    )

            # Generating splitted ligand files
            with open("3_docking_job/job/split.sh", "w") as fileout:
                fileout.writelines(
                    "bash splitMols.sh {ligands_file} {cpus} ligands/split\n".format(
                        ligands_file=ligands, cpus=cpus_docking
                    )
                )

        if protocol == "score":

            path_splitMols = "3_docking_job/rdock_score/{}/splitMols.sh".format(
                self.ligand_score
            )
            path_split = "3_docking_job/rdock_score/{}/split.sh".format(
                self.ligand_score
            )

            if not os.path.isfile(path_splitMols):
                with open(path_splitMols, "w") as filein:

                    filein.writelines(
                        "#!/bin/bash\n"
                        "#Usage: splitMols.sh <input> #Nfiles <outputRoot>\n"
                        "fname=$1\n"
                        "nfiles=$2\n"
                        "output=$3\n"
                        "molnum=$(grep -c '$$$$' $fname)\n"
                        'echo " - $molnum molecules found"\n'
                        "echo \" - Dividing '$fname' into $nfiles files\"\n"
                        'echo " "\n'
                        "rawstep=`echo $molnum/$nfiles | bc -l`\n"
                        "let step=$molnum/$nfiles\n"
                        "if [ ! `echo $rawstep%1 | bc` == 0 ]; then\n"
                        "        let step=$step+1;\n"
                        "fi;\n"
                        "sdsplit -$step -o$output $1\n"
                    )

                # Generating splitted ligand files
                with open(path_split, "w") as fileout:
                    fileout.writelines(
                        "bash splitMols.sh {ligands_file} {cpus} split\n".format(
                            ligands_file=ligands, cpus=cpus_docking
                        )
                    )

    def _rdockRunFilesGenerator(
        self, cpus_docking, protocol, output_models, queue, time
    ):
        """
        Generate all the necessary runs to run an individual
        rDock simulation, to prepare the rDock simulation, and
        to run all the individual simulations.

        Parameters
        ==========
        cpus_docking : str
            Number of cpus to be used in the docking.
        reference_ligand : str
            File name of the ligand used to generate cavity
            and grid for rDock.
        protocol : str
            Name of the protocol to be used, it can be either dock or score
        output_models : int
            Number of output conformations in the docking.
        """

        # Generating folders
        if not os.path.isdir("3_docking_job/job/ligands"):
            os.makedirs("3_docking_job/job/ligands", exist_ok=True)

        if not os.path.isdir("3_docking_job/job/results"):
            os.makedirs("3_docking_job/job/results", exist_ok=True)

        if protocol == "dock":
            # Generating run files
            for i in range(1, cpus_docking + 1):
                with open("3_docking_job/job/run{}".format(i), "w") as fileout:
                    fileout.writelines(
                        "#!/bin/sh\n"
                        "#SBATCH --job-name=rdock" + str(i) + " \n"
                        "#SBATCH --time=" + time + "\n"
                        "#SBATCH --ntasks=1\n"
                        "#SBATCH --output=rdock.out\n"
                        "#SBATCH --error=rdock.err\n"
                        "#SBATCH --qos=" + queue + "\n"
                        "#SBATCH --account=bsc72\n"
                        "\n"
                        "module load rdock\n"
                        "\n"
                        "rbdock -i ligands/split{val}.sd -o results/split{val}_out -r parameter_file.prm -p dock.prm -n {out} -allH\n".format(
                            val=i, out=output_models
                        )
                    )

            with open("3_docking_job/job/prepare_rDock_run.sh", "w") as fileout:
                fileout.writelines(
                    "#!/bin/bash\n"
                    "#SBATCH --job-name=rdock1\n"
                    "#SBATCH --time=00:30:00\n"
                    "#SBATCH --ntasks=1\n"
                    "#SBATCH --output=rdock.out\n"
                    "#SBATCH --error=rdock.err\n"
                    "#SBATCH --qos=gp_debug\n"
                    "#SBATCH --account=bsc72\n"
                    "\n"
                    "# Run grid.sh\n\n"
                    "echo ' '\n"
                    "echo ' - Generating grid and cavity'\n"
                    "echo ' - Loading rDock module:'\n"
                    "echo ' '\n"
                    "source grid.sh\n"
                    "echo ' '\n\n"
                    "# Run split.sh\n"
                    "echo ' - Splitting ligands'\n"
                    "source split.sh\n"
                )

            with open("3_docking_job/job/rDock_run.sh", "w") as fileout:
                fileout.writelines(
                    "#!/bin/bash\n" "for d in run*; do echo ${d}; sbatch ${d}; done"
                )

        elif protocol == "score":
            # Generating run files
            for i in range(1, cpus_docking + 1):
                with open(
                    "3_docking_job/rdock_score/{}/run{}".format(self.ligand_score, i),
                    "w",
                ) as fileout:
                    fileout.writelines(
                        "#!/bin/sh\n"
                        "#SBATCH --job-name=rdock" + str(i) + " \n"
                        "#SBATCH --time=" + time + "\n"
                        "#SBATCH --ntasks=1\n"
                        "#SBATCH --output=rdock.out\n"
                        "#SBATCH --error=rdock.err\n"
                        "#SBATCH --qos=" + queue + "\n"
                        "#SBATCH --account=bsc72\n"
                        "\n"
                        "module load rdock/rdock-24.04.204-legacy-gcc\n"
                        "module load ANACONDA/2019.10\n"
                        "module load intel\n"
                        "module load mkl\n"
                        "module load impi\n"
                        "module load gcc\n"
                        "module load boost/1.64.0\n"
                        "\n"
                        "\n"
                        "rbdock -i ligand.sdf -o ligand_out -r parameter_file.prm -p score.prm -allH\n".format(
                            val=i
                        )
                    )

            with open(
                "3_docking_job/rdock_score/{}/prepare_rDock_run.sh".format(
                    self.ligand_score
                ),
                "w",
            ) as fileout:
                fileout.writelines(
                    "#!/bin/bash\n"
                    "# Run grid.sh\n\n"
                    "echo ' '\n"
                    "echo ' - Generating grid and cavity'\n"
                    "echo ' - Loading rDock module:'\n"
                    "echo ' '\n"
                    "source grid.sh\n"
                    "echo ' '\n\n"
                    "# Run split.sh\n"
                    "echo ' - Splitting ligands'\n"
                    "source split.sh\n"
                )

            with open(
                "3_docking_job/rdock_score/{}/rDock_run.sh".format(self.ligand_score),
                "w",
            ) as fileout:
                fileout.writelines(
                    "#!/bin/bash\n" "for d in run*; do echo ${d}; sbatch ${d}; done"
                )

        print(" - Job generated to be sent to MN4 machine.")
        print(
            " - RDock docking job generated successfully to run with {} cpu(s).".format(
                cpus_docking
            )
        )
        print(
            " - Once in the MN4, first run: bash prepare_rDock_run and after that\n   run: bash rDock_run.sh"
        )

    def _equibindReceptorFormatChecker(self, receptor):
        """
        Check receptor's format and transform it if necessary to pdb.

        Parameters
        ==========
        receptor : str
            File name of the receptor.
        """

        if len(receptor.split(".")) == 2:
            receptor_name, receptor_format = receptor.split(".")

        else:
            raise Exception(
                "ReceptorNameError: Receptor name should only have one dot separating the name and the format."
            )

        if receptor_format != "pdb":
            print(" - Changing receptor's format to pdb.")

            receptor_generator = pybel.readfile(
                receptor_format, "1_input_files/receptor/" + receptor
            )
            receptor_molecule = next(receptor_generator)
            receptor_molecule.write(
                "pdb", "3_docking_job/{}.pdb".format(receptor_name), overwrite=True
            )

        else:
            shutil.copy("1_input_files/receptor/" + receptor, "3_docking_job")

    def _equibindSplitLigands(self, ligands):
        """
        Split ligprep's output sdf into individual molecules to
        perform an equibind docking on all of them.

        Parameters
        ==========
        ligands : str
            File name of the ligprep's output.
        """

        with open("2_ligprep_job/job/{}".format(ligands), "r") as f:
            content = f.read().strip()

        records = content.split("$$$$")

        for record in records:
            if record.strip():

                # Extract the value of the property for naming the output file
                variant_value = None
                lines = record.split("\n")
                for line in lines:
                    if "<s_lp_Variant>" in line:
                        variant_value = lines[lines.index(line) + 1]
                        break

                if variant_value:
                    output_file = variant_value + ".sdf"

                    if not os.path.isdir("3_docking_job/job/equibind_calculations"):
                        os.mkdir("3_docking_job/job/equibind_calculations")

                    if not os.path.isdir("3_docking_job/job/equibind_results"):
                        os.mkdir("3_docking_job/job/equibind_results")

                    if not os.path.isdir(
                        "3_docking_job/job/equibind_calculations/{}".format(
                            variant_value
                        )
                    ):
                        os.mkdir(
                            "3_docking_job/job/equibind_calculations/{}".format(
                                variant_value
                            )
                        )

                    # Write the record to the output file
                    with open(
                        "3_docking_job/job/equibind_calculations/{folder}/{output}".format(
                            folder=variant_value, output=output_file
                        ),
                        "w",
                    ) as f:
                        if record.startswith("\n"):
                            record = record[1:]

                        f.write(record + "$$$$")

    def _equibindFolderPreparation(self, receptor):
        """
        Prepare the folders to send an equibind docking job.

        Parameters
        ==========
        receptor : str
            File name of the receptor.
        """

        for folder in os.listdir("3_docking_job/job/equibind_calculations"):
            folder_path = os.path.join(
                "3_docking_job/job/equibind_calculations", folder
            )

            if os.path.isdir(folder_path):
                destination_protein = os.path.join(folder_path, folder + "_protein.pdb")
                destination_ligands = os.path.join(folder_path, folder + "_ligand.sdf")

                shutil.copyfile(
                    "3_docking_job/{}".format(receptor), destination_protein
                )
                os.rename(
                    "3_docking_job/job/equibind_calculations/{name}/{name}.sdf".format(
                        name=folder
                    ),
                    destination_ligands,
                )

    def _equibindFilesPreparation(self):
        """
        Write run file as well as inference.yml to be able to
        send an equibind docking job.
        """

        with open("3_docking_job/job/run", "w") as fileout:
            fileout.writelines(
                "#!/bin/bash\n"
                "#SBATCH --job-name=equi\n"
                "#SBATCH --time=1:00:00\n"
                "#SBATCH --gres gpu:1\n"
                "#SBATCH --cpus-per-task=40\n"
                "#SBATCH --ntasks=1\n"
                "#SBATCH --output=equi.out\n"
                "#SBATCH --error=equi.err\n"
                "\n"
                "module purge\n"
                "module load anaconda3/2020.02\n"
                "module list\n"
                "\n"
                'eval "$(conda shell.bash hook)"\n'
                "conda activate /apps/ANACONDA3/2020.02/envs/ESM-EquiBind-DiffDock\n"
                "\n"
                "cp -r /gpfs/apps/POWER9/ANACONDA3/2020.02/envs/ESM-EquiBind-DiffDock/modules/EquiBind/runs .\n"
                "python /gpfs/apps/POWER9/ANACONDA3/2020.02/envs/ESM-EquiBind-DiffDock/modules/EquiBind/inference.py --config=inference.yml\n"
            )

        with open("3_docking_job/job/inference.yml", "w") as fileout:
            fileout.writelines(
                "run_dirs:\n"
                '  - flexible_self_docking # the resulting coordinates will be saved here as tensors in a .pt file (but also as .sdf files if you specify an "output_directory" below)\n'
                "inference_path: 'equibind_calculations' # this should be your input file path as described in the main readme\n"
                "\n"
                "test_names: timesplit_test\n"
                "output_directory: 'equibind_results' # the predicted ligands will be saved as .sdf file here\n"
                "run_corrections: True\n"
                "use_rdkit_coords: False # generates the coordinates of the ligand with rdkit instead of using the provided conformer. If you already have a 3D structure that you want to use as initial conformer, then leave this as False\n"
                "save_trajectories: False\n"
                "\n"
                "num_confs: 1 # usually this should be 1\n"
            )

        print(" - Job created to be sent to CTE-POWER")
        print(" - Equibind docking job created successfully.")

    def setGlideDocking(self, grid_file, forcefield="OPLS_2005", output_models=50):
        """
        Prepare the job folder to send a Glide docking job.

        Parameters
        ==========
        grid_file : str
            Name of the grid file (.zip) to dock ligands.
        forcefield : str
            Name of the forcefield to be used in the Glide docking.
        output_models : int
            Number of output conformations.
        """

        self.grid_file = grid_file
        self.docking_tool = "glide"
        self.protocol = "dock"
        self.output_models = output_models

        protocol = self.protocol

        self._glidePrepareJob(grid_file, forcefield, protocol, output_models)

    def setRdockDocking(
        self,
        ligands,
        cpus_docking,
        reference_ligand=None,
        sphere=None,
        output_models=50,
        queue="gp_debug",
        time="01:00:00",
    ):
        """
        Prepare the job folder to send an rDock docking job.

        Parameters
        ==========
        ligands : str
            Name of the outputted ligands by ligprep.
        cpus_docking : int
            Number of cpus to use for the rDock docking.
        reference_ligand : str
            File name of the ligand used to generate cavity
            and grid for rDock.
        sphere : list[str, int]
            List with a string with the coordinates of the sphere's center and the radius in Angstrom (i.e. ['(0,0,0)', 10]).
        output_models : int
            Number of output conformations.
        queue : str
            Queue of the machine you want to put in the run files (bsc_ls or debug).
        time : str
            Time you request to use the machine you send the job to. The format is:
            hours:minutes:seconds.
        """

        self.ligands_to_dock = os.path.basename(ligands)
        self.reference_ligand = reference_ligand
        self.sphere = sphere
        self.docking_tool = "rdock"
        protocol = "dock"
        self.output_models = output_models

        if not (bool(reference_ligand) ^ bool(sphere)):  # XOR condition
            raise ValueError(
                "You must provide either `sphere` or `reference_ligand`, but not both."
            )

        # Function implementation goes here
        if reference_ligand:
            print(f"Using reference ligand method")
        else:
            print(" - Using sphere method")

        self._rdockReceptorFormatChecker(self.receptor)
        self._rdockFileCopier(reference_ligand)

        self._rdockParamFilesWriter(
            self.receptor, self.reference_ligand, self.sphere, protocol
        )
        self._rdockGridGenerator(protocol)
        self._rdockJobSplitter(ligands, cpus_docking, protocol)
        self._rdockRunFilesGenerator(cpus_docking, protocol, output_models, queue, time)

    def setEquibindDocking(self, ligands, receptor):
        """
        Prepare the job folder to send an Equibind docking job.

        Parameters
        ==========
        ligands : str
            Name of the outputted ligands by ligprep.
        receptor : str
            File name of the receptor.
        """

        self.docking_tool = "equibind"

        self._equibindReceptorFormatChecker(receptor)
        self._equibindSplitLigands(ligands)
        self._equibindFolderPreparation(receptor)
        self._equibindFilesPreparation()

    def rdockRescore(self, complete_structure):
        """
        Prepare the job folder to send an rDock rescoring job.

        Parameters
        ==========
        complete_structure : str
            File of the structure to rescore.
        """

        protocol = "score"
        cpus_docking = 1
        self.docking_tool = "rdock"

        self._rdockRescorePreparation(complete_structure)

        self._rdockParamFilesWriter(self.receptor, self.reference_ligand, protocol)

        self._rdockGridGenerator(protocol)
        self._rdockJobSplitter(self.ligands, cpus_docking, protocol)
        self._rdockRunFilesGenerator(cpus_docking, protocol, None)

    def glideRescore(self, grid_file, ligand):
        """
        Prepare the job folder to send a Glide rescoring job

        Parameters
        ==========
        grid_file : str
            Name of the grid file (.zip) to dock ligands.
        ligand : str
            Ligand to rescore
        """

        protocol = "score"
        self.docking_tool = "glide"
        self.ligand_score = ligand
        forcefield = "OPLS2005"

        self._glidePrepareJob(grid_file, forcefield, protocol)
