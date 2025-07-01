import os
import pandas as pd
import shutil
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import re
from glob import glob


class PELEAnalyzer:
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
    experimental_data : pandas.DataFrame
        Dataframe with the csv with experimental data
    calculated_data : pandas.DataFrame
        Dataframe with all the information form the simulations
    _path_simulation : str
        Path to store the simulations.
    _path_analysis : str
        Path to store the analysis made with the simulations.
    all_data : pandas.DataFrame
        Dataframe with calculated and experimental data.
    equibind_data : pandas.DataFrame
        Dataframe with the selection of best poses according to
        the binding energy minimum.

    Methods
    =======
    experimentalDataCollector(self, path_experimental_data=None)
        Retrieve the lready inputed experimental data or manage newly inputed experimental data.
    equibindDataTrimming(self, df)
        Manage all the data from all the tautomer simulations to only keep thos with best PELE score.
    PELEDataCollector(self)
        Collect all the data from te simulations: options chosen and the scores calculated for all the
        metrics.
    correlationPlotter(self, x_label, y_label, sampling, df=None)
        Plot correlations between the y_label and x_label according to the data in df.
    simulationAnalyzer(self, protocol, system)
        Generates plots of all the interesting combinations of metrics in the reports of
        concrete simulations.

    Hidden Methods
    ==============
    _PELEDownloadsAssemble(self)
        Prepare necessary folders to store data and move the data to them.
    _PELEFoldersToAnalyze(self)
        Determine which directories have to be checked for simulations.
    _energyCalculator(self, dataset_location, system, path_system)
        Calculate all the scores of all the metric of all the simulations performed for
        each rescoring method.
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
        self.experimental_data = None
        self.calculated_data = None
        self._path_simulation = None
        self._path_analysis = None
        self.all_data = None
        self.equibind_data = None

    def _PELEDownloadsAssemble(self):
        """
        Assembles all the data downloaded into a hierarchical folder
        structure.
        """

        # Generating paths
        docking_tools = ["glide", "rdock", "equibind"]
        path_pele_analysis = "5_pele_analysis"
        path_simulations = os.path.join(path_pele_analysis, "simulations")
        path_analysis = os.path.join(path_pele_analysis, "analysis")

        # Generating folders
        if not os.path.isdir(path_pele_analysis):
            os.mkdir(path_pele_analysis)

        if not os.path.isdir(path_analysis):
            os.mkdir(path_analysis)

        if not os.path.isdir(path_simulations):
            os.mkdir(path_simulations)

        methods = []

        # Moving the data
        for folder in [
            x for x in os.listdir(".") if (os.path.isdir(x) and x in docking_tools)
        ]:
            methods.append(folder)
            shutil.move(folder, path_simulations)

        self.docking_tool = methods
        self._path_analysis = path_analysis
        self._path_simulation = path_simulations

    def _PELEFoldersToAnalyze(self):
        """
        Determines which folders with a rescoring method name in it have to be
        taken into account to retrieve all the data from the reports.

        Returns
        =======
        folders_to_check : list
            List of paths with the directories to be checked.
        """

        path_simulation = self._path_simulation
        path_analysis = self._path_analysis

        folders_to_check = []
        rescorings = ["xshort", "short", "long", "xlong", "xxlong"]

        for root, _, _ in os.walk(path_simulation):

            if os.path.basename(root) in rescorings:
                new_path = root.replace(path_simulation, path_analysis)
                folders_to_check.append(root)

                if not os.path.isdir(new_path):
                    os.makedirs(new_path)

        return folders_to_check

    def _energyCalculator(self, dataset_location, system, path_system, sample):
        """
        First, retrieves the characteristics of all the simulations performed (docking protocol, forcefield, protein part,
        perturbation level and rescoring). Second, for each individual simulation calculates scores of the different
        metrics calculated in the PELE simulation. Third, unites all the data into a single dataframe to work with.

        Scores calculated per simulation:
        1. Binding Energy
            be_min -> binding energy minimum
            be_bz -> binding energy boltzmann weighted
            be_p5 -> binding energy 5th percentile
            be_p10 -> binding energy 10th percentile
            be_p25 -> binding energy 25th percentile

        2. Total energy
            te_min -> total energy mimimum
            te_bz -> total energy boltzmann weighted
            te_p5 -> total energy 5th percentile
            te_p10 -> total energy 10th percentile
            te_p25 -> total energy 25th percentile

        3. SASA
            sasa_min -> sasa minimum
            sasa_bz -> sasa boltzmann weighted
            sasa_av -> sasa average
            sasa_max -> sasa maximum

        4. RMSD
            rmsd_max -> rmsd maximum
            rmsd_av -> rmsd average

        Parameters
        ==========
        dataset_location : str
            Slice of the path of the whole dataset to retrieve all the characteristics of
            the simulation.
        system : str
            Name of the ligand containing the simulation that we are interested in.
        path_system : str
            Path to the simulation inside a dataset we are interested in.

        Returns
        ==========
        simulation_data_dict : dict
            Dictionary with all the data of all the and all the scores calculated.
        """

        def _PELEOptionsRetriever(dataset_location):
            """
            Retrieves the characteristics of all the simulations performed (docking protocol, forcefield, protein part,
            perturbation level and rescoring).

            Parameters
            ==========
            dataset_location : str
                Slice of the path of the whole dataset to retrieve all the characteristics of
                the simulation.

            Returns
            ==========
            docking_param : str
                Docking tool used to generate initial positions for the simulation
            forcefield_param : str
                Force field used in the simulation.
            truncated_param : str
                Protein part used in the simulation.
            perturbation_param : str
                Perturbation protocol used in the simulation.
            rescorings_param : str
                Sampling method used in the simulation.
            """

            # PELE options
            docking_tools = ["glide", "rdock", "equibind"]
            forcefield = ["opls", "openff"]
            truncated = ["full", "truncated"]
            perturbation = ["minimization", "if", "refinement"]
            rescorings = ["xshort", "short", "long", "xlong", "xxlong"]

            location = dataset_location.split("/")

            # Initializing parameters to store
            docking_param = None
            forcefield_param = "opls"
            truncated_param = "truncated"
            perturbation_param = "refinement"
            rescorings_param = None

            for directory in location:

                if directory in docking_tools:
                    docking_param = directory
                elif directory in forcefield:
                    forcefield_param = directory
                elif directory in truncated:
                    truncated_param = directory
                elif directory in perturbation:
                    perturbation_param = directory
                elif directory in rescorings:
                    rescorings_param = directory

            self.sampling = rescorings_param
            self.protein_portion = truncated_param
            self.forcefield = forcefield_param
            self.perturbation_protocol = perturbation_param

            return (
                docking_param,
                forcefield_param,
                truncated_param,
                perturbation_param,
                rescorings_param,
            )

        def _energyInSimulation(path_system, sample):
            """
            For each individual simulation calculates scores of the different
            metrics calculated in the PELE simulation.

            Parameters
            ==========
            path_system : str
                Path to the simulation inside a dataset we are interested in.

            Returns
            ==========
            be_list : list
                List with the different scores calculated for the binding energy.
            te_list : list
                List with the different scores calculated for the total energy.
            sasa_list : list
                List with the different scores calculated for the SASA.
            rmsd_list : list
                List with the different scores calculated for the RMSD.
            """

            def minimum(vector):
                return min(vector)

            def boltzmann(vector, te, T):
                """
                Function
                ----------
                Calculates boltzmann weighted energy.

                Parameters
                ----------
                - vector : list
                    Binding energies of all the simulation.
                - te : list
                    Total energies of all the simulation.
                - T : float
                    Temperature to perform the Boltzmann weights with.

                Returns
                ----------
                - ene_bz : float
                    Value of the boltzmann weighted energy.
                """

                R = 1.985e-3  # Gas constant in kcal/mol
                vector = np.array(vector)

                exp_bz = np.exp(-te / (R * T))
                nominator = vector.dot(exp_bz)
                denominator = np.sum(exp_bz)
                ene_bz = nominator / denominator

                return ene_bz

            def p5(vector):
                return np.average(np.percentile(vector, 5))

            def p10(vector):
                return np.average(np.percentile(vector, 10))

            def p25(vector):
                return np.average(np.percentile(vector, 25))

            def maximum(vector):
                return max(vector)

            def average(vector):
                return np.average(vector)

            te = []
            be = []
            sasa = []
            rmsd = []

            list_epochs_og = [
                x
                for x in os.listdir(path_system)
                if os.path.isdir(os.path.join(path_system, x))
            ]

            list_epochs = [str(x) for x in sorted([int(x) for x in list_epochs_og])]

            if sample == "all":
                pass
            else:
                try:
                    if sample == "xshort" or sample == "short":
                        list_epochs = list_epochs[0]
                    elif sample == "long" or sample == "xlong":
                        list_epochs = list_epochs[0:6]
                except IndexError:
                    print(
                        "SampleError: The subsampling introduced is not compatible with the simulation carried out.\
                             Please check that you have chosen a smaller sampling than the one in the simulation."
                    )

            if len(list_epochs) != 0:
                for epoch in list_epochs:
                    path_epoch = os.path.join(path_system, epoch)
                    list_reports_og = [
                        x for x in os.listdir(path_epoch) if x.startswith("report_")
                    ]
                    list_report = sorted(
                        list_reports_og, key=lambda x: int(x.split("_")[1])
                    )

                    if sample == "all":
                        pass
                    else:
                        try:
                            if sample == "xshort":
                                list_report = list_report[0:8]
                            elif sample == "short":
                                list_report = list_report[0:16]
                            elif sample == "long":
                                list_report = list_report[0:32]
                            elif sample == "xlong":
                                list_report = list_report[0:48]
                        except IndexError:
                            print(
                                "SampleError: The subsampling introduced is not compatible with the simulation carried out.\
                             Please check that you have chosen a smaller sampling than the simulation."
                            )

                    for report in list_report:

                        path_report = os.path.join(path_epoch, report)
                        with open(path_report, "r") as filein:
                            cont = 0
                            for line in filein:

                                cont += 1

                                if cont != 1:
                                    sline = line.split()

                                    te.append(float(sline[3]))
                                    be.append(float(sline[4]))
                                    sasa.append(float(sline[5]))
                                    rmsd.append(float(sline[6]))

                # Normalized Total Energy
                min_energy = min(te)
                te_n = np.array(te) - min_energy

                # Binding Energy
                be_min = minimum(be)
                be_bz = boltzmann(be, te_n, 298.0)
                be_p5 = p5(be)
                be_p10 = p10(be)
                be_p25 = p25(be)

                be_list = [be_min, be_bz, be_p5, be_p10, be_p25]

                # Total energy
                te_min = minimum(te)
                te_bz = boltzmann(te, te_n, 298.0)
                te_p5 = p5(te)
                te_p10 = p10(te)
                te_p25 = p25(te)

                te_list = [te_min, te_bz, te_p5, te_p10, te_p25]

                # SASA
                sasa_min = minimum(sasa)
                sasa_bz = boltzmann(sasa, te_n, 298.0)
                sasa_av = average(sasa)
                sasa_max = maximum(sasa)

                sasa_list = [sasa_min, sasa_bz, sasa_av, sasa_max]

                # RMSD
                rmsd_max = maximum(rmsd)
                rmsd_av = average(rmsd)

                rmsd_list = [rmsd_max, rmsd_av]

            else:

                be_list = 5 * [np.NaN]
                te_list = 5 * [np.NaN]
                sasa_list = 4 * [np.NaN]
                rmsd_list = 2 * [np.NaN]

            return be_list, te_list, sasa_list, rmsd_list

        def _resultsToDictionary(
            docking_param,
            forcefield_param,
            truncated_param,
            perturbation_param,
            rescorings_param,
            system,
            be_list,
            te_list,
            sasa_list,
            rmsd_list,
            sample,
        ):
            """
            Joins all the data into a single dataframe with all the information.

            Parameters
            ==========
            docking_param : str
                Docking tool used to generate initial positions for the simulation
            forcefield_param : str
                Force field used in the simulation.
            truncated_param : str
                Protein part used in the simulation.
            perturbation_param : str
                Perturbation protocol used in the simulation.
            rescorings_param : str
                Sampling method used in the simulation.
            system : str
                Name of the folder we are analyzing (name of the ligand).
            be_list : list
                List with the different scores calculated for the binding energy.
            te_list : list
                List with the different scores calculated for the total energy.
            sasa_list : list
                List with the different scores calculated for the SASA.
            rmsd_list : list
                List with the different scores calculated for the RMSD.

            Returns
            ==========
            simulation_data_dict : dict
                Dictionary with all the data of all the and all the scores calculated.
            """

            be_min, be_bz, be_p5, be_p10, be_p25 = be_list
            te_min, te_bz, te_p5, te_p10, te_p25 = te_list
            sasa_min, sasa_bz, sasa_av, sasa_max = sasa_list
            rmsd_max, rmsd_av = rmsd_list

            if sample is not "all":
                rescorings_param = sample

            simulation_data_dict = {
                "docking_tool": docking_param,
                "forcefield": forcefield_param,
                "protein_part": truncated_param,
                "perturbation": perturbation_param,
                "sampling": rescorings_param,
                "system": system,
                "be_min": be_min,
                "be_bz": be_bz,
                "be_p5": be_p5,
                "be_p10": be_p10,
                "be_p25": be_p25,
                "te_min": te_min,
                "te_bz": te_bz,
                "te_p5": te_p5,
                "te_p10": te_p10,
                "te_p25": te_p25,
                "sasa_min": sasa_min,
                "sasa_bz": sasa_bz,
                "sasa_av": sasa_av,
                "sasa_max": sasa_max,
                "rmsd_max": rmsd_max,
                "rmsd_av": rmsd_av,
            }

            return simulation_data_dict

        (
            docking_param,
            forcefield_param,
            truncated_param,
            perturbation_param,
            rescorings_param,
        ) = _PELEOptionsRetriever(dataset_location)
        be_list, te_list, sasa_list, rmsd_list = _energyInSimulation(
            path_system, sample
        )
        simulation_data_dict = _resultsToDictionary(
            docking_param,
            forcefield_param,
            truncated_param,
            perturbation_param,
            rescorings_param,
            system,
            be_list,
            te_list,
            sasa_list,
            rmsd_list,
            sample,
        )

        return simulation_data_dict

    def experimentalDataCollector(self, path_experimental_data=None):
        """
        Retrieves the experimental energetic data to be able to plot correlations.

        Parameters
        ==========
        path_experimental_data : str
            Path to the experimental data in case it has not been inputted previously
            in the pipeline.
        """

        path_input_data = "1_input_files/experimental_energies"

        if path_experimental_data is not None:
            if os.path.isfile(path_experimental_data):
                df = pd.read_csv(path_experimental_data)

                if not os.path.isdir(path_input_data):
                    os.mkdir(path_input_data)

                shutil.move(path_experimental_data, path_input_data)
            else:
                print("is not file")
                print(" - Path inputed is not valid. Please check the path.")
                print(" - Checking if experimental data has been inputted previously.")

        elif path_experimental_data is None:
            if os.path.isdir(path_input_data):
                experimental_data_file = [
                    x for x in os.listdir(path_input_data) if x.endswith(".csv")
                ][0]
                path_input_experimental_data = os.path.join(
                    path_input_data, experimental_data_file
                )
                df = pd.read_csv(path_input_experimental_data, index_col=0)
                print(
                    " -     Experimental data found at {}".format(
                        path_input_experimental_data
                    )
                )
            else:
                raise Exception(
                    "MissingExperimentalData: The inputted path is not valid and experimental data has not been inputed previously. Please use the method .experimentalDataCollector to input the data."
                )

        self.experimental_data = df

    def equibindDataTrimming(self, df):
        """
        Selects the best scores for all the simulations performed with all the
        tautomers/stereoisomers from ligPrep. Selects simulations with best minimum
        Binding Energy.

        Parameters
        ==========
        df : pandas.DataFrame
            Dataframe with the information of all the simulations.
        """

        def equibindBestPoses(df):
            """
            Selects the best scores for all the simulations performed with all the
            tautomers/stereoisomers from ligPrep.

            Parameters
            ==========
            df : pandas.DataFrame
                Dataframe with the information of all the simulations.
            """

            df.sort_values(by="be_min", inplace=True)
            df.drop_duplicates(subset="ligand", keep="first", inplace=True)
            df.sort_values(by="ligand", inplace=True)
            df.reset_index(drop=True, inplace=True)

            print(" - Best poses stored in 5_pele_analysis/Equibind_best_poses.csv")
            df.to_csv("5_pele_analysis/Equibind_best_poses.csv", index=False)

        docking_tool = df["docking_tool"].iloc[0]

        if docking_tool == "equibind":

            # Sorting data
            df[["ligand", "conformer"]] = (
                df["system"].str.split("-", expand=True).astype(int)
            )

            groups = df.groupby(
                [
                    "docking_tool",
                    "forcefield",
                    "protein_part",
                    "perturbation",
                    "sampling",
                ]
            )

            subset_list = [group for _, group in groups if len(group) > 1]
            dfs = []
            nan_lost = 0

            for subset in subset_list:
                length_df = len(subset)
                df = subset.dropna()
                length_df_wo_NaN = len(df)

                nan_lost += length_df - length_df_wo_NaN

                df_sorted = df.sort_values(by="ligand")
                df = df_sorted.drop_duplicates("ligand")
                final_df = df.sort_values("ligand")
                dfs.append(final_df)

            print(" - {} simulations have failed.".format(nan_lost))

            combined_df = pd.concat(dfs, ignore_index=True)
            self.all_data = combined_df.copy()

            equibindBestPoses(combined_df)

        else:
            raise Exception(
                "DockingToolError: The docking tool used is {} not equibind".format(
                    docking_tool
                )
            )

        print(
            " - Information of all the simulation stored in 5_pele_analysis/Equibind_all_dataset.csv"
        )

        combined_df.to_csv("5_pele_analysis/Equibind_all_dataset.csv", index=False)
        self.equibind_data = combined_df

    def PELEDataCollector(self, sample="all"):
        """
        Generates all the folder hierarchy and calculates all the scores for all the
        metrics in the PELE simulations. Then joins all this information into a dataframe.
        Lastly, merges the experimental data into this dataframe.
        """

        self._PELEDownloadsAssemble()
        folders_to_check = self._PELEFoldersToAnalyze()
        self.experimentalDataCollector()

        all_data_dict = {}
        cont = 0

        for dataset in folders_to_check:
            dataset_location = "/".join(dataset.split("/")[2:])

            for system in [
                x
                for x in os.listdir(dataset)
                if os.path.isdir(os.path.join(dataset, x))
            ]:
                cont += 1
                path_system = os.path.join(dataset, system)

                if "-" in system:
                    simulation_data_dict = self._energyCalculator(
                        dataset_location, system, path_system, sample
                    )
                    all_data_dict[cont] = simulation_data_dict
                else:
                    simulation_data_dict = self._energyCalculator(
                        dataset_location, int(system), path_system, sample
                    )
                    all_data_dict[cont] = simulation_data_dict

        if sample is not "all":
            print(
                " -     Subsampling of the simulation successful: from {sampling} to {subsampling}.".format(
                    sampling=self.sampling, subsampling=sample
                )
            )

        # Dataframe managing
        df_experimental = self.experimental_data
        df_all_data = (
            pd.DataFrame.from_dict(all_data_dict, orient="index")
            .sort_values(by="system")
            .reset_index(drop=True)
        )
        df_calculated = df_all_data.copy()

        docking_tool = df_all_data["docking_tool"].iloc[0]

        # Distinguishing docking tools
        if docking_tool == "equibind":

            # Trimming equibind data
            self.equibindDataTrimming(df_all_data)
            df = self.all_data

            for index, row in df.iterrows():
                system_value = row["ligand"]

                if system_value in df_experimental.index:
                    dg_value = df_experimental.loc[system_value, "dG"]
                    df.loc[index, "dG"] = dg_value

            self.all_data = df.copy()

        else:
            for index, row in df_all_data.iterrows():
                system_value = row["system"]

                if system_value in df_experimental.index:
                    dg_value = df_experimental.loc[system_value, "dG"]
                    df_all_data.loc[index, "dG"] = dg_value

            self.all_data = df_all_data

        self.calculated_data = df_calculated

    def correlationPlotter(
        self, x_label, y_label, sampling, x_range=None, y_range=None, df=None
    ):
        """
        Makes plots with the x label and y label inputed by the user. By
        default uses the dataframe with all the data from all the simulations.

        Parameters
        ==========
        x_label : str
            Label of the dataframe's column you want to take as the x axis.
        y_label : str
            Label of the dataframe's column you want to take as the x axis.
        sampling : str
            Sampling used for these results.
        x_range : list
            List with the range of the x axis of the plot
        y_range : list
            List with the range of the y axis of the plot
        df : pandas.DataFrame
            Dataframe with the information the user wants to analyze.
        """

        def plotter(df, x_label, y_label, sampling):
            """
            Generates two plots: one with the original data and one with z-score.

            Parameters
            ==========
            x_label : str
                Label of the dataframe's column you want to take as the x axis.
            y_label : str
                Label of the dataframe's column you want to take as the x axis.
            sampling : str
                Sampling used for these results.
            """

            x = df[x_label].to_numpy()
            y = df[y_label].to_numpy()

            # Calculate z-scores for x and y
            z_x = (x - np.mean(x)) / np.std(x)
            z_y = (y - np.mean(y)) / np.std(y)

            m_z, n_z, r_z, p_z, _ = linregress(z_x, z_y)

            plt.figure()
            plt.scatter(z_x, z_y)

            # Set labels and title
            plt.xlabel("Z score experimental")
            plt.ylabel("Z score calculated")
            plt.title(
                "{x} vs. {y}: Z-score correlation {res_m}".format(
                    x=x_label, y=y_label, res_m=sampling
                )
            )
            plt.plot(
                z_x,
                m_z * np.array(z_x) + n_z,
                color="orange",
                label="r = {:.2f}\np = {:.2f}\nn = {}".format(r_z, p_z, len(x)),
            )
            plt.legend(loc="best")
            plt.savefig(
                "5_pele_analysis/images/{res_m}_{x}_{y}_zscore_correlation.png".format(
                    res_m=sampling, x=x_label, y=y_label
                ),
                format="png",
            )

            m, n, r, p, _ = linregress(x, y)

            plt.figure()
            plt.scatter(x, y)

            # Set labels and title
            plt.xlabel("Experimental")
            plt.ylabel("Calculated")
            plt.title(
                "{x} vs. {y}: correlation {res_m}".format(
                    x=x_label, y=y_label, res_m=sampling
                )
            )
            plt.plot(
                x,
                m * np.array(x) + n,
                color="orange",
                label="r = {:.2f}\np = {:.2f}\nn = {}".format(r, p, len(x)),
            )
            plt.legend(loc="best")
            plt.savefig(
                "5_pele_analysis/images/{res_m}_{x}_{y}_correlation.png".format(
                    res_m=sampling, x=x_label, y=y_label
                ),
                format="png",
            )

        def data_selection(df, x_label, y_label, x_range, y_range):
            """
            Trims the dataframe to make the data fit certain ranges.

            Parameters
            ==========
            df : pandas.DataFrame
                Dataframe with the information the user wants to analyze.
            x_label : str
                Label of the dataframe's column you want to take as the x axis.
            y_label : str
                Label of the dataframe's column you want to take as the x axis.
            x_range : str
                Range where there can be data in the x-axis
            y_range : str
                Range where there can be data in the y-axis
            """

            if x_range is not None:

                length_before_x = len(df)
                df = df[(df[x_label] >= x_range[0]) & (df[x_label] < x_range[1])]
                length_after_x = len(df)

                print(
                    " - {} data points deleted after x-axis trimming.".format(
                        length_before_x - length_after_x
                    )
                )

            if y_range is not None:

                length_before_y = len(df)
                df = df[(df[y_label] >= y_range[0]) & (df[y_label] < y_range[1])]
                length_after_y = len(df)
                print(
                    " - {} data points deleted after y-axis trimming.".format(
                        length_before_y - length_after_y
                    )
                )

            return df

        if df is None:

            # Dropping NaN rows.
            df_original = self.all_data
            length_df = len(df_original)
            df = df_original.dropna()
            length_df_wo_NaN = len(df)
            NaN_rows = length_df - length_df_wo_NaN

            if NaN_rows != 0:
                print(
                    " - Warning: {} rows with NaN have been deleted.".format(NaN_rows)
                )

            if (x_range is not None) or (y_range is not None):
                df = data_selection(df, x_label, y_label, x_range, y_range)

            # Creating a folder to store plots
            if not os.path.isdir("5_pele_analysis/images"):
                os.mkdir("5_pele_analysis/images")

            plotter(df, x_label, y_label, sampling)

        else:

            # Dropping NaN rows.
            df_original = df
            length_df = len(df_original)
            df = df_original.dropna()
            length_df_wo_NaN = len(df)
            NaN_rows = length_df - length_df_wo_NaN

            if NaN_rows != 0:
                print(
                    " - Warning: {} rows with NaN have been deleted.".format(NaN_rows)
                )

            if (x_range is not None) or (y_range is not None):
                df = data_selection(df, x_label, y_label, x_range, y_range)

            # Creating a folder to store plots
            if not os.path.isdir("5_pele_analysis/images"):
                os.mkdir("5_pele_analysis/images")

            plotter(df, x_label, y_label, sampling)

    def simulationAnalyzer(self, protocol, system):
        """
        Makes plots with of the specific simulation we are interested in.

        Parameters
        ==========
        protocol : list
            List with all the protocol from which we want to analyze the
            PELE simulation (i.e. [truncated, xshort]). The order of the
            options should be: [forcefield, truncated, rescoring, xshort]
        system : str
            Name of the system protein ligand we want to analyze.
        """

        def directoryManager(protocol, system):
            """
            Creates necessary directories for storing the analysis results.

            Parameters
            ==========
            protocol : list
                List with all the protocol from which we want to analyze the
                PELE simulation (i.e. [truncated, xshort]). The order of the
                options should be: [forcefield, truncated, rescoring, xshort]
            system : str
                Name of the system protein ligand we want to analyze.

            Returns
            =======
            path_system : str
                The path to the directory where the analysis results will be stored.
            """

            path_pele_analysis = "5_pele_analysis"
            path_analysis = os.path.join(path_pele_analysis, "analysis")
            path_system = os.path.join(path_analysis, "/".join(protocol), system)

            # Generating folders
            if not os.path.isdir(path_pele_analysis):
                os.mkdir(path_pele_analysis)

            if not os.path.isdir(path_analysis):
                os.mkdir(path_analysis)

            if not os.path.isdir(path_system):
                os.makedirs(path_system)

            return path_system

        def simulationChecker(protocol, system):
            """
            Checks the existence of the simulation results directory.

            Parameters
            ==========
            protocol : list
                List with all the protocol from which we want to analyze the
                PELE simulation (i.e. [truncated, xshort]). The order of the
                options should be: [forcefield, truncated, rescoring, xshort]
            system : str
                Name of the system protein ligand we want to analyze.

            Returns
            =======
            path_simulation : str
                The path to the directory containing the simulation results for the specified system and protocol.
            """

            docking_tools = ["glide", "rdock", "equibind"]
            path = "5_pele_analysis/simulations"
            docking_method = [x for x in os.listdir(path) if x in docking_tools][0]

            path_simulation = os.path.join(
                path, docking_method, "/".join(protocol), system
            )

            if os.path.isdir(path_simulation):
                print(" - Analyzing the simulation in \n\t{}".format(path_simulation))
            else:
                raise Exception(
                    "WrongPathError: The path {} does not exist. Please check how\
                                the variable protocol has to be inputed.".format(
                        path_simulation
                    )
                )

            return path_simulation

        def datasetRetriever(path_system):
            """
            Retrieves data from PELE simulation reports.

            Parameters
            ==========
            path_system : str
                The path to the directory containing the PELE simulation reports.

            Returns
            =======
            te : list
                List of Total Energy values extracted from the simulation reports.
            be : list
                List of Binding Energy values extracted from the simulation reports.
            sasa : list
                List of SASA (Solvent Accessible Surface Area) values extracted from the simulation reports.
            rmsd : list
                List of RMSD (Root Mean Square Deviation) values extracted from the simulation reports.
            step : list
                List of integers with the step value of each MC step.
            """

            te = []
            be = []
            sasa = []
            rmsd = []
            step = []

            step_collector = {}

            list_epochs_og = [
                x
                for x in os.listdir(path_system)
                if os.path.isdir(os.path.join(path_system, x))
            ]

            list_epochs = [str(x) for x in sorted([int(x) for x in list_epochs_og])]

            if len(list_epochs) != 0:
                for epoch in list_epochs:
                    path_epoch = os.path.join(path_system, epoch)
                    list_reports_og = [
                        x for x in os.listdir(path_epoch) if x.startswith("report_")
                    ]
                    list_report = sorted(
                        list_reports_og, key=lambda x: int(x.split("_")[1])
                    )

                    for report in list_report:
                        if epoch == "0":
                            cont_steps = 1
                        else:
                            cont_steps = (
                                step_collector[
                                    os.path.join(str(int(epoch) - 1), report)
                                ]
                                + 1
                            )

                        path_report = os.path.join(path_epoch, report)
                        with open(path_report, "r") as filein:
                            cont = 0
                            for line in filein:

                                cont += 1

                                if cont != 1:
                                    sline = line.split()

                                    te.append(float(sline[3]))
                                    be.append(float(sline[4]))
                                    sasa.append(float(sline[5]))
                                    rmsd.append(float(sline[6]))
                                    step.append(cont_steps)

                                    step_collector[
                                        os.path.join(str(int(epoch)), report)
                                    ] = cont_steps
                                    cont_steps += 1

            return te, be, sasa, rmsd, step

        def plotter(path_store, protocol, system, te, be, sasa, rmsd, step):
            """
            Generates and stores scatter plots for the specific simulation data.

            Parameters
            ==========
            path_store : str
                The path to the directory where the generated plots will be stored.
            protocol : list
                List with all the protocol from which we want to analyze the
                PELE simulation (i.e. [truncated, xshort]). The order of the
                options should be: [forcefield, truncated, rescoring, xshort]
            system : str
                Name of the system protein ligand we want to analyze.
            te : list
                List of Total Energy values.
            be : list
                List of Binding Energy values.
            sasa : list
                List of SASA (Solvent Accessible Surface Area) values.
            rmsd : list
                List of RMSD (Root Mean Square Deviation) values.
            step : list
                List of step values for continuity between files.
            """

            def create_scatter_plot(
                x, y, c=None, x_var_name="", y_var_name="", color_var_name=""
            ):
                """
                Creates and stores a scatter plot with optional colorbar.

                Parameters
                ==========
                x : list
                    List of x-axis data values.
                y : list
                    List of y-axis data values.
                c : list or None, optional
                    List of color data values. If provided, a colorbar will be added to the plot. Default is None.
                x_var_name : str, optional
                    Name of the x-axis variable used for constructing the filename. Default is an empty string.
                y_var_name : str, optional
                    Name of the y-axis variable used for constructing the filename. Default is an empty string.
                color_var_name : str, optional
                    Name of the color variable used for constructing the filename when a colorbar is present. Default is an empty string.
                """
                plt.figure()
                if c is None:
                    plt.scatter(x, y)
                    filename = f"{x_var_name}_{y_var_name}.png"
                else:
                    plt.scatter(x, y, c=c, cmap="viridis", marker="o")
                    colorbar = plt.colorbar()
                    colorbar.set_label("Step", fontsize=12)
                    filename = f"{x_var_name}_{y_var_name}_{color_var_name}.png"

                plt.ylabel(
                    "Binding Energy"
                    if y is be
                    else "SASA" if y is sasa else "Total Energy"
                )
                plt.xlabel(
                    "Total Energy"
                    if x is te
                    else (
                        "Binding Energy" if x is be else "RMSD" if x is rmsd else "Step"
                    )
                )
                plt.title(
                    "System {ligand} from {prot}: {y_var} vs {x_var}".format(
                        prot=protocol[-1],
                        ligand=system,
                        y_var=(
                            "Binding Energy"
                            if y is be
                            else "SASA" if y is sasa else "Total Energy"
                        ),
                        x_var=(
                            "Total Energy"
                            if x is te
                            else "Binding Energy" if x is be else "RMSD"
                        ),
                    )
                )
                if y is be:
                    plt.ylim(top=0)
                plt.tight_layout()
                plt.savefig(f"{path_store}/{filename}")
                plt.show()
                plt.close()

            # Individual scatter plots
            create_scatter_plot(te, be, x_var_name="te", y_var_name="be")
            create_scatter_plot(te, sasa, x_var_name="te", y_var_name="sasa")
            create_scatter_plot(be, sasa, x_var_name="be", y_var_name="sasa")
            create_scatter_plot(rmsd, te, x_var_name="rmsd", y_var_name="te")
            create_scatter_plot(rmsd, be, x_var_name="rmsd", y_var_name="be")
            create_scatter_plot(rmsd, sasa, x_var_name="rmsd", y_var_name="sasa")
            create_scatter_plot(step, be, x_var_name="step", y_var_name="be")
            create_scatter_plot(step, te, x_var_name="step", y_var_name="te")
            create_scatter_plot(
                te, be, rmsd, x_var_name="te", y_var_name="be", color_var_name="rmsd"
            )
            create_scatter_plot(
                rmsd,
                be,
                c=step,
                x_var_name="rmsd",
                y_var_name="be",
                color_var_name="step",
            )
            create_scatter_plot(
                te, be, c=step, x_var_name="te", y_var_name="be", color_var_name="step"
            )

            print(" - Storing images at {}".format(path_store))

        path_store_system = directoryManager(protocol, system)
        simulation_path = simulationChecker(protocol, system)
        te, be, sasa, rmsd, step = datasetRetriever(simulation_path)
        plotter(path_store_system, protocol, system, te, be, sasa, rmsd, step)

    def PELEExtractData(self, folder, output_folder="analysis"):
        """Extracts and organizes data from PELE report files into CSVs inside output_folder."""

        def read_report(file_path):
            """Reads a report_* file using the first line as headers, returns a DataFrame."""

            try:
                with open(file_path, "r") as f:
                    lines = f.readlines()
                    if not lines:
                        print(f"  Skipped empty file: {file_path}")
                        return pd.DataFrame()

                    headers = lines[0].strip().split()
                    data = []
                    for line in lines[1:]:
                        parts = line.strip().split()
                        if len(parts) != len(headers):
                            print(f"  Skipping malformed line: {line.strip()}")
                            continue
                        data.append(parts)

                df = pd.DataFrame(data, columns=headers)
                # Convert numeric columns where possible
                for col in df.columns:
                    try:
                        df[col] = pd.to_numeric(df[col])
                    except ValueError:
                        print(f"  Column '{col}' could not be converted to numeric")
                return df
            except Exception as e:
                print(f"  Error reading {file_path}: {e}")
                return pd.DataFrame()

        # Create output folder if it doesn't exist
        os.makedirs(output_folder, exist_ok=True)
        print(f"Output folder is: {output_folder}")

        equilibration_paths = []
        numbered_paths = []

        # Identify equilibration and numbered folders
        for root, dirs, _ in os.walk(folder):
            for d in dirs:
                full_path = os.path.join(root, d)
                if d.startswith("equilibration"):
                    equilibration_paths.append(full_path)
                elif re.match(r"^\d+$", d):
                    numbered_paths.append(full_path)

        # Process equilibration folders
        for eq_path in equilibration_paths:
            report_files = glob(os.path.join(eq_path, "report_*"))
            print(
                f"Processing equilibration path: {eq_path} with {len(report_files)} reports"
            )
            if not report_files:
                continue

            dfs = []
            print("Reading reports...")
            for report in report_files:
                df = read_report(report)
                if df.empty:
                    print(f"  Skipping empty DataFrame from {report}")
                    continue
                match = re.search(r"report_(\d+)", os.path.basename(report))
                trajectory = int(match.group(1)) if match else -1
                df["trajectory"] = trajectory
                dfs.append(df)

            if dfs:
                full_df = pd.concat(dfs, ignore_index=True)
                match = re.search(r"equilibration_(\d+)", eq_path)
                epoch_label = f"e{match.group(1)}" if match else "e0"
                full_df["epoch"] = epoch_label
                cols = ["epoch", "trajectory"] + [
                    col for col in full_df.columns if col not in ("epoch", "trajectory")
                ]
                full_df = full_df[cols]

                output_path = os.path.join(
                    output_folder,
                    f"equilibration_data_{match.group(1) if match else '0'}.csv",
                )
                full_df.to_csv(output_path, index=False)
                print(f"Wrote equilibration CSV to: {output_path}")

        # Process numbered simulation folders
        final_dfs = []

        for sim_path in numbered_paths:
            report_files = glob(os.path.join(sim_path, "report_*"))
            print(
                f"Processing simulation path: {sim_path} with {len(report_files)} reports"
            )
            if not report_files:
                continue

            dfs = []
            for report in report_files:
                df = read_report(report)
                if df.empty:
                    print(f"  Skipping empty DataFrame from {report}")
                    continue
                match = re.search(r"report_(\d+)", os.path.basename(report))
                trajectory = int(match.group(1)) if match else -1
                df["trajectory"] = trajectory
                dfs.append(df)

            if dfs:
                folder_name = os.path.basename(sim_path)
                full_df = pd.concat(dfs, ignore_index=True)
                full_df["epoch"] = folder_name
                cols = ["epoch", "trajectory"] + [
                    col for col in full_df.columns if col not in ("epoch", "trajectory")
                ]
                full_df = full_df[cols]
                final_dfs.append(full_df)

        if final_dfs:
            merged_df = pd.concat(final_dfs, ignore_index=True)
            output_path = os.path.join(output_folder, "simulation_data.csv")
            merged_df.to_csv(output_path, index=False)
            print(f"Wrote simulation CSV to: {output_path}")
        else:
            print("No simulation data was collected.")
