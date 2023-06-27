import os
import pandas as pd
import shutil
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import seaborn as sns

class PELEAnalyzer:

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

    def _PELEDownloadsAssemble(self):

        # Generating paths
        docking_tools = ['glide', 'rdock', 'equibind']
        path_pele_analysis = '5_pele_analysis'
        path_simulations = os.path.join(path_pele_analysis, 'simulations')
        path_analysis = os.path.join(path_pele_analysis, 'analysis')

        # Generating folders
        if not os.path.isdir(path_pele_analysis):
            os.mkdir(path_pele_analysis)

        if not os.path.isdir(path_analysis):
            os.mkdir(path_analysis)

        if not os.path.isdir(path_simulations):
            os.mkdir(path_simulations)

        methods = []

        # Moving the data
        for folder in [x for x in os.listdir('.') if (os.path.isdir(x) and x in docking_tools)]:
            methods.append(folder)
            shutil.move(folder, path_simulations)

        self.docking_tool = methods
        self._path_analysis = path_analysis
        self._path_simulation = path_simulations

    def _PELEFoldersToAnalyze(self):

        path_simulation = self._path_simulation
        path_analysis = self._path_analysis

        folders_to_check = []
        rescorings = ['xshort','short','long','xlong']
        
        for root, _, _ in os.walk(path_simulation):

            if os.path.basename(root) in rescorings:
                new_path = root.replace(path_simulation, path_analysis)
                folders_to_check.append(root)
                
                if not os.path.isdir(new_path):
                    os.makedirs(new_path)

        return folders_to_check
    
    def _energyCalculator(self, dataset_location, system, path_system):

        def _PELEOptionsRetriever(dataset_location):

            # PELE options
            docking_tools = ['glide','rdock','equibind']
            forcefield = ['opls','openff']
            truncated = ['full', 'truncated']
            perturbation = ['minimization', 'if', 'refinement']
            rescorings = ['xshort','short','long','xlong']

            location = dataset_location.split('/')

            # Initializing parameters to store
            docking_param = None
            forcefield_param = 'opls'
            truncated_param = 'truncated'
            perturbation_param = 'refinement'
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

            return docking_param, forcefield_param, truncated_param, perturbation_param, rescorings_param

        def _energyInSimulation(path_system):

            def minimum(vector):
                return min(vector)
            
            def boltzmann(vector,te,T):
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
                - steps : list
                    Steps associated to poses for all the simulation.

                Returns
                ----------
                - ene_bz : float
                    Value of the boltzmann weighted energy.
                """

                R = 1.985e-3 # Gas constant in kcal/mol
                vector = np.array(vector)

                exp_bz = np.exp(-te/(R*T))
                nominator = vector.dot(exp_bz)
                denominator = np.sum(exp_bz)
                ene_bz = nominator/denominator

                return ene_bz
            
            def p5(vector):
                return np.average(np.percentile(vector,5))

            def p10(vector):
                return np.average(np.percentile(vector,10))

            def p25(vector):
                return np.average(np.percentile(vector,25))

            def maximum(vector):
                return max(vector)
        
            def average(vector):
                return np.average(vector)
            

            te = []
            be = []
            sasa = []
            rmsd = []

            list_epochs = [x for x in os.listdir(path_system) if os.path.isdir(os.path.join(path_system,x))]

            if len(list_epochs) != 0:
                for epoch in list_epochs:
                    path_epoch = os.path.join(path_system,epoch)
                    for report in [x for x in os.listdir(path_epoch) if x.startswith('report')]:
                        path_report = os.path.join(path_epoch,report)
                        with open(path_report, 'r') as filein:
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
                be_bz = boltzmann(be,te_n,298.)
                be_p5 = p5(be)
                be_p10 = p10(be)
                be_p25 = p25(be)

                be_list = [be_min,be_bz,be_p5,be_p10,be_p25]

                # Total energy
                te_min = minimum(te)
                te_bz = boltzmann(te,te_n,298.)
                te_p5 = p5(te)
                te_p10 = p10(te)
                te_p25 = p25(te)

                te_list = [te_min,te_bz,te_p5,te_p10,te_p25]

                # SASA
                sasa_min = minimum(sasa)
                sasa_bz = boltzmann(sasa,te_n,298.)
                sasa_p5 = p5(sasa)
                sasa_p10 = p10(sasa)
                sasa_p25 = p25(sasa)

                sasa_list = [sasa_min,sasa_bz,sasa_p5,sasa_p10,sasa_p25]

                # RMSD
                rmsd_max = maximum(rmsd)
                rmsd_av = boltzmann(rmsd,te_n,298.)

                rmsd_list = [rmsd_max,rmsd_av]

            else: 

                be_list = 5*[np.NaN]
                te_list = 5*[np.NaN]
                sasa_list = 5*[np.NaN]
                rmsd_list = 2*[np.NaN]

            return be_list, te_list, sasa_list, rmsd_list

        def _resultsToDictionary(docking_param, forcefield_param, truncated_param, 
                                perturbation_param, rescorings_param, system, 
                                be_list, te_list, sasa_list, rmsd_list):
            
            be_min, be_bz, be_p5, be_p10, be_p25 = be_list
            te_min, te_bz, te_p5, te_p10, te_p25 = te_list
            sasa_min, sasa_bz, sasa_p5, sasa_p10, sasa_p25 = sasa_list
            rmsd_max, rmsd_av = rmsd_list

            simulation_data_dict = {'docking_tool': docking_param,
                                    'forcefield': forcefield_param,
                                    'protein_part': truncated_param,
                                    'pertrurbation': perturbation_param,
                                    'sampling': rescorings_param,
                                    'system': system,
                                    'be_min': be_min,
                                    'be_bz': be_bz,
                                    'be_p5': be_p5,
                                    'be_p10': be_p10,
                                    'be_p25': be_p25,
                                    'te_min': te_min,
                                    'te_bz': te_bz,
                                    'te_p5': te_p5,
                                    'te_p10': te_p10,
                                    'te_p25': te_p25,
                                    'sasa_min': sasa_min,
                                    'sasa_bz': sasa_bz,
                                    'sasa_p5': sasa_p5,
                                    'sasa_p10': sasa_p10,
                                    'sasa_p25': sasa_p25,
                                    'rmsd_max': rmsd_max,
                                    'rmsd_av': rmsd_av}
            
            return simulation_data_dict
            
        docking_param, forcefield_param, truncated_param, perturbation_param, rescorings_param = _PELEOptionsRetriever(dataset_location)
        be_list, te_list, sasa_list, rmsd_list = _energyInSimulation(path_system)
        simulation_data_dict = _resultsToDictionary(docking_param, forcefield_param, truncated_param, 
                                perturbation_param, rescorings_param, system, 
                                be_list, te_list, sasa_list, rmsd_list)
    
        return simulation_data_dict

    def experimentalDataCollector(self, path_experimental_data=None):
        
        path_input_data = '1_input_files/experimental_energies'

        if os.path.isdir(path_input_data):
            experimental_data_file = [x for x in os.listdir(path_input_data) if x.endswith('.csv')][0]
            path_input_experimental_data = os.path.join(path_input_data,experimental_data_file)

        if os.path.isfile(path_input_experimental_data):
            df = pd.read_csv(path_input_experimental_data, index_col=0)
        elif (path_experimental_data is not None) and (os.path.isfile(path_experimental_data)):
            df = pd.read_csv(path_experimental_data)
        else:
            raise Exception('ExperimentalDataMissing: Experimental data should be located at {path_input} or be inputted with this method.'.format(path_input_experimental_data))
            
        self.experimental_data = df

    def PELEDataCollector(self):

        self._PELEDownloadsAssemble()
        folders_to_check = self._PELEFoldersToAnalyze()
        self.experimentalDataCollector()

        all_data_dict = {}
        cont = 0

        for dataset in folders_to_check:
            dataset_location = '/'.join(dataset.split('/')[2:])

            for system in [x for x in os.listdir(dataset) if os.path.isdir(os.path.join(dataset,x))]:
                cont += 1
                path_system = os.path.join(dataset,system)
                simulation_data_dict = self._energyCalculator(dataset_location, int(system), path_system)
                all_data_dict[cont] = simulation_data_dict

        # Dataframe managing        
        df_experimental = self.experimental_data
        df_all_data = pd.DataFrame.from_dict(all_data_dict, orient='index').sort_values(by='system').reset_index(drop=True)

        df_total = df_all_data.merge(df_experimental['dG'], left_index=True, right_index=True)

        self.calculated_data = df_all_data
        self.all_data = df_total

    def correlationPlotter(self, x_label, y_label, df=None):

        def plotter(df, x_label, y_label):

            x = df[x_label].to_numpy()
            y = df[y_label].to_numpy()

            # Calculate z-scores for x and y
            z_x = (x - np.mean(x)) / np.std(x)
            z_y = (y - np.mean(y)) / np.std(y)

            m_z, n_z, r_z, p_z, _ = linregress(z_x, z_y)

            plt.figure()
            plt.scatter(z_x, z_y)

            # Set labels and title
            plt.xlabel('Z score experimental')
            plt.ylabel('Z score calculated')
            plt.title('{x} vs. {y}: Z-score correlation'.format(x=x_label, y=y_label))
            plt.plot(z_x, m_z*np.array(z_x) + n_z, color='orange',
                    label='r = {:.2f}\np = {:.2f}\nn = {}'.format(r_z, p_z, len(x)))
            plt.legend(loc='best')
            plt.savefig(
                '5_pele_analysis/images/{x}_{y}_zscore_correlation.png'.format(x=x_label, y=y_label), format='png')

            m, n, r, p, _ = linregress(x, y)

            plt.figure()
            plt.scatter(x, y)

            # Set labels and title
            plt.xlabel('Experimental')
            plt.ylabel('Calculated')
            plt.title('{x} vs. {y}: correlation'.format(x=x_label, y=y_label))
            plt.plot(x, m*np.array(x) + n, color='orange',
                    label='r = {:.2f}\np = {:.2f}\nn = {}'.format(r, p, len(x)))
            plt.legend(loc='best')
            plt.savefig(
                '5_pele_analysis/images/{x}_{y}_correlation.png'.format(x=x_label, y=y_label), format='png')

        if df is None:

            # Dropping NaN rows.
            df_original = self.all_data
            length_df = len(df_original)
            df = df_original.dropna()
            length_df_wo_NaN = len(df)

            print(' - Warning: {} rows with NaN have been deleted.'.format(length_df-length_df_wo_NaN))

            # Creating a folder to store plots
            if not os.path.isdir('5_pele_analysis/images'):
                os.mkdir('5_pele_analysis/images')

            plotter(df, x_label, y_label)

        else: 

            # Dropping NaN rows.
            df_original = df
            length_df = len(df_original)
            df = df_original.dropna()
            length_df_wo_NaN = len(df)

            print(' - Warning: {} rows with NaN have been deleted.'.format(length_df-length_df_wo_NaN))

            # Creating a folder to store plots
            if not os.path.isdir('5_pele_analysis/images'):
                os.mkdir('5_pele_analysis/images')

            plotter(df, x_label, y_label)

        

        
