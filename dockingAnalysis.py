import os
import shutil
import matplotlib.pyplot as plt
from scipy.stats import linregress
import numpy as np
import pandas as pd
import seaborn as sns
import csv


class DockingAnalyzer:
    """
    Attributes
    ==========
    receptor : str
        Name of the receptor's file (sdf, mol2 or pdb).
    ligands : str
        Name of the file were ligands in a csv with SMILES are located.
    docking_tool : str
        Docking software wanted.
    experimental_data : str
        File containing the experimental data sorted by ligand name.
    calculated_data : str
        File containing the docking data sorted by ligand name.

    Methods
    =======
    glideAnalysis(self, experimental_data, column_name)
        Calculate energetic correlation between Glide's predictions and
        experimental data. Also it plots the distribution of time spent
        per ligand.
    rdockAnalysis(self,  experimental_data, column_name)
        Calculate energetic correlation between rDock's predictions and
        experimental data. 


    Hidden Methods
    ==============
    _correlationPlotter(self, x, y, docking_method)
        Plotts x and y in a z_score format and stores 
        the image.
    _glideDockingResultsChecker(self)
        Checks if the results from the Glide's docking
        are in place.
    _glideDataFrameRetriever(self)
        Retrieves csv generated by Glide's docking, trims,
        and adds important information.
    _glideCorrelation(self, experimental_data, column_name)
        Generates the directory to store plots and obtains
        x and y vectors to pass onto _correlationPlotter.
    _glideTimePlotter(self)
        Plots a histogram with Glide's time data.
    _rdockDockingResultsChecker(self) 
        Checks if the results from the rDock's docking
        are in place.  
    _rdockDataFrameGenerator(self)
        Generates a dataframe with all the important information
        from all the fils generated with the rDock docking.
    _rdockDataFrameTrimmer(self)
        Modifies and trims the original dataframe to obtain the 
        one that we are interested in.
    _rdockCorrelation(self, experimental_data, column_name)
        Generates the directory to store plots and obtains
        x and y vectors to pass onto _correlationPlotter.
    """

    def __init__(self):
        """
        Initialize object and assign atributes.
        """

        self.receptor = os.listdir('1_input_files/receptor')[0]
        self.ligands = os.listdir('1_input_files/ligands')[0]
        self.docking_tool = None
        self.experimental_data = None
        self.calculated_data = None

    def _correlationPlotter(self, x, y, docking_method):
        """
        Makes a scatter plot of the two vectors' z-score
        and finds the correlation between them

        Parameters
        ==========
        x : np.array
            Array with the x-axis values
        y : np.array
            Array with the y-axis values
        docking_method : str
            Method with which the values have been obtained.
        """

        # Creating a folder to store plots
        if not os.path.isdir('3_docking_job/images'):
            os.mkdir('3_docking_job/images')

        # Calculate z-scores for x and y
        z_x = (x - np.mean(x)) / np.std(x)
        z_y = (y - np.mean(y)) / np.std(y)

        m_z, n_z, r_z, p_z, _ = linregress(z_x, z_y)

        plt.figure()
        plt.scatter(z_x, z_y)

        # Set labels and title
        plt.xlabel('Z score experimental')
        plt.ylabel('Z score calculated')
        plt.title('{} Z-score correlation'.format(docking_method))
        plt.plot(z_x, m_z*np.array(z_x) + n_z, color='orange',
                 label='r = {:.2f}\np = {:.2f}\nn = {}'.format(r_z, p_z, len(x)))
        plt.legend(loc='best')
        plt.savefig(
            '3_docking_job/images/{}_zscore_correlation.png'.format(docking_method), format='png')

        m, n, r, p, _ = linregress(x, y)

        plt.figure()
        plt.scatter(x, y)

        # Set labels and title
        plt.xlabel('Experimental')
        plt.ylabel('Calculated')
        plt.title('{} correlation'.format(docking_method))
        plt.plot(x, m*np.array(x) + n, color='orange',
                 label='r = {:.2f}\np = {:.2f}\nn = {}'.format(r, p, len(x)))
        plt.legend(loc='best')
        plt.savefig(
            '3_docking_job/images/{}_correlation.png'.format(docking_method), format='png')

    def _glideDockingResultsChecker(self):
        """
        Checks if the results obtained with glide have been downloaded 
        in the correct path.
        """

        path_docking = '3_docking_job/job'
        path_results = os.path.join(
            path_docking, [x for x in os.listdir(path_docking) if x.endswith('.csv')][0])

        if not os.path.isfile(path_results):
            raise Exception(
                'ResultsMissingError: Before initializing the object the results must be downloaded and at {}'.format(path_docking))

        print(' - Glide docking results found')

    def _glideDataFrameRetriever(self):
        """
        Retrieves the data frame generated with the Glide docking.
        It also modifies certain columns and values to end up
        having a dataframe with the values that we are interested
        in such as SMILES, name of the ligand, ligand number,
        time, and score. It also stores only the best conformation
        per ligand to retrieve the original number of ligands of
        the dataset (since ligprep generates conformers).
        """

        path_docking = '3_docking_job/job'
        path_results = os.path.join(
            path_docking, [x for x in os.listdir(path_docking) if x.endswith('.csv')][0])

        # Keeping important columns
        df_og = pd.read_csv(path_results)
        columns_to_keep = ['SMILES', 'title', 'i_i_glide_lignum',
                           'r_glide_cpu_time', 'r_i_glide_gscore']
        df = df_og[columns_to_keep].copy()

        # Adding molecule number to the dataframe
        prev_title = None
        prev_value = None
        modified_i_i_glide_lignum = np.zeros(df.shape[0], dtype=int)

        for i, row in df.iterrows():
            if row['title'] != prev_title:
                prev_title = row['title']
                prev_value = row['i_i_glide_lignum']
            modified_i_i_glide_lignum[i] = row['i_i_glide_lignum'] - prev_value

        df.insert(2, 'conformation', modified_i_i_glide_lignum + 1)

        df.to_csv('3_docking_job/Glide_whole_dataset.csv')

        # Sorting by energies and keeping only one per molecule
        df.sort_values(by="r_i_glide_gscore", inplace=True)
        df.drop_duplicates(subset="title", keep="first", inplace=True)
        df.sort_values(by="title", inplace=True)
        df.reset_index(drop=True, inplace=True)

        df.to_csv('3_docking_job/Glide_dataset.csv')

        print(' - Csv information imported and sorted (self.calculated_data)')

        self.calculated_data = df

    def _glideCorrelation(self, experimental_data, column_name):
        """
        Uses _correlationPlotter to plot the calculated vs the
        experimental values of the energies involved in a 
        glide docking.

        Parameters
        ==========
        experimental_data : str
            Name of the csv file with the experimental data.
        column_name : str
            Name of the column where the data in the csv is stored.
        """

        file_name = experimental_data.split('/')[-1]

        # Move experimental data to input data
        if not os.path.isdir('1_input_files/experimental_energies'):
            os.mkdir('1_input_files/experimental_energies')
            shutil.move(file_name, '1_input_files/experimental_energies/')

        df_experimental = pd.read_csv(os.path.join(
            '1_input_files/experimental_energies/', file_name), index_col=0)
        df_calculated = self.calculated_data

        self.experimental_data = df_experimental

        x = df_experimental[column_name].to_numpy()
        y = df_calculated.r_i_glide_gscore.to_numpy()

        self._correlationPlotter(x, y, 'glide')

        print(' - Correlation image generated succesfully')
        print(' - Image stored at 3_docking_job/images\n')

    def _glideTimePlotter(self):
        """
        Makes a histogram plot to show the distribution of times 
        invested in performing the docking.
        """

        df = self.calculated_data

        plt.figure()
        plt.hist(df['r_glide_cpu_time'], bins=10,
                 alpha=0.3, color='blue', density=True)

        kde = sns.kdeplot(df['r_glide_cpu_time'], color='red')

        x_max = kde.get_lines()[0].get_data()[
            0][kde.get_lines()[0].get_data()[1].argmax()]

        plt.axvline(x_max, color='black', linestyle='--',
                    label='Max KDE: {:.2f}'.format(x_max))

        plt.xlabel("Glide's cpu time (s)")
        plt.ylabel("Density")
        plt.title("Glide's time distribution")
        plt.xlim(0, df['r_glide_cpu_time'].max())
        plt.legend()
        plt.savefig(
            '3_docking_job/images/glide_time_distribution.png', format='png')

        print(' - Time distribution figure plotted correctly.')
        print()

    def _rdockDockingResultsChecker(self):
        """
        Checks if the results obtained with rdock have been downloaded 
        in the correct path.
        """

        path_docking = '3_docking_job/job/results'
        path_results = [x for x in os.listdir(
            path_docking) if x.endswith('.sd')][0]

        if len(path_results) == 0 and path_results[0].endswith('sd'):
            raise Exception(
                'ResultsMissingError: Before initializing the object the results must be downloaded and located at {}'.format(path_docking))

    def _rdockDataFrameGenerator(self):
        """
        It generates a dataframe with all the important information
        stored in the multiple files, such as: file name, location index,
        ligand name, conformer, and rdock score.
        """

        # Folder path containing the files
        folder_path = '3_docking_job/job/results'
        storage_path = '3_docking_job/'

        data = []

        for filename in [x for x in os.listdir(folder_path) if x.startswith('split')]:
            file_path = os.path.join(folder_path, filename)

            counter = 1
            score_bool = False
            conformer_bool = False

            # Open the file
            with open(file_path, 'r') as file:
                for line in file:
                    if score_bool:
                        score = line.split()[0]
                    if conformer_bool:
                        ligand, conformer = line.split('-')
                        data.append(
                            [filename, counter, ligand, conformer, score])
                    if '$$$$' in line:
                        counter += 1
                    if '>  <SCORE>' in line:
                        score_bool = True
                    else:
                        score_bool = False
                    if '>  <s_lp_Variant>' in line:
                        conformer_bool = True
                    else:
                        conformer_bool = False

        # Write the extracted data to a CSV file
        output_file = 'rDock_data.csv'
        with open(os.path.join(storage_path, output_file), 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['file_name', 'index', 'ligand',
                             'conformer', 'rdock_score'])
            writer.writerows(data)

        print(' - rDock data extraction completed.')
        print(' - Data saved in {}'.format(os.path.join(storage_path, output_file)))

    def _rdockDataFrameTrimmer(self):
        """
        Modifies and trims the rdock dataframe to have only 
        the best score per ligand sorted by name of the ligand
        to end up having the same number of ligands than in the 
        original dataset (since ligprep generates multiple conformers).
        """

        df = pd.read_csv('3_docking_job/rDock_data.csv')
        sorted_df = df.sort_values('rdock_score')
        unique_df = sorted_df.drop_duplicates('ligand')
        final_df = unique_df.sort_values('ligand')
        final_df.to_csv('3_docking_job/rDock_best_poses.csv', index=False)

        print(' - Csv generated at 3_docking_job/rDock_best_poses.csv with best poses.')

        self.calculated_data = final_df

    def _rdockCorrelation(self, experimental_data, column_name):
        """
        Uses _correlationPlotter to plot the calculated vs the
        experimental values of the energies involved in a 
        rdock docking.

        Parameters
        ==========
        experimental_data : str
            Name of the csv file with the experimental data.
        column_name : str
            Name of the column where the data in the csv is stored.
        """

        file_name = experimental_data.split('/')[-1]

        # Move experimental data to input data
        if not os.path.isdir('1_input_files/experimental_energies'):
            os.mkdir('1_input_files/experimental_energies')
            shutil.move(file_name, '1_input_files/experimental_energies/')

        df_experimental = pd.read_csv(os.path.join(
            '1_input_files/experimental_energies/', file_name), index_col=0)
        df_calculated = self.calculated_data

        self.experimental_data = df_experimental

        x = df_experimental[column_name].to_numpy()
        y = df_calculated.rdock_score.to_numpy()

        self._correlationPlotter(x, y, 'rdock')

        print(' - Correlation image generated succesfully')
        print(' - Image stored at 3_docking_job/images\n')

    def glideAnalysis(self, experimental_data, column_name):
        """
        Uses different hidden methods to retrieve all the data 
        from the glide docking simulation and generate an 
        energy correlation plot as well as a histogram of the 
        ditribution of time spent per ligand.

        Parameters
        ==========
        experimental_data : str
            Name of the csv file with the experimental data.
        column_name : str
            Name of the column where the data in the csv is stored.
        """

        self._glideDockingResultsChecker()
        self._glideDataFrameRetriever()
        self._glideCorrelation(experimental_data, column_name)
        self._glideTimePlotter()

    def rdockAnalysis(self, experimental_data, column_name):
        """
        Uses different hidden methods to retrieve all the data 
        from the rdock docking simulation and generate an 
        energy correlation plot.

        Parameters
        ==========
        experimental_data : str
            Name of the csv file with the experimental data.
        column_name : str
            Name of the column where the data in the csv is stored.
        """

        self._rdockDockingResultsChecker()
        self._rdockDataFrameGenerator()
        self._rdockDataFrameTrimmer()
        self._rdockCorrelation(experimental_data, column_name)
