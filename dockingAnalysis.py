import os
import shutil
import matplotlib.pyplot as plt
from scipy.stats import linregress
import numpy as np
import pandas as pd
import seaborn as sns
import csv

class DockingAnalyzer:

    def __init__(self):
        
        self.receptor = os.listdir('1_input_files/receptor')[0]
        self.ligands = os.listdir('1_input_files/ligands')[0]
        self.docking_tool = None
        self.experimental_data = None
        self.calculated_data = None

    def _correlationPlotter(self, x, y, docking_method):

        # Creating a folder to store plots
        if not os.path.isdir('3_docking_job/images'):
            os.mkdir('3_docking_job/images')

        # Calculate z-scores for x and y
        z_x = (x - np.mean(x)) / np.std(x)
        z_y = (y - np.mean(y)) / np.std(y)

        m, n, r, p, _ = linregress(z_x,z_y)
        plt.scatter(z_x,z_y)

        # Set labels and title
        plt.xlabel('Z score experimental')
        plt.ylabel('Z score calculated')
        plt.title('{} Z-score correlation'.format(docking_method))
        plt.plot(z_x,m*np.array(z_x) + n, color='orange',label='r = {:.2f}\np = {:.2f}\nn = {}'.format(r, p, len(x)))
        plt.legend(loc='best')
        plt.savefig('3_docking_job/images/{}_correlation.png'.format(docking_method), format='png')

    def _glideDockingResultsChecker(self):
        
        path_docking = '3_docking_job/job'
        path_results = os.path.join(path_docking,[x for x in os.listdir(path_docking) if x.endswith('.csv')][0])

        if not os.path.isfile(path_results):
            raise Exception('ResultsMissingError: Before initializing the object the results must be downloaded and at {}'.format(path_docking))
        
        print(' - Glide docking results found')
        
    def _glideDataFrameRetriever(self):

        path_docking = '3_docking_job/job'
        path_results = os.path.join(path_docking,[x for x in os.listdir(path_docking) if x.endswith('.csv')][0])

        # Keeping important columns
        df_og = pd.read_csv(path_results)
        columns_to_keep = ['SMILES','title','i_i_glide_lignum','r_glide_cpu_time','r_i_glide_gscore']
        df = df_og[columns_to_keep]
        
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

        # Sorting by energies and keeping only one per molecule
        df.sort_values(by="r_i_glide_gscore", inplace=True)    
        df.drop_duplicates(subset="title", keep="first", inplace=True)
        df.sort_values(by="title", inplace=True)
        df.reset_index(drop=True, inplace=True)

        print(' - Csv information imported and sorted (self.calculated_data)')

        self.calculated_data = df

    def _glideCorrelation(self, experimental_data, column_name):
    
        # Move experimental data to input data
        if not os.path.isdir('1_input_files/experimental_energies'):
            os.mkdir('1_input_files/experimental_energies')
            
        file_name = experimental_data.split('/')[-1]
        shutil.move(file_name, '1_input_files/experimental_energies/')
        
        df_experimental = pd.read_csv(os.path.join('1_input_files/experimental_energies/', file_name), index_col=0)
        df_calculated = self.calculated_data

        x = df_experimental[column_name].to_numpy()
        y = df_calculated.r_i_glide_gscore.to_numpy() 

        self._correlationPlotter(x,y,'glide')

        print(' - Correlation image generated succesfully')
        print(' - Image stored at 3_docking_job/images')

    def _glideTimePlotter(self):

        # Creating a folder to store plots
        if not os.path.isdir('3_docking_job/images'):
            os.mkdir('3_docking_job/images')
        
        df = self.calculated_data

        plt.hist(df['r_glide_cpu_time'], bins=10, alpha=0.3, color='blue', density=True)

        kde = sns.kdeplot(df['r_glide_cpu_time'], color='red')

        x_max = kde.get_lines()[0].get_data()[0][kde.get_lines()[0].get_data()[1].argmax()]

        plt.axvline(x_max, color='black', linestyle='--', label='Max KDE: {:.2f}'.format(x_max))

        plt.xlabel("Glide's cpu time (s)")
        plt.ylabel("Density")
        plt.title("Glide's time distribution")
        plt.xlim(0, df['r_glide_cpu_time'].max())
        plt.legend()
        plt.savefig('3_docking_job/images/glide_time_distribution.png', format='png')

        print(' - Time distribution figure plotted correctly.')
        print()

    def _rdockDockingResultsChecker(self):

        path_docking = '3_docking_job/job/results'
        path_results = [x for x in os.listdir(path_docking) if x.endswith('.sd')][0]

        if len(path_results) == 0 and path_results[0].endswith('sd'):
            raise Exception('ResultsMissingError: Before initializing the object the results must be downloaded and located at {}'.format(path_docking))
        
    def _rdockDataFrameGenerator(self):

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
                        data.append([filename, counter, ligand, conformer, score])
                    if '$$$$' in line:
                        counter += 1
                    if '>  <SCORE>' in line:
                        score_bool = True
                    else: score_bool = False
                    if '>  <s_lp_Variant>' in line:
                        conformer_bool = True
                    else: conformer_bool = False
                    
                    
        # Write the extracted data to a CSV file
        output_file = 'rDock_data.csv'
        with open(os.path.join(storage_path,output_file), 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['file_name', 'index', 'ligand', 'conformer', 'rdock_score'])
            writer.writerows(data)

        print(' - rDock data extraction completed.') 
        print(' - Data saved in {}'.format(os.path.join(storage_path,output_file)))   

    def _rdockDataFrameTrimmer(self):

        df = pd.read_csv('3_docking_job/rDock_data.csv')
        sorted_df = df.sort_values('rdock_score')
        unique_df = sorted_df.drop_duplicates('ligand')
        final_df = unique_df.sort_values('ligand')

        self.calculated_data = final_df

    def _rdockCorrelation(self, experimental_data, column_name):

        # Move experimental data to input data
        if not os.path.isdir('1_input_files/experimental_energies'):
            os.mkdir('1_input_files/experimental_energies')
            
        file_name = experimental_data.split('/')[-1]
        shutil.move(file_name, '1_input_files/experimental_energies/')
        
        df_experimental = pd.read_csv(os.path.join('1_input_files/experimental_energies/', file_name), index_col=0)
        df_calculated = self.calculated_data

        x = df_experimental[column_name].to_numpy()
        y = df_calculated.rdock_score.to_numpy() 

        self._correlationPlotter(x,y,'rdock')

        print(' - Correlation image generated succesfully')
        print(' - Image stored at 3_docking_job/images\n')

    def glideAnalysis(self,experimental_data, column_name):
        
        self._glideDockingResultsChecker()
        self._glideDataFrameRetriever()
        self._glideCorrelation(experimental_data, column_name)
        self._glideTimePlotter()

    def rdockAnalysis(self, experimental_data, column_name):

        self._rdockDockingResultsChecker()
        self._rdockDataFrameGenerator()
        self._rdockDataFrameTrimmer()
        self._rdockCorrelation(experimental_data, column_name)

