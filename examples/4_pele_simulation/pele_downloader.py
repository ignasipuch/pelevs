# -*- coding: utf-8 -*-
"""
This script is designed to retrieve all the important data
in a docking protocol pipeline.
"""

__author__ = "Ignasi Puch-Giner "
__maintainer__ = "Ignasi Puch-Giner"
__email__ = "ignasi.puchginer@bsc.es"

import sys
import argparse
import os
import shutil
import pandas as pd

def parse_args(args):
    """
    Function
    ----------
    It parses the command-line arguments.

    Parameters
    ----------
    - args : list[str]
        List of command-line arguments to parse

    Returns
    ----------
    - parsed_args : argparse.Namespace
        It contains the command-line arguments that are supplied by the user
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-a", "--analysis_folder_name", type=str, dest="analysis_folder_name",
                        default=None, help="Name of the folder we want to analyze.")
    parser.add_argument("-p", "--protocol_name", type=str, dest="protocol_name",
                        default=None, help="Name of the folder that contains all the simulations with all the ligands.")
    
    parser.add_argument("-o", "--ouput_name", type=str, dest="output_name",
                        default=None, help="Name of the folder to be created.")
    
    parsed_args = parser.parse_args(args)

    return parsed_args

def pele_reports_retriever(analysis_folder_name,
                           protocol_name,
                           output_name):
    """
    Function
    ----------
    Retrieves all the reports done in all the simulations belonging to single
    pele_simulation folder.

    Parameters
    ----------
    - output_name : str
        Name of the folder where the data is going to be stored.
    """   

    def folders_with_simulations(analysis_folder_name,
                                 protocol_name,
                                 output_name):
        """
        Function
        ----------
        Stores all the directory trees that have at some point a 
        rescoring method name: xshort, short, long or xlong. It stores 
        them in a variable named folders_to_check.

        Parameters
        ----------
        - output_name : str
            Name of the folder where the data is going to be stored.

        Returns
        ----------
        - folders_to_check : list
            List of paths that are going to contain simulations of the 
            whole dataset.
        """  

        if not os.path.isdir(output_name):
            os.mkdir(output_name)

        folders_to_check = []
        rescorings = ['xshort','short','long','xlong','xxlong']

        if protocol_name is not None:
            rescorings.append(protocol_name)

        if analysis_folder_name == 'standard':
            analysis_folder_name = 'pele_simulation/'
        
        for root, _, _ in os.walk('{}/'.format(analysis_folder_name)):

            if os.path.basename(root) in rescorings:
                new_path = root.replace(analysis_folder_name, output_name)
                folders_to_check.append(root)
                
                if not os.path.isdir(new_path):
                    os.makedirs(new_path)

        return folders_to_check

    def dataset_retriever(analysis_folder_name, 
                          output_name, 
                          folders_to_check):
        """
        Function
        ----------
        Searches all the directories inputed with the list folders_to_check
        and copies all the reports into the newly created folder to download.

        Parameters
        ----------
        - output_name : str
            Name of the folder where the data is going to be stored.
        - folders_to_check : list
            List of paths that are going to contain simulations of the 
            whole dataset.

        Returns
        ----------
        - failed_simulations : list
            List of simulation's paths where the simulations have failed.
        """  

        def simulation_retriever(analysis_folder_name,
                                 output_name, 
                                 simulation_path):
            """
            Function
            ----------
            Copies all the reports into the wanted directories and generates 
            the failed simulations_list.

            Parameters
            ----------
            - output_name : str
                Name of the folder where the data is going to be stored.
            - simulation_path : str
                System folder to check for simulations and with that reports.

            Returns
            ----------
            - rescoring_path : str
                Path of that is being checked.
            - system : str  
                System that is being checked.
            - failed_bool : bool
                Boolean indicating whether the simulation has failed or not.
            """  

            if analysis_folder_name is None:
                analysis_folder_name = 'pele_simulation'

            system = simulation_path.split('/')[-1]
            rescoring_path = '/'.join(simulation_path.split('/')[:-1])
            adaptive_simulation_path = '6_adaptive_pele_simulation/'
            complex_name = [x for x in os.listdir(os.path.join(simulation_path,system,adaptive_simulation_path)) if x.startswith('complex')][0]
            nbd_suite_path = '6_adaptive_pele_simulation/{}/output'.format(complex_name)
            ligand_folder = simulation_path.replace(analysis_folder_name, output_name)

            whole_path = os.path.join(simulation_path,system,nbd_suite_path)
            
            if os.path.isdir(whole_path):

                failed_bool = False

                for epoch in [x for x in os.listdir(whole_path) if x.isdigit()]:

                    epoch_path = os.path.join(whole_path,epoch) 
                    new_epoch_path = os.path.join(ligand_folder,epoch)
                    
                    if not os.path.isdir(new_epoch_path):
                        os.mkdir(new_epoch_path)

                    for report in [x for x in os.listdir(epoch_path) if x.startswith('report')]:
                        report_path = os.path.join(epoch_path,report)
                        
                        shutil.copy(report_path,new_epoch_path)

            else: failed_bool = True

            return rescoring_path, system, failed_bool

        failed_simulations = []

        if analysis_folder_name is None:
            analysis_folder_name = 'pele_simulation/' 
               
        for folder in folders_to_check:
            for simulation in os.listdir(folder):
                simulation_path = os.path.join(folder,simulation)

                if os.path.isdir(simulation_path):
                    ligand_folder = simulation_path.replace(analysis_folder_name, output_name)

                    if not os.path.isdir(ligand_folder):
                        os.mkdir(ligand_folder)
                        
                    rescoring_path, system, failed_bool = simulation_retriever(analysis_folder_name,
                                                                               output_name, 
                                                                               simulation_path) 
                    
                    if failed_bool:
                        failed_simulations.append(os.path.join(rescoring_path,system))
        
        print(' - Process complete.')

        return failed_simulations

    def failed_simulations_writer(output_name, failed_simulations):
        """
        Function
        ----------
        Searches all the directories inputed with the list folders_to_check
        and copies all the reports into the newly created folder to download.

        Parameters
        ----------
        - output_name : str
            Name of the folder where the data is going to be stored.
        - failed_simulations : list
            List of paths that where the simulation has failed.
        """  

        csv_path = os.path.join(output_name,'failed_simulations.csv')
        df = pd.DataFrame()
        df['failed_simulations'] = failed_simulations
        df.to_csv(csv_path, index=False)

        print(' - Total number of failed simulations: {}.'.format(len(failed_simulations)))
        print(' - All data about failed simulations stored in failed_simulations.csv.')

    folders_to_check = folders_with_simulations(analysis_folder_name,
                                                protocol_name,
                                                output_name)
    failed_simulations = dataset_retriever(analysis_folder_name,
                                           output_name,
                                           folders_to_check)
    failed_simulations_writer(output_name,
                              failed_simulations)
    
def main(args):
    """
    Function
    ----------
    It reads the command-line arguments and runs pele_reports_retriever.

    Parameters
    ----------
    - args : argparse.Namespace
        It contains the command-line arguments that are supplied by the user
    """

    pele_reports_retriever(analysis_folder_name=args.analysis_folder_name,
                           protocol_name=args.protocol_name,
                           output_name=args.output_name)

if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)