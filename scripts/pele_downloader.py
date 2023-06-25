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

    parser.add_argument("-o", "--ouput_name", type=list, dest="output_name",
                        default=None, help="Name of the folder to be created.")
    
    parsed_args = parser.parse_args(args)

    return parsed_args

def pele_reports_retriever(output_name):
    
    def folder_hierarchy(output_name):

        source_directory = 'pele_simulation'
        destination_directory = output_name

        if not os.path.isdir(output_name):
            os.mkdir(output_name)

        shutil.copytree(source_directory, destination_directory)

    def folders_with_simulations():

        folders_to_check = []

        rescoring_methods = ['xshort','short','long','xlong']
        for direct in os.walk('4_pele_simulation/'):
            if os.path.basename(direct[0]) in rescoring_methods:
                folders_to_check.append(direct[0])

        return folders_to_check

    folder_hierarchy(output_name)
    folders_to_check = folders_with_simulations()

    
def main(args):
    """
    Function
    ----------
    It reads the command-line arguments and runs lice_results.

    Parameters
    ----------
    - args : argparse.Namespace
        It contains the command-line arguments that are supplied by the user
    """

    pele_reports_retriever(output_name=args.output_name)

if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)