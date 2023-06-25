# -*- coding: utf-8 -*-
"""
This script is designed to convert glide's maegz to pdbs.
"""

__author__ = "Ignasi Puch-Giner "
__maintainer__ = "Ignasi Puch-Giner"
__email__ = "ignasi.puchginer@bsc.es"

import sys
import argparse
import pandas as pd
import os
import pathlib
from schrodinger import structure

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

    parser.add_argument("-jn", "--job_name", type=str, dest="job_name",
                        default=None, help="Name of the docking job.")
     
    parsed_args = parser.parse_args(args)

    return parsed_args

def glide_pdb(job_name):

    def indeces_maegz(path_csv):
        df = pd.read_csv(path_csv)
        df_sorted = df.sort_values('r_i_glide_gscore')
        df_result = df_sorted.drop_duplicates(['title', 'i_i_glide_lignum'])
        sorted_df = df_result.sort_values(['title','i_i_glide_lignum'])

        sorted_df.to_csv('3_docking_job/job/{}.csv'.format('best_scores_conformers'))

        sorted_df = sorted_df.sort_values('r_i_glide_gscore')
        sorted_df = sorted_df.drop_duplicates('title')
        sorted_df = sorted_df.sort_values('title')

        sorted_df.to_csv('3_docking_job/job/{}.csv'.format('best_scores_ligands'))

        indeces = []
        molecules = []

        for i in sorted_df.title.unique():
            df_select = sorted_df[sorted_df['title'] == i]
            indeces.append(df_select.index[0] + 1)
            molecules.append(i)

        return indeces, molecules

    # Setting paths
    path = str(pathlib.Path().absolute())
    path_input = os.path.join(path, '3_docking_job/job/{}'.format(job_name + '.csv'))
    path_maegz = os.path.join(path, '3_docking_job/job/{}'.format(job_name + '_pv.maegz'))
    output_directory = '3_docking_job/job/output_pdb_files/'

    # Creating output directory
    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    # Calculating indeces needed from the maegz
    indeces, molecules = indeces_maegz(path_input)

    structure_reader = structure.StructureReader(path_maegz)

    for i, st in enumerate(structure_reader):
        if i in indeces:
            index = indeces.index(i)
            pdb_filename = f'{output_directory}{str(molecules[index])}.pdb'
            structure_writer = structure.StructureWriter(pdb_filename)
            structure_writer.append(st)
            structure_writer.close()
    
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

    glide_pdb(job_name=args.job_name)

if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)