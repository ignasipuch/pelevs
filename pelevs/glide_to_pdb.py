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
    parser.add_argument("--multiple_poses", action='store_true', dest="multiple_poses",
                        default=False, help="Bool to know whether there are multiple poses per ligand.")
     
    parsed_args = parser.parse_args(args)

    return parsed_args

def glide_pdb(job_name, multiple_poses):

    def indeces_maegz(path_csv, path_select_data, multiple_poses):

        if not multiple_poses:

            # Sorting csv to match maegz structures indeces and sorting by energy
            df = pd.read_csv(path_csv)
            df_csv_sort = df.sort_values('r_i_docking_score').reset_index(drop=True)

            # Droping duplicates and only keeping best energies
            df_result = df_csv_sort.drop_duplicates(['title', 'i_i_glide_lignum'])
            sorted_df = df_result.sort_values(['title','i_i_glide_lignum'])

            # Saving best conformers
            sorted_df.to_csv('3_docking_job/job/{}.csv'.format('best_scores_conformers'))

            # Only keeping one conformer per ligand
            sorted_df = sorted_df.sort_values('r_i_docking_score')
            sorted_df = sorted_df.drop_duplicates('title')
            sorted_df = sorted_df.sort_values('title')

            # Saving best ligands
            sorted_df.to_csv('3_docking_job/job/{}.csv'.format('best_scores_ligands'))

            indeces = []
            molecules = []

            for i in sorted_df.title.unique():
                df_select = sorted_df[sorted_df['title'] == i]
                indeces.append(df_select.index[0] + 1)
                molecules.append(i)

        if multiple_poses:

            df2 = pd.read_csv(path_select_data)
            indeces_with_protein = df2['original_index'].to_numpy() + 1
            indeces = list(indeces_with_protein)
            molecules = df2['title'].to_list()

        return indeces, molecules

    # Setting paths
    path = str(pathlib.Path().absolute())
    path_input = os.path.join(path, '3_docking_job/job/{}'.format(job_name + '.csv'))
    path_input_select_data = os.path.join(path, '3_docking_job/Glide_dataset.csv')
    path_maegz = os.path.join(path, '3_docking_job/job/{}'.format(job_name + '_pv.maegz'))
    output_directory = '3_docking_job/job/output_pdb_files'

    # Creating output directory
    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    # Calculating indeces needed from the maegz
    indeces, molecules = indeces_maegz(path_input, path_input_select_data, multiple_poses)

    structure_reader = structure.StructureReader(path_maegz)

    cont = 0 
    for i, st in enumerate(structure_reader):
        if i in indeces:
            cont += 1 
            if multiple_poses:
                index = indeces.index(i)
                pdb_filename = os.path.join(output_directory, '{}_{}.pdb'.format(str(molecules[index]),str(cont)))
                structure_writer = structure.StructureWriter(pdb_filename)
                structure_writer.append(st)
                structure_writer.close()
            else:
                index = indeces.index(i)
                pdb_filename = os.path.join(output_directory, '{}.pdb'.format(str(molecules[index])))
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

    glide_pdb(job_name=args.job_name,
              multiple_poses=args.multiple_poses)

if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)