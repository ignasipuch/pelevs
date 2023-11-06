# -*- coding: utf-8 -*-
"""
This script is designed to optimize and calculate resp charges. It yeilds the necessary
files for peleffy to generate a DataLocal file.
"""

__author__ = "Ignasi Puch-Giner "
__maintainer__ = "Ignasi Puch-Giner"
__email__ = "ignasi.puchginer@bsc.es"

from schrodinger import structure
from schrodinger.application.jaguar.input import JaguarInput
from schrodinger.structure import StructureReader
from schrodinger.job.jobcontrol import launch_job

import argparse
import sys
import os
import shutil

# 1. Optimization -> 2. Charges calculation -> 3. Ligand preparation -> 4. Generation of folder with needed files

def parse_args(args):
    """
    Function
    ----------
    It parses the command-line arguments.
    Parameters

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

    parser.add_argument("-f", "--file", type=str, dest="input_file",
                        default=None, required=True, help="Name of the file corresponding to the isolated ligand with connectivity.")
    
    parsed_args = parser.parse_args(args)

    return parsed_args

def multiplicity(input_file):
    """
    Function
    ----------
    Calculate or retrieve the number of electrons and multiplicity.

    Parameters
    ----------
    - input_file : str
        Ligand file in pdb format.

    Returns
    ----------
    - electrons : int
        Number of electrons calculated with rdkit (without double bonds)  
    - parity_bool : bool
        Bool that indicates if the number of electrons is odd or even.
    """

    outfile = input_file.split('.')[0] + '_qm_charges.out'
    electrons = 0
    job_finished = False

    if os.path.isfile(outfile):

        print('\n -     Reading .out file.')

        error_bool = False

        with open(outfile) as filein:
            for line in filein:
                if 'ERROR' in line:
                    error_bool = True

        if error_bool:

            with open(outfile) as filein:
                for line in filein:
                    if 'Total' in line and 'number' in line and 'electrons' in line:
                        line = line.split()
                        electrons = int(line[4])

            if electrons%2 == 0:
                parity_bool = True
            else: parity_bool = False

        else:

            job_finished = True
            parity_bool = False

    elif job_finished == False and not os.path.isfile(outfile):

        from rdkit import Chem

        print("\n -   This code assumes there are no double bonds.\n     If that is not the case, please re-run this code.")

        m = Chem.MolFromPDBFile(input_file)
        m1 = Chem.AddHs(m)

        for atom in m1.GetAtoms():
            electrons += atom.GetAtomicNum()

        if electrons%2 == 0:
            parity_bool = True
        else: parity_bool = False

    return electrons, parity_bool, job_finished

def jaguar_input(input_file,
                 electrons,
                 parity_bool):
    """
    Function
    ----------
    Prepare input for the schrödinger job.

    Parameters
    ----------
    - input_file : str
        Ligand file in pdb format.
    - electrons : int
        Number of electrons calculated with rdkit (without double bonds)  
    - parity_bool : bool
        Bool that indicates if the number of electrons is odd or even.

    Returns
    ----------
    - jaguar_input_file : str
        Name of the output file the geometry optimization will yield.
    """

    print(" -     Writing .in file.")

    if parity_bool:

        print("     -     Spin multiplicity set to 1 due to even number of electrons (%a)." %electrons)

        options = { 
               "isolv" : "2",                # With Poisson-Boltzmann solvation model
               "maxitg" : "10",              # Optimization number of steps
               "basis" : "6-31g**",          # QM basis set
               "igeopt" : "1",               # Indicating an optimization
               "dftname" : "B3LYP-D3",       # Density Functional Method 
               "icfit" : "1",                # Monopoles centered at the atomic center
               "nogas" : "1"                 # Skip gas phase geometry optimization
              }

    else: 

        print("     -     Spin multiplicity set to 2 due to odd number of electrons (%a)." %electrons)

        options = { 
               "isolv" : "2",
               "maxitg" : "10",
               "basis" : "6-31g**",
               "igeopt" : "1",
               "dftname" : "B3LYP-D3",
               "icfit" : "1",
               "nogas" : "1",
               "multip" : "2"                # Multiplicity of the ligand singlet/doublet
              }

    jaguar_input_file = str(input_file.split('.')[0] + "_qm_charges.in")

    st = next(StructureReader(input_file))
    ji = JaguarInput(structure=st, genkeys=options)
    ji.saveAs(jaguar_input_file)

    return jaguar_input_file

def jaguar_job(jaguar_input):
    """
    Function
    ----------
    Launches the schrödinger job.

    Parameters
    ----------
    - jaguar_input : str
        Name of the job's input.
    """

    run_command = ["jaguar", "run", jaguar_input]
    print("\n -   Running Jaguar optimization under jobcontrol...")
    job = launch_job(run_command)
    job.wait()

def jaguar_charges(optimization_file,
                   charges_file):
    """
    Function
    ----------
    Prepare input for charge calculation and launch schrödinger job.

    Parameters
    ----------
    - optimization_file : str
        Name of the out file of the schrödinger optimization.
    - charges_file : str
        Name of the input file used in the charge calculation.  
    """

    with open(optimization_file, 'r') as filein:
        with open(charges_file, 'w') as fileout:
            for line in filein:
                if 'igeopt=1' in line:
                    fileout.write('igeopt=0\n')
                else:
                    fileout.write(line)

    run_command = ["jaguar", "run", charges_file]
    print("\n -   Running Jaguar charges under jobcontrol...")
    job = launch_job(run_command)
    job.wait()

def protein_preparation(input_file,
                        prep_file):

    print('\n -   Preparing ligand to have connectivity.\n')
    os.system('$SCHRODINGER/utilities/prepwizard -nohtreat -noepik -noprotassign -noimpref -noccd -delwater_hbond_cutoff 0 -NOJOBID ' + input_file + ' ' +  prep_file)

def jaguar_output(output_file):
    """
    Function
    ----------
    Prepare directory to parametrize the ligand with peleffy.

    Parameters
    ----------
    - output_file : str
        Name of the output file of the schrödinger optimization.
    """

    def file_generation(system,
                        pdb_file,
                        mae_file):
        """
        Function
        ----------
        Write the run file for MN4.

        Parameters
        ----------
        - system : str
            Name of the PDB the ligand is extracted from.
        - pdb_file : str
            Name of the PDB we are going to use to parametrize the ligand.
        - mae_file : str
            Name of the outputfile we are going to be using for the charges.

        """

        with open('run', 'w') as fileout:

            fileout.writelines(
                '#!/bin/bash\n'
                '#SBATCH -J ' + system + '_jag\n'
                '#SBATCH --output=jag.out\n'
                '#SBATCH --error=jag.err\n'
                '#SBATCH --ntasks=1\n'
                '#SBATCH --qos=debug\n'
                '#SBATCH --time=00-00:10:00\n'
                '\n'
                'module purge\n'
                'module load ANACONDA/2019.10\n'
                '\n'
                'eval "$(conda shell.bash hook)"\n'
                '\n'
                'conda activate /gpfs/projects/bsc72/conda_envs/platform/1.6.3\n'
                '\n'
                'python -m peleffy.main ' + pdb_file + ' -f OPLS2005 --charges_from_file ' + mae_file + ' --with_solvent --as_datalocal\n'
            )

    def file_copying(system,
                     output_file,
                     pdb_file,
                     mae_file):
        """
        Function
        ----------
        Copies all the important files into input anfd output.

        Parameters
        ----------
        - system : str
            Name of the PDB the ligand is extracted from.
        - output_file : str 
            Name of the output file of the schrödinger optimization.
        - pdb_file : str
            Name of the PDB we are going to use to parametrize the ligand.
        - mae_file : str
            Name of the outputfile we are going to be using for the charges.

        Returns
        -----------
        
        - path_results : str
            Path where alla the results are kept
        """

        system = output_file.split('_qm_charges_POP.01.mae')[0]

        path_results = system + '_jag'
        path_parametrization = os.path.join(path_results,system + '_param')

        if os.path.exists(path_results) == False:
            os.mkdir(path_results)
        if os.path.exists(path_parametrization) == False:
            os.mkdir(path_parametrization)

        shutil.copyfile('run', os.path.join(path_parametrization,'run'))    
        shutil.copyfile(pdb_file, os.path.join(path_parametrization,pdb_file))
        shutil.copyfile(mae_file, os.path.join(path_parametrization,mae_file))

        path_input = os.path.join(path_results,'input')
        optimization_input = system + '_qm_charges.in'
        charges_input = system + '_qm_charges_POP.in'
        original_pdb = system + '.pdb'
        mae_optimization = output_file.split('_POP.01.mae')[0] + '.mae'
        mae_charges = output_file.split('_POP.01.mae')[0] + '.01.mae'

        if os.path.exists(path_input) == False:
            os.mkdir(path_input)

        shutil.copyfile(optimization_input, os.path.join(path_input,optimization_input))    
        shutil.copyfile(charges_input, os.path.join(path_input,charges_input))
        shutil.copyfile(original_pdb, os.path.join(path_input,original_pdb))
        shutil.copyfile(mae_charges, os.path.join(path_input,mae_charges))
        shutil.copyfile(mae_optimization, os.path.join(path_input,mae_optimization))

    def deletion(system):
        """
        Function
        ----------

        Deletes all files in the working directory.
        """

        for file in os.listdir():
            if os.path.isfile(file):
                os.remove(file)
            elif os.path.isdir('{file}-001'.format(file=system)):
                shutil.rmtree('{file}-001'.format(file=system))
            else:
                continue

    system = output_file.split('_qm_charges_POP.01.mae')[0]
    path_results = system + '_jag'
    pdb_file = system + '_prep.pdb'
    mae_file = output_file

    print('\n -     Preparing job in ' + str(path_results) + '.')
    print(' -     Cleaning directory.\n')

    file_generation(system,pdb_file,mae_file)
    file_copying(system,output_file,pdb_file,mae_file)
    deletion(system)

def assemble(input_file):
    """
    Function
    ----------
    Function that assembles all the other functions to easy-to-tead format.

    Parameters
    ----------
    - input_file : str
        Ligand file in pdb format.
    """

    prep_file = input_file.split('.pdb')[0] + '_prep.pdb'
    optimization_file = input_file.split('.')[0] + '_qm_charges.01.in'
    charges_file = input_file.split('.')[0] + '_qm_charges_POP.in'
    output_file = input_file.split('.')[0] + '_qm_charges_POP.01.mae'

    if not os.path.isfile(optimization_file):

        electrons, parity_bool, job_finished = multiplicity(input_file)

        if not job_finished:

            jaguar_input_file = jaguar_input(input_file,
                                         electrons,
                                         parity_bool)
            jaguar_job(jaguar_input_file)

    if not os.path.isfile(output_file):

        jaguar_charges(optimization_file,
                       charges_file)

    if not os.path.isfile(prep_file):

        protein_preparation(input_file,prep_file)

    jaguar_output(output_file)

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

    assemble(args.input_file)

if __name__ == '__main__':

    args = parse_args(sys.argv[1:])
    main(args)

