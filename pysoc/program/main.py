#!/usr/bin/env python3.6
'''Spin-Orbit Coupling(SOC) with interface to
   td-dft(Gaussian)/td-dftb(DFTB+, Thomas Niehaus)
   by Xing Gao, MPI@Muelheim, 2015
   
   Rewritten for python3 by Oliver Lee
'''

# General imports.
import sys
import argparse
from pathlib import Path
from tabulate import tabulate
import io
import csv

# PySOC imports.
import pysoc
from pysoc.io.SOC import Calculator

# Program setup.

# Printable name of this program.
NAME = "PySOC spin-orbit coupling"
DESCRIPTION = "Calculate spin-orbit coupling"
EPILOG = "{} V{}. Last updated {}.".format(NAME, pysoc.version, pysoc.last_updated.strftime("%d/%m/%Y"))

def main(
        calc_file,
        num_singlets = None,
        num_triplets = None,
        QM_program = None,
        calculation = None,
        SOC_scale = None,
        output = None,
        include_ground = None,
        CI_coefficient_threshold = None,
        print_csv = None,
        **aux_files):
    """
    Main program function for PySOC controller program.
    
    :param calc_file: The main QM output file (.log for Gaussian, .xyz for DFTB+). Other required QM output files will be found automatically based on the location of this file.
    :param QM_program: A string identifying the QM program to interface with (currently, one of either 'Gaussian' or 'DFTB+'.
    :param num_singlets: The number of singlet excited states to calculate SOC for. This should not exceed the number of singlets calculated by the QM program.
    :param num_triplets: The number of triplet excited states to calculate SOC for. This should not exceed the number of triplets calculated by the QM program.
    :param soc_scale: Scaling factor for Zeff.
    :param output: Path to a directory where intermediate files will be written. If none is given, a temporary directory will be used (in which case these intermediate files will be unavailable to the user).
    """
    # Get our controlling object.
    SOC = Calculator(calc_file = calc_file, calculation = calculation, num_singlets = num_singlets, num_triplets = num_triplets, QM_program = QM_program)
    
    # Compute SOC.
    SOC.calculate(output = output, SOC_scale = SOC_scale, include_ground = include_ground, CI_coefficient_threshold = CI_coefficient_threshold)
    
    # Finally, output results depending on requested format.
    if print_csv:
        # CSV.
        string_file = io.StringIO()
        
        # Get our CSV writer.
        writer = csv.writer(string_file)
        
        # Write.
        writer.writerows(SOC.soc_td.table)
        output_string = string_file.getvalue()
    else:
        # Text table.
        output_string = tabulate(
            SOC.soc_td.table,
            headers = "firstrow",
            numalign = "decimal",
            stralign = "center",
            floatfmt = "7.4f"
        )
        
    # Print.
    print(output_string)


# If we've been invoked as a program, call main().    
if __name__ == '__main__':
    # Get our args from the command line.
    parser = argparse.ArgumentParser(
        description = DESCRIPTION,
        epilog = EPILOG)
    
    # Add arguments.
    parser.add_argument("calc_file", help = "QM calculation output file to calculate SOC from", type = str)
    #parser.add_argument("-s", "--singlets", dest = "singlets", help = "The singlet excited states to calculate SOC for (eg, '-s 1 2 3' to calculate SOC for S1, S2 and S3).", type = int, nargs = "*", default = ())
    #parser.add_argument("-t", "--triplets", dest = "triplets", help = "The triplet excited states to calculate SOC for.", type = int, nargs = "*", default = ())
    parser.add_argument("--rwf_file", dest = "rwf_file_name", help = "Gaussian .rwf file. If this is not given explicitly, the location of this file will be guessed from the given QM calculation output file.", default = None, type = str)
    parser.add_argument("--fitted_basis", dest = "fitted_basis_directory", help = "Directory containing fitted basis set for use with DFTB+", type = str, default = None)
    parser.add_argument("-c", "--CSV", dest = "print_csv", help = "Output in CSV format", action = "store_true")
    parser.add_argument("-s", "--singlets", dest = "num_singlets", help = "The number of singlet excited states to calculate SOC for. The default is all available states.", type = int, default = None)
    parser.add_argument("-t", "--triplets", dest = "num_triplets", help = "The number of triplet excited states to calculate SOC for. The default is all available states.", type = int, default = None)
    parser.add_argument("-p", "--program", dest = "QM_program", help = "The QM program that the excited states were calculated with. If None is given, the program will be guessed from the given calc_file.", choices = ["Gaussian", "DFTB+"], default = None)
    parser.add_argument("-T", "--calculation", help = "The type of SOC calculation to perform, see the molsoc manual for more information. one: one-electron SOC; two: SOC with the full Breit–Pauli operator; zeff: one-electron SOC with the screened-nuclear charge method; auto: zeff if supported for the atoms in the given molecule, one otherwise", choices = ("one", "two", "zeff", "auto"), default = "auto")
    parser.add_argument("-S", "--SOC_scale", help = "Scaling factor for Zeff. If None is give, a default (1.0) will be used.", type = float, default = None)
    parser.add_argument("-o", "--output", help = "Path to a directory to write intermediate molsoc files to. If None is given, intermediate files will not be saved.", type = Path, default = None)
    parser.add_argument("-n", "--no_ground", dest = "include_ground", help = "Don't include the ground state in the SOC calculation.", action = "store_false")
    parser.add_argument("-C", "--CI_threshold", dest = "CI_coefficient_threshold", help = "Threshold for CI (CIS) coefficients.", type = float, default = None)    
    parser.add_argument("-v", "--version", action = "version", version = str(pysoc.version))
    
    # Call the main function, passing command line arguments as keyword arguments.
    sys.exit(
        main(
            **vars(parser.parse_args())
        )
    )
    
    