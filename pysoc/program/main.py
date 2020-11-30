#!/usr/bin/env python3.6
'''Spin-Orbit Coupling(SOC) with interface to
   td-dft(Gaussian)/td-dftb(DFTB+, Thomas Niehaus)
   by Xing Gao, MPI@Muelheim, 2015
   
   Rewritten for python3 by Oliver Lee
'''

# General imports.
import sys
import tempfile
import argparse
from pathlib import Path
from tabulate import tabulate

# Set path to molsoc.
#molsoc_path = '/home/oliver/ownCloud/Chemistry/St. Andrews PhD/PySOC/PySOC Src/bin/molsoc0.1.exe'

# PySOC imports.
import pysoc
from pysoc.io.gaussian import Gaussian_parser
from pysoc.io.dftb_plus import DFTB_plus_parser
from pysoc.io.soc_td import Soc_td
import io
import csv

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
    # TODO: Because of issue #1, we don't allow selecting exact singlets/triplets, we just ask how many.
    singlets = list(range(1, num_singlets +1)) if num_singlets is not None else None
    triplets = list(range(1, num_triplets +1)) if num_triplets is not None else None
    
    # Set defaults if not given.
    if include_ground is None:
        include_ground = True
        
    if CI_coefficient_threshold is None:
        CI_coefficient_threshold = 1.0e-5
        
    if print_csv is None:
        print_csv = False
    
    # First, get a temp dir if we need one.
    with tempfile.TemporaryDirectory() as tempdir:
        # Only use if necessary.
        if output is None:
            output = tempdir
            
        # If we weren't told what QM_program to use, guess from the calc_file.
        calc_file = Path(calc_file)
        if QM_program is None:
            # Try and guess from the input file type.
            if calc_file.suffix.lower() == ".log":
                QM_program = "Gaussian"
            elif calc_file.suffix.lower() == ".xyz":
                QM_program = "DFTB+"
            else:
                raise Exception("Could not guess input program type from file '{}'; try specifying explicitly with '--program'".format(calc_file))
            
    
        # Now we need to parse the output from our QM program.
        # Get an appropriate parser.
        if QM_program == 'Gaussian':
            # Get our calculation parser.
            molsoc = Gaussian_parser.from_output_files(calc_file, requested_singlets = singlets, requested_triplets = triplets, **aux_files)
            
            # Keywords for molsoc
            keywords = ('ANG', 'Zeff', 'DIP')
        
        elif QM_program == 'DFTB+':
            # Get our calculation parser.
            molsoc = DFTB_plus_parser.from_output_files(calc_file, requested_singlets = singlets, requested_triplets = triplets, **aux_files)
            
            # Keywords for molsoc
            keywords = ('ANG', 'Zeff', 'DIP', 'TDB')
        
        else:
            # We were given something random.
            raise Exception("Unknown or unrecognised program name '{}'".format(QM_program))
        
        # Parse and prepare input for molsoc.
        molsoc.parse()
        molsoc.prepare(keywords, SOC_scale, output)
        # Run molsoc.
        molsoc.run()
        
        
        
        # Prepare input for soc_td.
        soc_td = Soc_td(molsoc)
        soc_td.prepare(keywords, include_ground, CI_coefficient_threshold)
        
        
        # Now call soc_td.
        soc_td.run()
        
        # Finally, output results depending on requested format.
        if print_csv:
            # CSV.
            string_file = io.StringIO()
            
            # Get our CSV writer.
            writer = csv.writer(string_file)
            
            # Write.
            writer.writerows(soc_td.table)
            output_string = string_file.getvalue()
        else:
            # Text table.
            output_string = tabulate(
                soc_td.table,
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
    parser.add_argument("-c", "--CSV", dest = "print_csv", help = "Output in CSV format", action = "store_true")
    parser.add_argument("-s", "--singlets", dest = "num_singlets", help = "The number of singlet excited states to calculate SOC for. The default is all available states.", type = int, default = None)
    parser.add_argument("-t", "--triplets", dest = "num_triplets", help = "The number of triplet excited states to calculate SOC for. The default is all available states.", type = int, default = None)
    parser.add_argument("-p", "--program", dest = "QM_program", help = "The QM program that the excited states were calculated with. If None is given, the program will be guessed from the given calc_file.", choices = ["Gaussian", "DFTB+"], default = None)
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
    
    