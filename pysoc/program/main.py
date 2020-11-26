#!/usr/bin/env python3.6
'''Spin-Orbit Coupling(SOC) with interface to
   td-dft(Gaussian)/td-dftb(DFTB+, Thomas Niehaus)
   by Xing Gao, MPI@Muelheim, 2015
   
   Rewritten for python3 by Oliver Lee
'''

# General imports.
import sys
import subprocess
import tempfile
import argparse
from pathlib import Path

# Set path to molsoc.
molsoc_path = '/home/oliver/ownCloud/Chemistry/St. Andrews PhD/PySOC/PySOC Src/bin/molsoc0.1.exe'

# PySOC imports.
import pysoc
from pysoc.io.file import read_file, write_file
from pysoc.io.gaussian import Gaussian_parser
from pysoc.io.dftb_plus import DFTB_plus_parser

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
        #parser = parser_from_program(QM_program)
        if QM_program == 'Gaussian':
            # Get our calculation parser.
            #parser = Gaussian_parser(qm_out[0], qm_out[1], num_singlets, num_triplets)
            parser = Gaussian_parser.from_output_files(calc_file, requested_singlets = singlets, requested_triplets = triplets, **aux_files)
            #parser.parse()
            
            # Keywords for molsoc
            keywords = ('ANG', 'Zeff', 'DIP')
            
            #molsoc_input_file = parser.prepare_molsoc_input(keywords, SOC_scale, output)[0]
        
        elif QM_program == 'DFTB+':
            # Get our calculation parser.
            parser = DFTB_plus_parser.from_output_files(calc_file, requested_singlets = singlets, requested_triplets = triplets, **aux_files)
            
            # Keywords for molsoc
            keywords = ('ANG', 'Zeff', 'DIP', 'TDB')
            
            #parser.parse(qm_out, keywords, SOC_scale, molsoc_input, geom_xyz, dir_para_basis)
            #molsoc_input_file = parser.prepare_molsoc_input(keywords, SOC_scale, output)[0]
        
        else:
            # We were given something random.
            raise Exception("Unknown or unrecognised program name '{}'".format(QM_program))
        
        
        # Parse and prepare.
        parser.parse()
        molsoc_input_file = parser.prepare_molsoc_input(keywords, SOC_scale, output)[0]
        
        
        ########################################################################
        #save file on disk
        write_file(parser.MO_energies, Path(output, 'mo_ene.dat'), '{}  ')
        write_file(parser.MOA_coefficients, Path(output, 'mo_coeff.dat'), '{}  ')
        write_file(parser.AO_overlaps, Path(output, 'ao_overlap.dat'), '{}  ')
        write_file(parser.CI_coefficients, Path(output, 'ci_coeff.dat'), '{}  ')
        
        ########################################################################
        
#         molsoc_retval = subprocess.call([molsoc_path, molsoc_input_file])
#         if molsoc_retval != 0:
#             raise Exception("molsoc failed with error code {}".format(molsoc_retval))
        subprocess.run(
            (molsoc_path, molsoc_input_file.resolve()),
            check = True,
            universal_newlines = True,
            cwd = output
        )
        
        if QM_program == 'DFTB+':
            nt = parser.ao_ncart[0] **2
        else:
            nt = parser.num_orbitals **2
            
        soint = []
        hso = dict(hso_x='X COMPONENT', hso_y='Y COMPONENT', hso_z='Z COMPONENT')
        for i, cont0 in sorted(hso.items()):
            hso[i] = read_file(Path(output, 'soint'), cont0, nt, 0)
            hso[i] = [cont1 for k, cont1 in enumerate(hso[i]) if (k-2)%3==0]
            #soint = soint + hso[i] + ['\n']
            soint = soint + hso[i] 
        
        write_file(soint, Path(output, 'soc_ao.dat'), '{}  ')
        #print hso['hso_z']
        #print f_num.findall(hso_x)
        ########################################################################
        s_matr = read_file(Path(output, 'molsoc_overlap.dat'), 'AO_overlap', nt, 0)
        write_file(s_matr, Path(output, 's_matr.dat'), '{}  ')
        ########################################################################
        #read dipole moment on atomic basis
        dipole_flag = ['False']
        if 'DIP' in keywords:
            dipole_flag = ['True']
            dipint = []
            dip = dict(dip_x='DIM=1', dip_y='DIM=2', dip_z='DIM=3')
            for i, cont0 in sorted(dip.items()):
                dip[i] = read_file(Path(output, 'molsoc_dipole.dat'), cont0, nt, 0)
                #dipint = dipint + dip[i] + ['\n']
                dipint = dipint + dip[i] 
        
        write_file(dipint, Path(output, 'dip_ao.dat'), '{}  ')
        ########################################################################
        #td soc
        #
        
        dat_out = dict(a0_code = ["gauss_tddft" if  QM_program == "Gaussian" else "tddftb"],
                       a0_dipole = dipole_flag,
                       a0_num_g = ["True" if include_ground else "False"], 
                       a0_thresh = [CI_coefficient_threshold],
                       a0_total_s = [len(parser.requested_singlets)],
                       a1_num_s = parser.requested_singlets,
                       a1_order_s = parser.singlet_levels,
                       a2_ene_s = parser.singlet_energies,
                       a0_total_t = [len(parser.requested_triplets)], 
                       b1_num_t = parser.requested_triplets,
                       b1_order_t = parser.triplet_levels,
                       b2_ener_t = parser.triplet_energies, 
                       c_num_bov = [parser.num_orbitals, parser.num_occupied_orbitals, parser.num_virtual_orbitals])
        if QM_program == 'DFTB+':
            dat_out.update({'d_basis_ncart': parser.ao_ncart}) 
        
        with open(Path(output, 'soc_td_input.dat'), 'w') as f:
            for i, line in sorted(dat_out.items()):
                for k in line:
                    f.write('{:<15}'.format(k))
                f.write('{:>4}\n'.format(i))
        
        # Now call soc_td.
        subprocess.run(
            ("soc_td",),
            cwd = output,
            check = True,
            universal_newlines =True
        )
        
        # Finally, read and output the SOC values.
        with open(Path(output, "soc_out.dat"), "r") as soc_file:
            print(soc_file.read())


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
    parser.add_argument("-s", "--singlets", dest = "num_singlets", help = "The number of singlet excited states to calculate SOC for. The default is all available states.", type = int, default = None)
    parser.add_argument("-t", "--triplets", dest = "num_triplets", help = "The number of triplet excited states to calculate SOC for. The default is all available states.", type = int, default = None)
    parser.add_argument("-p", "--program", dest = "QM_program", help = "The QM program that the excited states were calculated with. If None is given, the program will be guessed from the given calc_file.", choices = ["Gaussian", "DFTB+"], default = None)
    parser.add_argument("-S", "--SOC_scale", help = "Scaling factor for Zeff. If None is give, a default (1.0) will be used.", type = float, default = None)
    parser.add_argument("-o", "--output", help = "Path to a directory to write intermediate molsoc files to. If None is given, intermediate files will not be saved.", type = Path, default = None)
    parser.add_argument("-n", "--no_ground", dest = "include_ground", help = "Don't include the ground state in the SOC calculation.", action = "store_false")
    parser.add_argument("-C", "--CI_threshold", dest = "CI_coefficient_threshold", help = "Threshold for CI (CIS) coefficients.", type = float, default = None)    
    parser.add_argument("-v", "--version", action = "version", version = str(pysoc.version))
    
    sys.exit(
        main(
            **vars(parser.parse_args())
        )
    )
    
    