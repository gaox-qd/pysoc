#!/usr/bin/env python3.6
'''Spin-Orbit Coupling(SOC) with interface to
   td-dft(Gaussian)/td-dftb(DFTB+, Thomas Niehaus)
   by Xing Gao, MPI@Muelheim, 2015
'''

# General imports.
import sys
import subprocess
sys.path.append('.')
from init import *

# Set path to molsoc.
molsoc_path = '/home/oliver/ownCloud/Chemistry/St. Andrews PhD/PySOC/PySOC Src/bin/molsoc0.1.exe'

# PySOC imports.
from pysoc.io.file import read_file, write_file, check_file_exist
from pysoc.io.gaussian import Gaussian_parser
from pysoc.io.dftb_plus import DFTB_plus_parser
                   
def main():
    """
    Main program function for PySOC controller program.
    
    
    """
    # Now we need to parse the output from our QM program.
    # Get an appropriate parser.
    #parser = parser_from_program(QM_program)
    if QM_program == 'Gaussian':
        # Get our calculation parser.
        parser = Gaussian_parser(qm_out[0], qm_out[1], n_s, n_t)
        parser.parse()
        #parser.prepare_molsoc_input(molsoc_input[0], molsoc_input[1], soc_key, soc_scal)
        parser.prepare_molsoc_input(soc_key, soc_scal, molsoc_input[0], molsoc_input[1], "ao_basis.dat")
    
    elif QM_program == 'DFTB+':
        # Get our calculation parser.
        parser = DFTB_plus_parser(n_s, n_t)
        parser.parse(qm_out, soc_key, soc_scal, molsoc_input, geom_xyz, dir_para_basis)
    
    else:
        # We were given something random.
        raise Exception("Unknown or unrecognised program name '{}'".format(QM_program))
    
    
    ########################################################################
    #save file on disk
    write_file(parser.MO_energies, 'mo_ene.dat', '{}  ')
    write_file(parser.MOA_coefficients, 'mo_coeff.dat', '{}  ')
    write_file(parser.AO_overlaps, 'ao_overlap.dat', '{}  ')
    write_file(parser.CI_coefficients, 'ci_coeff.dat', '{}  ')
    
    ########################################################################
    
    molsoc_retval = subprocess.call([molsoc_path, molsoc_input[0]])
    if molsoc_retval != 0:
        raise Exception("molsoc failed with error code {}".format(molsoc_retval))
    
    if QM_program == 'DFTB+':
        nt = parser.ao_ncart[0] **2
    else:
        nt = parser.num_orbitals **2
        
    soint = []
    hso = dict(hso_x='X COMPONENT', hso_y='Y COMPONENT', hso_z='Z COMPONENT')
    for i, cont0 in sorted(hso.items()):
        hso[i] = read_file('soint', cont0, nt, 0)
        hso[i] = [cont1 for k, cont1 in enumerate(hso[i]) if (k-2)%3==0]
        #soint = soint + hso[i] + ['\n']
        soint = soint + hso[i] 
    
    write_file(soint, 'soc_ao.dat', '{}  ')
    #print hso['hso_z']
    #print f_num.findall(hso_x)
    ########################################################################
    s_matr = read_file('molsoc_overlap.dat', 'AO_overlap', nt, 0)
    write_file(s_matr, 's_matr.dat', '{}  ')
    ########################################################################
    #read dipole moment on atomic basis
    dipole_flag = ['False']
    if 'DIP' in soc_key:
        dipole_flag = ['True']
        dipint = []
        dip = dict(dip_x='DIM=1', dip_y='DIM=2', dip_z='DIM=3')
        for i, cont0 in sorted(dip.items()):
            dip[i] = read_file('molsoc_dipole.dat', cont0, nt, 0)
            #dipint = dipint + dip[i] + ['\n']
            dipint = dipint + dip[i] 
    
    write_file(dipint, 'dip_ao.dat', '{}  ')
    ########################################################################
    #td soc
    #
    
    dat_out = dict(a0_code = ["gauss_tddft" if  QM_program == "Gaussian" else "tddftb"],
                   a0_dipole = dipole_flag,
                   a0_num_g = n_g, 
                   a0_thresh = cicoeff_thresh,
                   a0_total_s = [len(n_s)],
                   a1_num_s = n_s,
                   a1_order_s = parser.singlet_levels,
                   a2_ene_s = parser.singlet_energies,
                   a0_total_t = [len(n_t)], 
                   b1_num_t = n_t,
                   b1_order_t = parser.triplet_levels,
                   b2_ener_t = parser.triplet_energies, 
                   c_num_bov = [parser.num_orbitals, parser.num_occupied_orbitals, parser.num_virtual_orbitals])
    if QM_program == 'DFTB+':
        dat_out.update({'d_basis_ncart': parser.ao_ncart}) 
    
    with open('soc_td_input.dat', 'w') as f:
        for i, line in sorted(dat_out.items()):
            for k in line:
                f.write('{:<15}'.format(k))
            f.write('{:>4}\n'.format(i))
    
    #
    soc_cis = './soc_td'


# If we've been invoked as a program, call main().    
if __name__ == '__main__':
    sys.exit(main())