#!/usr/bin/env python2
'''Spin-Orbit Coupling(SOC) with interface to
   td-dft(Gaussian)/td-dftb(DFTB+, Thomas Niehaus)
   by Xing Gao, MPI@Muelheim, 2015
'''

# General imports.
import sys
import subprocess
from pysoc.io.dftb_plus import DFTB_plus_parser
sys.path.append('.')
from init import *

# Set path to molsoc.
molsoc_path = '/home/oliver/ownCloud/Chemistry/St. Andrews PhD/PySOC/PySOC Src/bin/molsoc0.1.exe'

# PySOC imports.
from pysoc.io.file import read_file, write_file, check_file_exist
from pysoc.io.gaussian import Gaussian_parser
                   
########################################################################
#                                                                      #
#                           main prepare                               #
#                                                                      #
########################################################################

print "QM_code", QM_code
if QM_ex_flag and QM_code == 'gauss_tddft':
   	'''prepare gaussian tddft calculation
      a. check gaussian input file *.com 
      b. read the geometry for the geom file by newton-X
      call gaussian tddft calculation
   	'''
   	pass

elif QM_ex_flag and QM_code == 'tddftb':
	pass

elif QM_ex_flag:
	print 'Only Gaussian tddft or TD-DFTB+ code are available now,,,'


#check the existence of Output file from third party quantum chemistry code
# like gaussian.log and gaussian.rwf

check_file_exist(qm_out)
#r = Read_QM_out()


if QM_code == 'gauss_tddft':

    try:
        singlet_energies, triplet_energies, singlet_levels, triplet_levels, nbov, MO_energies, MOA_coefficients, AO_overlaps, CI_coefficients = Gaussian_parser(qm_out[0], qm_out[1], n_s, n_t).parse(molsoc_input[0], molsoc_input[1], soc_key, soc_scal)
        #r.read_gauss(qm_out, soc_key, soc_scal, molsoc_input)
        
        # Get our calculation parser.
        #parser = Gaussian_parser(qm_out[0], qm_out[1], n_s, n_t)
        #parser.parse(molsoc_input[0], molsoc_input[1], soc_key, soc_scal)

    except Exception as e:
        raise Exception("Error when reading Gaussian output")

elif QM_code == 'tddftb':
    print "output from td-dftb+ is called..." 
    check_file_exist(geom_xyz)
    try:
        singlet_energies, triplet_energies, singlet_levels, triplet_levels, nbov, ao_ncart,\
        MO_energies, MOA_coefficients, AO_overlaps, CI_coefficients = \
        DFTB_plus_parser(n_s, n_t).parse(qm_out, soc_key, soc_scal, molsoc_input, geom_xyz, dir_para_basis)
        #r.read_tdtb(qm_out, soc_key, soc_scal, molsoc_input, geom_xyz,\
        #           dir_para_basis)

    except Exception as e:
        raise Exception("Error when reading td-dftb+ output")
       

########################################################################
#save file on disk
write_file(MO_energies, 'mo_ene.dat', '{}  ')
write_file(MOA_coefficients, 'mo_coeff.dat', '{}  ')
write_file(AO_overlaps, 'ao_overlap.dat', '{}  ')
write_file(CI_coefficients, 'ci_coeff.dat', '{}  ')
#

########################################################################
#effective single electron SOC calculation
#on atomic basis
#check the existence of molsoc.inp and molsoc_basis
########################################################################
check_file_exist(molsoc_input)

#molsoc_retval = subprocess.call(molsoc_path+' '+molsoc_input[0], shell=True)
molsoc_retval = subprocess.call([molsoc_path, molsoc_input[0]])
if molsoc_retval != 0:
    raise Exception("molsoc failed with error code {}".format(molsoc_retval))

if QM_code == 'tddftb':
    nt = ao_ncart[0]**2
else:
    nt = nbov[0]**2
    
soint = []
hso = dict(hso_x='X COMPONENT', hso_y='Y COMPONENT', hso_z='Z COMPONENT')
for i, cont0 in sorted(hso.iteritems()):
    hso[i] = read_file('soint', cont0, nt, 0)
    hso[i] = [cont1 for k, cont1 in enumerate(hso[i]) if (k-2)%3==0]
    #soint = soint + hso[i] + ['\n']
    soint = soint + hso[i] 
print "len(soint)", len(soint)
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
    for i, cont0 in sorted(dip.iteritems()):
        dip[i] = read_file('molsoc_dipole.dat', cont0, nt, 0)
        #dipint = dipint + dip[i] + ['\n']
        dipint = dipint + dip[i] 

write_file(dipint, 'dip_ao.dat', '{}  ')
########################################################################
#td soc
#
n_total_s = [len(n_s)] 
n_total_t = [len(n_t)]
dat_out = dict(a0_code=[QM_code], a0_dipole=dipole_flag, a0_num_g=n_g, 
               a0_thresh=cicoeff_thresh, a0_total_s=n_total_s,
               a1_num_s=n_s, a1_order_s=singlet_levels,
               a2_ene_s=singlet_energies, a0_total_t=n_total_t, 
               b1_num_t=n_t, b1_order_t=triplet_levels, b2_ener_t=triplet_energies, 
               c_num_bov=nbov)
if QM_code == 'tddftb':
    dat_out.update({'d_basis_ncart': ao_ncart}) 

with open('soc_td_input.dat', 'w') as f:
    for i, line in sorted(dat_out.iteritems()):
        for k in line:
            f.write('{:<15}'.format(k))
        f.write('{:>4}\n'.format(i))

#
soc_cis = './soc_td'

#p = subprocedss.Popen(call_soc_cis, shell=True, stdout=subprocess.PIPE)
#subprocess.call(soc_cis+' '+test, shell=True)
########################################################################
















