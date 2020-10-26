#!/usr/bin/env python2
'''Spin-Orbit Coupling(SOC) with interface to
   td-dft(Gaussian)/td-dftb(DFTB+, Thomas Niehaus)
   by Xing Gao, MPI@Muelheim, 2015
'''
import sys 
#import init as io
import os
import re
import subprocess
from itertools import islice
sys.path.append('.')
from init import *


#parameter
au_2_angs = (0.52917721092) # from wikipedia
angs_2_au = (1.88972612456) # 1/Bohr


class Prep_Basis:
   pass
   def __init__(self):
       print "prepare basis set"

#test = QM_ex()
'''read the output from QM calculation for the:
   a. # of basis set, virt orb, occ orb
   b. desired excitation energy 
   c. generate basis set
   d. molecular orbital, orbital energy
   e. CIS coefficients
'''
f_num = re.compile(r'[-+]?\d*\.\d+|\d+')

#######################################################################
def read_file(filein='fileinput', sign='where_start',
              nt=1, nr=1):
    dat_out = []
    with open(filein, 'r') as f:
        for nline, line in enumerate(f):
            if sign in line:
                break
        nf = nt - 1 
        #print filein, nf, nt, nline
        for nline, line in enumerate(f):
            dat_out.extend(line.split())
            if nline == nf:
                break
        #print 'endline', line
        if nr != 0:
            #print 'nr=', nr
            for line in f:
               print line.split()[0:nr]
               dat_out.extend(line.split()[0:nr])
               break
    return dat_out

#######################################################################
def write_file(dat_in=[], fileout='out.dat',style='{}'):
       with open(fileout, 'w' ) as f:
           for line in dat_in:
               f.write(style.format(line))

#######################################################################

class Read_QM_out:
   

   def read_gauss(self, outfile, soc_key, soc_scal, molsoc_input):
       ##gaussian.log
       print 'open file {}'.format(outfile[0])
       with open(outfile[0], 'r') as f:
           #nbov = [nbasis, nocc, nvirt]
           nbov = []
           ene_t, ene_s = [], []
           nbs = []
           ngeom, element, geom  = [], [], []
           for nf, line in enumerate(f):
               if all(s in line for s in ['Charge','Multiplicity']):
                   ngeom = nf + 1
               elif 'NAtoms=' in line:
                   natom = int(line.split()[1])
               elif 'NFC=' in line:
                   a = int(line.split()[1])
                   nfc = int(line.split()[7])
                   nbov.append(a)
               elif 'NVA=' in line:
                   a = int(line.split()[3])
                   b = int(line.split()[7])
                   nbov.extend((a,b))
               elif all(s in line for s in ['Excited State','Triplet']):
                   num = f_num.findall(line)
                   # state number and excitation energy for triplet
                   a = [int(num[0]), float(num[1])]
                   ene_t.append(a)
               elif all(s in line for s in ['Excited State','Singlet']):
                   # state number and excitation energy for singlet
                   num = f_num.findall(line)
                   a = [int(num[0]), float(num[1])]
                   ene_s.append(a) #root num and excited energy
               elif any(s in line for s in 
                        ['AO basis set','primitive gaussians']):
                   nbs.append(nf)
                   #break
       print ene_t
       print ene_s[0][0], ene_s[0][1] 

       ###prepare input for molsoc
       with open(outfile[0], 'r') as f:
           lines = list(islice(f, ngeom, ngeom+natom))
           element = [x.split()[0] for x in lines]
           geom = [map(float,x.split()[1:4]) for x in lines]
           lines.append('End')
       with open(molsoc_input[0], 'w') as f:
           f.write('#input for soc in atomic basis\n')
           print soc_key
           for k in soc_key:
               f.write('{}  '.format(k))
           f.write('\n')
           f.write('#\n')
           for i, xyz in enumerate(geom):
               f.write('{:<4}{:>15.5f}{:>15.5f}{:>15.5f}{:>7.2f}\n'
                      .format(element[i], xyz[0], xyz[1], xyz[2], soc_scal))
           f.write('End')
       
       ###read basis set info
       shell = dict(S=1, P=3, SP=4, D=6, F=10 )
       ao_basis = []
       with open(outfile[0], 'r') as f:
           lines = list(islice(f, nbs[0]+1, nbs[1]-1))
           lines.pop()
           lines.append('END')
           #print lines
       with open(molsoc_input[1], 'w') as f:
           i = 0
           for line in lines:
               if line.split()[0] in shell:
                  sub_orb = line.split()[0]
                  pair = [sub_orb, shell[sub_orb]]
                  ao_basis.append(pair)
               if len(line.split()) == 2 and re.search(r'\d \d', line):
                   line = '{}  {}\n'.format(element[i], line.split()[1])
                   i += 1
               f.write(line)
       print 'ao_basis', ao_basis
       nc_basis = len(ao_basis)
       n_basis = sum(ao_basis[k][1] for k in range(nc_basis)) 
       print "n_basis", n_basis
       if n_basis != nbov[0]:
           print "n_basis is not eq to nbov[0]!!", n_basis, nbov[0]
           sys.exit()
       with open('ao_basis.dat', 'w') as f:
           print "hello1"
           f.write('{}  {}\n'.format(nc_basis, n_basis))
           print "hello2"
           for i in range(nc_basis):
               f.write('{}  '.format(ao_basis[i][1]))
            

       ##gaussian.rwf
       print 'open file {}'.format(outfile[1])
       with open(outfile[1], 'r') as f:
          cmd0 = g09root+'/g09/bsd/g09.profile'
          cmd1 = g09root+'/g09/rwfdump gaussian.rwf'
          out_dat = dict(MO_energy=' 522R', MOA_coeffs=' 524R', 
                         MOB_coeffs=' 526R', XY_coeffs=' 635R',
                         AO_overlap=' 514R' )
          for key, cont in out_dat.iteritems():
               print key+cont
               subprocess.call(cmd0+';'+cmd1+' '+key+cont, shell=True)

       
       mo_ene, ao_overlap, mo_coeff, ci_coeff = [], [], [], []

       ###reading MO_energy
       nt = nbov[0] / 5
       nr = nbov[0] % 5
       mo_ene = read_file('MO_energy', 'Dump of file', nt, nr )
       #nbi, nbf = nbov[0]-nbov[1]-nbov[2], nbov[0]
       nbi, nbf = nfc, nbov[0]
       mo_ene = mo_ene[nbi:nbf]

       ###reading AO_overlap
       nao = nbov[0] * (nbov[0]+1) /2
       nt  = nao / 5
       nr  = nao % 5
       ao_overlap = read_file('AO_overlap', 'Dump of file', nt, nr)

       ###reading MOA_coeffs
       mo_coeff = read_file('MOA_coeffs', 'Dump of file', -2, 0 )
       #print 'mo_coeff1', len(mo_coeff)
       nbi, nbf = nbi*nbov[0], nbov[0]**2
       mo_coeff = mo_coeff[nbi:nbf]
       #print 'mo_coeff1', len(mo_coeff)
       
       ###reading XY_coeffs
       order1 = [ene_s[i-1][0] for i in n_s]
       order2 = [ene_t[i-1][0] for i in n_t]
       order0 = order1 + order2
       e_s = [ene_s[i-1][1] for i in n_s]
       e_t = [ene_t[i-1][1] for i in n_t]
       print 'e_s', e_s
       print 'e_t', e_t
       print 'total order for the roots:', order0
       max_root = max(order0)
       #print max_root
       dim = nbov[1] * nbov[2]
       with open('XY_coeffs', 'r') as f:
           for line in f:
               if 'Dump of file' in line:
                   dat_lenth = int(f_num.findall(line)[1])
                   break
       mseek = (dat_lenth-12) / (dim*4+1)
       print mseek
       nline = dim * 2 * (max_root+mseek) + 12
       nt = nline / 5
       nr = nline % 5
       ci_coeff = read_file('XY_coeffs', 'Dump of file', nt, nr)
       ci_xpy, ci_xmy = [], []
       for i in order0:   #singlets come first, then triplets
           np = 12 + dim * 2 * (i-1) 
           #ci_xpy += ci_coeff[np:np+dim] #X+Y(xpy)
           ci_xpy += ci_coeff[np:np+dim*2] #X+Y(xpy)
           np = 12 + dim * 2 * (mseek+i-1) #X-Y(xmy)
           #ci_xmy += ci_coeff[np:np+dim]
           ci_xmy += ci_coeff[np:np+dim*2]
       ci_coeff = ci_xpy + ci_xmy 
       print len(ci_xpy), len(ci_xmy)
       print len(ci_coeff)
       
       return e_s, e_t, order1, order2, nbov, \
              mo_ene, mo_coeff, ao_overlap, ci_coeff
            
   
########################################################################
   def read_tdtb(self, outfile, soc_key, soc_scal, molsoc_input, \
                 geom_xyz, dir_para_basis):
       ###reading MO_energy, # of mo, occ, virt
       print "reading MO_energy, # of mo, occ, virt..." 
       mo_ene, nbov = [], []
       var = []
       with open(outfile[0], 'r') as f:
           for nline, line in enumerate(f):
               #print nline
               if len(line.split()) > 1:
                   var = line.split()[0]
                   #print "var_temp", var
                   if f_num.match(var):
                       mo_ene.append(var)
                       if float(line.split()[1]) - 1.0 > 1.0E-5:
                           nb = nline
                           #print nb
           nbov.extend((nline-1, nb, nline-1-nb))
       print "mo_ene", mo_ene
       print "nbov", nbov
       ###reading singlet and triplet energies
       print "reading singlet and triplet energies..."
       ene_t, ene_s = [], []
       iit, iis = 0, 0
       with open(outfile[1], 'r') as f:
           for line in f:
               if len(line.split()) > 1:
                   var = line.split()[0]
                   if f_num.match(var):
                       if 'T' in line.split():
                           iit = iit + 1
                           a = [iit, float(var)]
                           ene_t.append(a)
                       elif 'S' in line.split():
                           iis = iis + 1
                           a = [iis, float(var)]
                           ene_s.append(a)
       print "ene_t", ene_t
       print "ene_s", ene_s
       order1 = [ene_s[i-1][0] for i in n_s]
       order2 = [ene_t[i-1][0] for i in n_t]
       order0 = order1 + order2
       e_s = [ene_s[i-1][1] for i in n_s]
       e_t = [ene_t[i-1][1] for i in n_t]
       print 'e_s', e_s
       print 'e_t', e_t
       print 'total order for the roots:', order0
       ###reading oversqr.dat (ao overlap)
       print "reading oversqr.dat (ao overlap)..."
       ao_overlap = []
       with open(outfile[2], 'r') as f:
           for line in f:
               if len(line.split()) == nbov[0]:
                   ao_overlap.extend(line.split())
       print "len(ao_overlap)", len(ao_overlap)
       ###reading eigenvec.out(mo coeffs)
       ###and save above basis set info  
       print "reading eigenvec.out(mo coeffs)..."
       mo_coeff = []
       ao_basis = []
       orb = ['s1', 'p1', 'p2', 'p3', \
              'd1', 'd2', 'd3', 'd4', 'd5', \
              'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7']
       max_shell = ['s1', 'p3', 'd5', 'f7']
       with open(outfile[3], 'r') as f:
           for line in f:
               a = line.split()
               if len(a) > 1:
                   if a[0] in orb:
                       mo_coeff.append(a[1])
       with open(outfile[3], 'r') as f:
           for line in f:
               a = line.split()
               if len(a) > 1:
                   if a[0] in max_shell:
                       if a[0] == 'd5':
                           ao_basis.append('6')
                       elif a[0] == 'f7':
                           ao_basis.append('10')
                       else:
                           ao_basis.append(list(a[0])[1])
                   elif all(s in a for s in ['Eigenvector:', '2']):
                       break
       
       n_basis = sum(int(io) for io in ao_basis) # the # of Cartesian AOs
       nc_basis = len(ao_basis)

       ao_virt =n_basis - nbov[1]
       ao_ncart  = [n_basis, nbov[1], ao_virt]
       print "ao_ncart", ao_ncart
       print "len,nbasis", nc_basis, n_basis
       print ao_basis       
       with open('ao_basis.dat', 'w') as f:
           f.write('{}  {}\n'.format(nc_basis, n_basis))
           for i in range(nc_basis):
               f.write('{}  '.format(ao_basis[i]))
       
       ###reading XplusY.DAT(X+Y coeffs) 
       ###and reorder XplusY according SPX.DAT 
       print "reading XplusY.DAT(X+Y coeffs)..."
       ci_coeff, ci_xpy = [], []
       ndim = nbov[1] * nbov[2]
       print "nidm", ndim
       with open(outfile[4], 'r') as f:
           for line in f:
               if(line.split()[1] not in ['S','T'] and 
                  float(line.split()[0]) != ndim): #6 data each line
                   ci_coeff.extend(line.split())
       print "ci_coeff", ci_coeff[0], ci_coeff[-1]
       trans_key = [] #occ->virt transition order
       with open(outfile[5], 'r') as f:
           for line in f:
               if len(line.split()) == 6: #6 data each line
                   #chr: convert num to character
                   trans_temp = [chr(int(line.split()[i])) for i in [3, 5]]
                   trans_key.append("".join(trans_temp))

       #ndim = nbov[1] * nbov[2]
       #print "nidm", ndim
       nt = len(ene_t) #num of triplets come before singlets
       for i in n_s:
           np = (nt+i-1) * ndim
           trans_dict = {trans_key[k]: ci_coeff[np+k] for k in range(ndim)}
           for k, v in sorted(trans_dict.iteritems()):
               ci_xpy.append(v)
               #print "trans_ci_singlet:",[ord(o) for o in k.split(' ')]
       for i in n_t:
           np = (i-1) * ndim
           trans_dict = {trans_key[k]: ci_coeff[np+k] for k in range(ndim)}
           for k, v in sorted(trans_dict.iteritems()):
               ci_xpy.append(v)
               #print "trans_ci_triplet:",[ord(o) for o in k.split(' ')]
       ci_coeff = ci_xpy
       print len(ci_xpy) 
       #print ci_xpy
       

       ###prepare input for molsoc
       print "prepare input for molsoc..."
       ngeom = 1 # coordinate after line 1 in xyz file
       print geom_xyz[0]
       with open(geom_xyz[0], 'r') as f:
           for line in f:
               natom = int(f_num.findall(line)[0])
               break
           print "natom", natom
           lines = list(islice(f, ngeom, ngeom+natom)) 
           element = [x.split()[0] for x in lines]
           geom = [map(float,x.split()[1:4]) for x in lines]
           #print geom
           print element

       with open(molsoc_input[0], 'w') as f:
           print "hello", molsoc_input[0]
           f.write('#input for soc in atomic basis\n')
           print soc_key
           for k in soc_key:
               f.write('{}  '.format(k))
           f.write('\n')
           f.write('#\n')
           for i, xyz in enumerate(geom):
               f.write('{:<4}{:>15.5f}{:>15.5f}{:>15.5f}{:>7.2f}\n'
                      .format(element[i], xyz[0], xyz[1], xyz[2], soc_scal))
           f.write('End')
               
       ###build basis set 
       print "build basis set..."
       with open(molsoc_input[1], 'w') as fw:
           print "dir_para_basis", dir_para_basis
           for i, el in enumerate(element):
               bas_f = dir_para_basis+'/'+el+'.basis'
               #print "bas_f", bas_f
               with open(bas_f, 'r') as fr:
                   lines = fr.readlines()
                   fw.writelines(lines)
               if i != natom-1:
                   fw.write(' ****\n')
               else:
                   fw.write('END')


       return e_s, e_t, order1, order2, nbov, ao_ncart, \
              mo_ene, mo_coeff, ao_overlap, ci_coeff
                   
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

def check_file_exist(filename=[], err_inf='not exist'):
   for f in filename: 
       isfile = os.path.isfile(f)
       if not isfile:
           print ("file {} is " + err_inf).format(f)
           sys.exit()

#check the existence of Output file from third party quantum chemistry code
# like gaussian.log and gaussian.rwf

check_file_exist(qm_out)
r = Read_QM_out()

if QM_code == 'gauss_tddft':

   try:
       e_s, e_t, order_s, order_t, nbov, \
       mo_ene, mo_coeff, ao_overlap, ci_coeff = \
       r.read_gauss(qm_out, soc_key, soc_scal, molsoc_input)

       print e_s
       print order_s
       print e_t
       print order_t
       print nbov
   except:
       print "Error when reading gaussian output"
   #
elif QM_code == 'tddftb':
   print "output from td-dftb+ is called..." 
   check_file_exist(geom_xyz)
   try:
       e_s, e_t, order_s, order_t, nbov, ao_ncart,\
       mo_ene, mo_coeff, ao_overlap, ci_coeff = \
       r.read_tdtb(qm_out, soc_key, soc_scal, molsoc_input, geom_xyz,\
                  dir_para_basis)

       print e_s
       print order_s
       print e_t
       print order_t
       print nbov
   except:
       print "Error when reading td-dftb+ output"
   #
       

########################################################################
#save file on disk
write_file(mo_ene, 'mo_ene.dat', '{}  ')
write_file(mo_coeff, 'mo_coeff.dat', '{}  ')
write_file(ao_overlap, 'ao_overlap.dat', '{}  ')
write_file(ci_coeff, 'ci_coeff.dat', '{}  ')
#

########################################################################
#effective single electron SOC calculation
#on atomic basis
#check the existence of molsoc.inp and molsoc_basis
########################################################################
check_file_exist(molsoc_input, 'not generated successfully,,,')
subprocess.call(molsoc_path+' '+molsoc_input[0], shell=True)

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
               a1_num_s=n_s, a1_order_s=order_s,
               a2_ene_s=e_s, a0_total_t=n_total_t, 
               b1_num_t=n_t, b1_order_t=order_t, b2_ener_t=e_t, 
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
















