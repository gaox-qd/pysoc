from pysoc.io import Output_parser
from itertools import islice


class DFTB_plus_parser(Output_parser):
    """
    Class for parsing the required output from Gaussian.
    """
    
    def __init__(self, n_s, n_t):
        self.n_s = n_s
        self.n_t = n_t
    
    def parse(self, outfile, soc_key, soc_scal, molsoc_input, geom_xyz, dir_para_basis):
        """
        Parse the output from DFTB+.
        
        This section has not yet been updated...
        """
        ###reading MO_energy, # of mo, occ, virt
        print "reading MO_energy, # of mo, occ, virt..." 
        MO_energies, nbov = [], []
        var = []
        with open(outfile[0], 'r') as f:
            for nline, line in enumerate(f):
                #print nline
                if len(line.split()) > 1:
                    var = line.split()[0]
                    #print "var_temp", var
                    if self.NUMBER_SEARCH.match(var):
                        MO_energies.append(var)
                        if float(line.split()[1]) - 1.0 > 1.0E-5:
                            nb = nline
                            #print nb
            nbov.extend((nline-1, nb, nline-1-nb))
        print "MO_energies", MO_energies
        print "nbov", nbov
        ###reading singlet and triplet energies
        print "reading singlet and triplet energies..."
        ene_t, ene_s = [], []
        iit, iis = 0, 0
        with open(outfile[1], 'r') as f:
            for line in f:
                if len(line.split()) > 1:
                    var = line.split()[0]
                    if self.NUMBER_SEARCH.match(var):
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
        order1 = [ene_s[i-1][0] for i in self.n_s]
        order2 = [ene_t[i-1][0] for i in self.n_t]
        order0 = order1 + order2
        singlet_energies = [ene_s[i-1][1] for i in self.n_s]
        triplet_energies = [ene_t[i-1][1] for i in self.n_t]
        print 'singlet_energies', singlet_energies
        print 'triplet_energies', triplet_energies
        print 'total order for the roots:', order0
        ###reading oversqr.dat (ao overlap)
        print "reading oversqr.dat (ao overlap)..."
        AO_overlaps = []
        with open(outfile[2], 'r') as f:
            for line in f:
                if len(line.split()) == nbov[0]:
                    AO_overlaps.extend(line.split())
        print "len(AO_overlaps)", len(AO_overlaps)
        ###reading eigenvec.out(mo coeffs)
        ###and save above basis set info  
        print "reading eigenvec.out(mo coeffs)..."
        MOA_coefficients = []
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
                        MOA_coefficients.append(a[1])
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
        CI_coefficients, ci_xpy = [], []
        ndim = nbov[1] * nbov[2]
        print "nidm", ndim
        with open(outfile[4], 'r') as f:
            for line in f:
                if(line.split()[1] not in ['S','T'] and 
                   float(line.split()[0]) != ndim): #6 data each line
                    CI_coefficients.extend(line.split())
        print "CI_coefficients", CI_coefficients[0], CI_coefficients[-1]
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
        for i in self.n_s:
            np = (nt+i-1) * ndim
            trans_dict = {trans_key[k]: CI_coefficients[np+k] for k in range(ndim)}
            for k, v in sorted(trans_dict.iteritems()):
                ci_xpy.append(v)
                #print "trans_ci_singlet:",[ord(o) for o in k.split(' ')]
        for i in self.n_t:
            np = (i-1) * ndim
            trans_dict = {trans_key[k]: CI_coefficients[np+k] for k in range(ndim)}
            for k, v in sorted(trans_dict.iteritems()):
                ci_xpy.append(v)
                #print "trans_ci_triplet:",[ord(o) for o in k.split(' ')]
        CI_coefficients = ci_xpy
        print len(ci_xpy) 
        #print ci_xpy
       

        ###prepare input for molsoc
        print "prepare input for molsoc..."
        ngeom = 1 # coordinate after line 1 in xyz file
        print geom_xyz[0]
        with open(geom_xyz[0], 'r') as f:
            for line in f:
                natom = int(self.NUMBER_SEARCH.findall(line)[0])
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


        return singlet_energies, triplet_energies, order1, order2, nbov, ao_ncart, \
              MO_energies, MOA_coefficients, AO_overlaps, CI_coefficients