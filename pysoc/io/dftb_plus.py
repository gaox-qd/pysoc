from itertools import islice
import pkg_resources
from pathlib import Path

from pysoc.io import Output_parser

class DFTB_plus_parser(Output_parser):
    """
    Class for parsing the required output from Gaussian.
    """

    @classmethod
    def get_default_fitted_basis(self):
        """
        Default path to our fitted basis set dir.
        """
        return Path(pkg_resources.resource_filename('pysoc', 'data/parameters/mio-1-1-fit'))
    
    def __init__(self,
                 xyz_file_name,
                 band_file_name,
                 excitations_file_name,
                 overlap_matrix_file_name,
                 eigenvectors_file_name,
                 x_plus_y_file_name,
                 SPX_file_name,
                 requested_singlets,
                 requested_triplets,
                 fitted_basis_directory = None):
        """
        Constructor for the DFTB+ parser.
        
        :param xyz_file_name: Name/path to the xyz geometry file.
        :param band_file_name: Name/path to the density of states (DOS) output file (band.out).
        :param excitations_file_name: Name/path to the excited states file (EXC.DAT).
        :param overlap_matrix_file_name: Name/path to the file containing Hamiltonian overlap matrices (oversqr.dat).
        :param eigenvectors_file_name: Name/path to the file containing the 'eigenvectors from the Hamiltonian' (eigenvec.dat).
        :param x_plus_y_file_name: Name/path to the XplusY.DAT file.
        :param SPX_file_name: Name/path to the file containing the single particle excitations (SPX.DAT).
        :param basis_parameters_directory: Name/path to the directory containing basis set parameters (used to fit STOs to GTOs).
        :param requested_singlets: The number of singlet states to calculate with PySOC (should be no greater than the number of singlet states calculated with DFTB+).
        :param requested_triplets: The number of triplet states to calculate with PySOC (should be no greater than the number of triplet states calculated with DFTB+).
        :param fitted_basis_directory: Path to the directory containing GTO fitted basis sets. If none is given, a default fitted mio-1-1 set is used.
        """
        self.xyz_file_name = xyz_file_name
        self.band_file_name = band_file_name
        self.excitations_file_name = excitations_file_name
        self.overlap_matrix_file_name = overlap_matrix_file_name
        self.eigenvectors_file_name = eigenvectors_file_name
        self.x_plus_y_file_name = x_plus_y_file_name
        self.SPX_file_name = SPX_file_name
        self.fitted_basis_directory = fitted_basis_directory
        self.requested_singlets = requested_singlets
        self.requested_triplets = requested_triplets
        
        # Lists of triplet and singlet energies. Each item is an iterable where the first item is the level (1, 2, 3 etc), and the second the energy (in eV).
        self.singlet_states = []
        self.triplet_states = []
        
    @classmethod
    def from_output_files(self,
            xyz_file_name,
            *,
            band_file_name = None,
            excitations_file_name = None,
            overlap_matrix_file_name = None,
            eigenvectors_file_name = None,
            x_plus_y_file_name = None,
            SPX_file_name = None,
            **kwargs):
        """
        Create a DFTB+ parser from an XYZ output file.
        
        This constructor is more intelligent than __init__() and will attempt to guess the location of other required output files from the location of the given xyz file.
        """
        # Convert to path.
        xyz_file_name = Path(xyz_file_name)
        
        # Return new object using normal constructor.
        return self(
            xyz_file_name =             xyz_file_name,
            # Use the specified aux files if given, otherwise guess from the xyz file and standard file names.
            band_file_name =            band_file_name              if band_file_name is not None else              xyz_file_name.with_name("band.out"),
            excitations_file_name =     excitations_file_name       if excitations_file_name is not None else       xyz_file_name.with_name("EXC.DAT"),
            overlap_matrix_file_name =  overlap_matrix_file_name    if overlap_matrix_file_name is not None else    xyz_file_name.with_name("oversqr.dat"),
            eigenvectors_file_name =    eigenvectors_file_name      if eigenvectors_file_name is not None else      xyz_file_name.with_name("eigenvec.out"),
            x_plus_y_file_name =        x_plus_y_file_name          if x_plus_y_file_name is not None else          xyz_file_name.with_name("XplusY.DAT"),
            SPX_file_name =             SPX_file_name               if SPX_file_name is not None else               xyz_file_name.with_name("SPX.DAT"),
            **kwargs
        )
        
    @classmethod
    def from_output_file(self, *args, **kwargs):
        """
        This method is an alias for from_xyz_file().
        """
        return self.from_xyz_file(*args, **kwargs)
        
    @property
    def fitted_basis_directory(self):
        """
        Path to the fitted basis set directory.
        """
        fitted_basis_directory = getattr(self, "_fitted_basis_directory")
        return fitted_basis_directory if fitted_basis_directory is not None else self.get_default_fitted_basis()
    
    @fitted_basis_directory.setter
    def fitted_basis_directory(self, value):
        """
        Path to the fitted basis set directory.
        """
        self._fitted_basis_directory = value
        
    def parse(self):
        """
        TODO: This should do something...
        """
        pass
    
    def prepare_molsoc_input(self, keywords, soc_scale, output):
        """
        Parse the output from DFTB+.
        
        This section has not yet been fully updated...
        """
        inp_file_name = Path(output, "molsoc.inp")
        basis_file_name = Path(output, "molsoc_basis")
        AO_basis_file_name = Path(output, "ao_basis.dat")
        
        # Use default scaling if none given.
        soc_scale = soc_scale if soc_scale is not None else 1
        
        ###reading MO_energy, # of mo, occ, virt
        MO_energies, nbov = [], []
        var = []
        with open(self.band_file_name, 'r') as f:
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

        ###reading singlet and triplet energies
        
        #ene_t, ene_s = [], []
        iit, iis = 0, 0
        with open(self.excitations_file_name, 'r') as f:
            for line in f:
                if len(line.split()) > 1:
                    var = line.split()[0]
                    if self.NUMBER_SEARCH.match(var):
                        if 'T' in line.split():
                            iit = iit + 1
                            a = [iit, float(var)]
                            #ene_t.append(a)
                            self.triplet_states.append([iit, float(var)])
                        elif 'S' in line.split():
                            iis = iis + 1
                            a = [iis, float(var)]
                            #ene_s.append(a)
                            self.singlet_states.append([iis, float(var)])
        
        #singlet_levels = [ene_s[i-1][0] for i in self.requested_singlets]
        #triplet_levels = [ene_t[i-1][0] for i in self.requested_triplets]
        
        #singlet_energies = [ene_s[i-1][1] for i in self.requested_singlets]
        #triplet_energies = [ene_t[i-1][1] for i in self.requested_triplets]
        
        ###reading oversqr.dat (ao overlap)
        
        AO_overlaps = []
        with open(self.overlap_matrix_file_name, 'r') as f:
            for line in f:
                if len(line.split()) == nbov[0]:
                    AO_overlaps.extend(line.split())
        
        ###reading eigenvec.out(mo coeffs)
        ###and save above basis set info  
        
        MOA_coefficients = []
        ao_basis = []
        orb = ['s1', 'p1', 'p2', 'p3', \
               'd1', 'd2', 'd3', 'd4', 'd5', \
               'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7']
        max_shell = ['s1', 'p3', 'd5', 'f7']
        with open(self.eigenvectors_file_name, 'r') as f:
            for line in f:
                a = line.split()
                if len(a) > 1:
                    if a[0] in orb:
                        MOA_coefficients.append(a[1])
        with open(self.eigenvectors_file_name, 'r') as f:
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
        
        with open(AO_basis_file_name, 'w') as f:
            f.write('{}  {}\n'.format(nc_basis, n_basis))
            for i in range(nc_basis):
                f.write('{}  '.format(ao_basis[i]))
       
        ###reading XplusY.DAT(X+Y coeffs) 
        ###and reorder XplusY according SPX.DAT 
        
        CI_coefficients, ci_xpy = [], []
        ndim = nbov[1] * nbov[2]
        
        with open(self.x_plus_y_file_name, 'r') as f:
            for line in f:
                if(line.split()[1] not in ['S','T'] and 
                   float(line.split()[0]) != ndim): #6 data each line
                    CI_coefficients.extend(line.split())
        
        trans_key = [] #occ->virt transition order
        with open(self.SPX_file_name, 'r') as f:
            for line in f:
                if len(line.split()) == 6: #6 data each line
                    #chr: convert num to character
                    trans_temp = [chr(int(line.split()[i])) for i in [3, 5]]
                    trans_key.append("".join(trans_temp))

        #ndim = nbov[1] * nbov[2]
        #print "nidm", ndim
        nt = len(self.triplet_states) #num of triplets come before singlets
        for i in self.requested_singlets:
            np = (nt+i-1) * ndim
            trans_dict = {trans_key[k]: CI_coefficients[np+k] for k in range(ndim)}
            for k, v in sorted(trans_dict.items()):
                ci_xpy.append(v)
                #print "trans_ci_singlet:",[ord(o) for o in k.split(' ')]
        for i in self.requested_triplets:
            np = (i-1) * ndim
            trans_dict = {trans_key[k]: CI_coefficients[np+k] for k in range(ndim)}
            for k, v in sorted(trans_dict.items()):
                ci_xpy.append(v)
                #print "trans_ci_triplet:",[ord(o) for o in k.split(' ')]
        CI_coefficients = ci_xpy
         
        #print ci_xpy
       

        ###prepare input for molsoc
        ngeom = 1 # coordinate after line 1 in xyz file
        with open(self.xyz_file_name, 'r') as f:
            for line in f:
                natom = int(self.NUMBER_SEARCH.findall(line)[0])
                break
            
            lines = list(islice(f, ngeom, ngeom+natom)) 
            element = [x.split()[0] for x in lines]
            geom = [list(map(float,x.split()[1:4])) for x in lines]
            #print geom

        with open(inp_file_name, 'w') as f:
            f.write('#input for soc in atomic basis\n')
            for k in keywords:
                f.write('{}  '.format(k))
            f.write('\n')
            f.write('#\n')
            for i, xyz in enumerate(geom):
                f.write('{:<4}{:>15.5f}{:>15.5f}{:>15.5f}{:>7.2f}\n'
                      .format(element[i], xyz[0], xyz[1], xyz[2], soc_scale))
            f.write('End')
               
        ###build basis set 
        with open(basis_file_name, 'w') as fw:
            for i, el in enumerate(element):
                #bas_f = self.fitted_basis_directory+'/'+el+'.basis'
                bas_f = Path(self.fitted_basis_directory, el + '.basis')
                #print "bas_f", bas_f
                with open(bas_f, 'r') as fr:
                    lines = fr.readlines()
                    fw.writelines(lines)
                if i != natom-1:
                    fw.write(' ****\n')
                else:
                    fw.write('END')


#         self.singlet_energies = singlet_energies
#         self.triplet_energies = triplet_energies
#         self.singlet_levels = singlet_levels
#         self.triplet_levels = triplet_levels
        self.num_orbitals = nbov[0]
        self.num_occupied_orbitals = nbov[1]
        self.num_virtual_orbitals = nbov[2]
        self.ao_ncart = ao_ncart
        self.MO_energies = MO_energies
        self.MOA_coefficients = MOA_coefficients
        self.AO_overlaps = AO_overlaps
        self.CI_coefficients = CI_coefficients
        
        # Return the filenames we wrote.
        return (inp_file_name, basis_file_name, AO_basis_file_name)
