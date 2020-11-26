from itertools import islice
import pkg_resources
from pathlib import Path

from pysoc.io import Output_parser

class DFTB_plus_parser(Output_parser):
    """
    Class for parsing the required output from Gaussian.
    """
    
    # Recognised orbital labels.
    ORBITALS = ['s1', 'p1', 'p2', 'p3', 'd1', 'd2', 'd3', 'd4', 'd5', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7']

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
        
        # Init attributes.
        # Orbital information.
        self.num_orbitals = 0
        self.num_occupied_orbitals = 0
        self.num_virtual_orbitals = 0

        # List of MO energies.
        self.MO_energies = []
        
        # Lists of triplet and singlet energies. Each item is an iterable where the first item is the level (1, 2, 3 etc), and the second the energy (in eV).
        self.singlet_states = []
        self.triplet_states = []
        
        # List of atomic orbital overlaps (not really sure what the format of this is).
        self.AO_overlaps = []
        
        # ao_basis is a list of len() == 2 iterables, where the first item is a shell label (S, P, SP etc), and the second is the corresponding occupancy? (1, 3, 4 etc).
        # Not sure what the purpose of ao_basis is.
        self.ao_basis = []
        
        # Alpha and beta MO coefficients.
        self.MOA_coefficients = []
        self.MOB_coefficients = []
        
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
    
    def parse_excited_states(self):
        """
        Parse excited states information.
        
        This method is called as part of parse() (you do not normally need to call this method yourself).
        """
        # TODO: This could be improved...
        triplet_level, singlet_level = 1, 1
        
        # Start reading the file.
        with open(self.excitations_file_name, 'r') as excitations_file:
            for line in excitations_file:
                
                # Split up the line on whitespace.
                parts = line.split()
                
                # Ignore comment or header lines which don't contain numbers.
                if len(parts) > 1 and self.NUMBER_SEARCH.match(parts[0]):
                    
                    # The first item (parts[0]) is the energy in eV.
                    energy = float(parts[0])
                    
                    if 'T' in line.split():
                        # Triplet excited state.
                        self.triplet_states.append([triplet_level, energy])
                        triplet_level +=1
                        
                    elif 'S' in line.split():
                        # Singlet excited state.
                        self.singlet_states.append([singlet_level, energy])
                        singlet_level += 1
                            
    def parse_MO_energies(self):
        """
        Parse MO energies.
        
        This method is called as part of parse() (you do not normally need to call this method yourself).
        """
        # Start reading.
        with open(self.band_file_name, 'r') as band_file:
            for line in band_file:
                
                # Split.
                parts = line.split()
                if len(parts) > 1 and self.NUMBER_SEARCH.match(parts[0]):
                    # The first part is the energy.
                    self.MO_energies.append(float(parts[0]))
                    
                    # The second is the occupancy.
                    occupancy = float(parts[1])
                    # For some reason, we take 0.00005 as the threshold for occupancy...
                    if occupancy - 1.0 > 1.0E-5:
                        self.num_occupied_orbitals += 1
                    else:
                        self.num_virtual_orbitals += 1
                        
                    # Increment out number of orbitals.
                    self.num_orbitals += 1
                            
    def parse_AO_overlap(self):
        """
        Parse atomic orbital overlaps.
        
        This method is called as part of parse() (you do not normally need to call this method yourself).
        """
        with open(self.overlap_matrix_file_name, 'r') as overlap_matrix_file:
            for line in overlap_matrix_file:
                if len(line.split()) == self.num_orbitals:
                    self.AO_overlaps.extend(line.split())
    
    def parse_MOA_coefficients(self):
        """
        Parse molecular orbital (alpha) coefficients.
        """
        
        max_shell = ['s1', 'p3', 'd5', 'f7']
        # Start reading.
        with open(self.eigenvectors_file_name, 'r') as eigenvectors_file:
            for line in eigenvectors_file:
                # Split on whitespace.
                parts = line.split()
                
                # Only save lines which contain coefficients.
                if len(parts) > 1 and parts[0] in self.ORBITALS:
                    self.MOA_coefficients.append(parts[1])
        
        # Start reading from the start again to determine which orbital shells we're dealing with.
        with open(self.eigenvectors_file_name, 'r') as eigenvectors_file:
            for line in eigenvectors_file:
                # Split on whitespace.
                parts = line.split()
                
                if len(parts) > 1:
                    # Save the number of shells once we've finished reading each orbital
                    if parts[0] in max_shell:
                        # For d and f orbitals we save 6 and 10 shells rather than the 5 and 7.
                        # Otherwise we just save the number of given shells.
                        if parts[0] == 'd5':
                            #self.ao_basis.append('6')
                            self.ao_basis.append(6)
                        elif parts[0] == 'f7':
                            self.ao_basis.append('10')
                            self.ao_basis.append(10)
                        else:
                            self.ao_basis.append(int(list(parts[0])[1]))
                            
                    # We only need to save shells once, so stop after the first orbital.
                    elif all(s in parts for s in ['Eigenvector:', '2']):
                        break
        
    def parse(self):
        """
        Parse required data from the specified files.
        """
        self.parse_excited_states()
        self.parse_MO_energies()
        self.parse_AO_overlap()
        self.parse_MOA_coefficients()
        
    @property
    def ao_ncart(self):
        """
        A 3 membered list consisting of:
         - The total number of basis sets.
         - The number of occupied orbitals.
         - The number of 'virtual' basis sets.
         
        Not currently clear why this is necessary...
        """
        return [self.ao_basis_sum, self.num_occupied_orbitals, self.ao_basis_sum - self.num_occupied_orbitals]
    
    
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
        
        
        with open(AO_basis_file_name, 'w') as f:
            f.write('{}  {}\n'.format(len(self.ao_basis), self.ao_basis_sum))
            for i in range(len(self.ao_basis)):
                f.write('{}  '.format(self.ao_basis[i]))
       
        ###reading XplusY.DAT(X+Y coeffs) 
        ###and reorder XplusY according SPX.DAT 
        
        CI_coefficients, ci_xpy = [], []
        
        #ndim = self.num_occupied_orbitals * self.num_virtual_orbitals
        
        with open(self.x_plus_y_file_name, 'r') as f:
            for line in f:
                if(line.split()[1] not in ['S','T'] and float(line.split()[0]) != self.ndim): #6 data each line
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
            np = (nt+i-1) * self.ndim
            trans_dict = {trans_key[k]: CI_coefficients[np+k] for k in range(self.ndim)}
            for k, v in sorted(trans_dict.items()):
                ci_xpy.append(v)
                #print "trans_ci_singlet:",[ord(o) for o in k.split(' ')]
        for i in self.requested_triplets:
            np = (i-1) * self.ndim
            trans_dict = {trans_key[k]: CI_coefficients[np+k] for k in range(self.ndim)}
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


        self.CI_coefficients = CI_coefficients
        
        # Return the filenames we wrote.
        return (inp_file_name, basis_file_name, AO_basis_file_name)
