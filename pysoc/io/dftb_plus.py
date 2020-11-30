import pkg_resources
from pathlib import Path

from pysoc.io import Molsoc

class DFTB_plus_parser(Molsoc):
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
        super().__init__(requested_singlets, requested_triplets)
        
        self.xyz_file_name = xyz_file_name
        self.band_file_name = band_file_name
        self.excitations_file_name = excitations_file_name
        self.overlap_matrix_file_name = overlap_matrix_file_name
        self.eigenvectors_file_name = eigenvectors_file_name
        self.x_plus_y_file_name = x_plus_y_file_name
        self.SPX_file_name = SPX_file_name
        self.fitted_basis_directory = fitted_basis_directory
        
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
    
    @property
    def num_transitions(self):
        """
        The number of transitions.
        
        This property may be misnamed (it was simply named 'nt' in the original code...)
        """
        return self.ao_ncart[0] **2
    
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
        
        This method is called as part of parse() (you do not normally need to call this method yourself).
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
    
    def parse_CI_coefficients(self):
        """
        Parse CI coefficients.
        
        This method is called as part of parse() (you do not normally need to call this method yourself).
        """
        # There are two parts to this method:
        # First we read in all available coefficients,
        # second, we prune to only retain those coefficients we're interested in.
        CI_coefficients = []
        
        # Read in all coefficients.
        with open(self.x_plus_y_file_name, 'r') as x_plus_y_file:
            for line in x_plus_y_file:
                if(line.split()[1] not in ['S','T'] and float(line.split()[0]) != self.ndim): #6 data each line
                    CI_coefficients.extend(line.split())
                    
        # Next we read in excited state transitions (?)
        # We use this to reorder our list later.
        transitions = [] #occ->virt transition order
        with open(self.SPX_file_name, 'r') as f:
            for line in f:
                if len(line.split()) == 6: #6 data each line
                    #chr: convert num to character
                    transitions.append(
                        "".join([chr(int(line.split()[i])) for i in [3, 5]])
                    )

        for i in self.requested_singlets:
            np = (len(self.triplet_states) +i -1) * self.ndim
            transition_coefficients = {transitions[k]: CI_coefficients[np+k] for k in range(self.ndim)}
            for transition, coefficient in sorted(transition_coefficients.items()):
                self.CI_coefficients.append(coefficient)
                
        for i in self.requested_triplets:
            np = (i-1) * self.ndim
            transition_coefficients = {transitions[k]: CI_coefficients[np+k] for k in range(self.ndim)}
            for transition, coefficient in sorted(transition_coefficients.items()):
                self.CI_coefficients.append(coefficient)
                
                    
    def parse_geometries(self):
        """
        Parse xyz geometry.
        """
        with open(self.xyz_file_name, 'r') as xyz_file:
            # Split into lines.
            lines = xyz_file.read().split("\n")
            
            # The first line contains the number of atoms.
            num_atoms = int(lines[0])
            
            # Remove the first two lines and any trailing lines.
            lines = lines[2:num_atoms +2]
            
            # Parse into our geometry.
            self.geometry = [(line.split()[0], float(line.split()[1]), float(line.split()[2]), float(line.split()[3])) for line in lines]
        
    
    def parse(self):
        """
        Parse required data from the specified files.
        """
        self.parse_excited_states()
        self.parse_MO_energies()
        self.parse_AO_overlap()
        self.parse_MOA_coefficients()
        self.parse_CI_coefficients()
        self.parse_geometries()
        
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
    
    def write_molsoc_basis(self, file_name = "molsoc_basis"):
        """
        Write the molsoc basis input file.
        
        :param file_name: The file name to write to.
        """
        with open(file_name, 'w') as basis_file:
            for i, geometry in enumerate(self.geometry):
                # Write basis functions for each element (we copy from a supplied directory).
                # First decide where we're copying from.
                element_basis = Path(self.fitted_basis_directory, geometry[0] + '.basis')
                
                # Copy each element basis to a single file.
                with open(element_basis, 'r') as element_basis_file:
                    lines = element_basis_file.readlines()
                    basis_file.writelines(lines)
                
                # Write out separator between basis, which changes if we're at the end of the file.
                if i != len(self.geometry)-1:
                    basis_file.write(' ****\n')
                else:
                    basis_file.write('END')
    
    