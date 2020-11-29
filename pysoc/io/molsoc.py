import re
from pathlib import Path
import subprocess
import pysoc
from logging import getLogger

class Molsoc(object):
    """
    Abstract class for output parsers.
    """
    
    # This regex is used to extract numbers from gaussian output.
    # Don't understand the purpose of the boolean or here...
    NUMBER_SEARCH = re.compile(r'[-+]?\d*\.\d+|\d+')
    
    MOLSOC_PATH = 'molsoc'
    
    def __init__(self,requested_singlets, requested_triplets):
        """
        """
        self.requested_singlets = requested_singlets
        self.requested_triplets = requested_triplets
        
        # Orbital information.
        self.num_orbitals = 0
        self.num_occupied_orbitals = 0
        self.num_virtual_orbitals = 0
        
        # List of MO energies.
        self.MO_energies = []
        
        # Lists of triplet and singlet energies. Each item is an iterable where the first item is the level (1, 2, 3 etc), and the second the energy (in eV).
        self.singlet_states = []
        self.triplet_states = []
        
        # This is a list of iterables of the form [element, x_coord, y_coord, z_coord]
        self.geometry = []
        
        # List of atomic orbital overlaps (not really sure what the format of this is).
        self.AO_overlaps = []
        
        # List of CI coefficients (configuration interaction coefficients?)
        self.CI_coefficients = []
        
        # Alpha and beta MO coefficients.
        self.MOA_coefficients = []
        self.MOB_coefficients = []
        
        # ao_basis is a list of len() == 2 iterables, where the first item is a shell label (S, P, SP etc), and the second is the corresponding occupancy? (1, 3, 4 etc).
        # Not sure what the purpose of ao_basis is.
        self.ao_basis = []
        
        self.output = None
    
    @property
    def ndim(self):
        """
        The product of the number of occupied and virtual orbitals.
        """
        return self.num_occupied_orbitals * self.num_virtual_orbitals
    
    @property
    def num_transitions(self):
        """
        The number of transitions.
        
        This property may be misnamed (it was simply named 'nt' in the original code...)
        """
        return self.num_orbitals **2
    
    @property
    def ao_basis_sum(self):
        """
        The total number of atomic orbitals.
        """
        #return sum(self.ao_basis[k][1] for k in range(len(self.ao_basis)))
        return sum(self.ao_basis[k] for k in range(len(self.ao_basis)))
    
    @property
    def requested_singlets(self):
        """
        The singlet excited states we have been asked to calculate SOC for.
        """
        if self._requested_singlets is None:
            # Assume all.
            return list(range(1, len(self.singlet_states) +1))
        else:
            return self._requested_singlets
        
    @requested_singlets.setter
    def requested_singlets(self, value):
        """
        """
        self._requested_singlets = value
        
    @property
    def requested_triplets(self):
        """
        The triplet excited states we have been asked to calculate SOC for.
        """
        if self._requested_triplets is None:
            # Assume all.
            return list(range(1, len(self.triplet_states) +1))
        else:
            return self._requested_triplets
        
    @requested_triplets.setter
    def requested_triplets(self, value):
        """
        """
        self._requested_triplets = value
    
    @property
    def singlet_energies(self):
        """
        A list of singlet state energies parsed from the Gaussian output files.
        """
        return [self.singlet_states[i-1][1] for i in self.requested_singlets]
    
    @property
    def singlet_levels(self):
        """
        A list of singlet state levels parsed from the Gaussian output files.
        
        These levels take into account both singlets and triplets (if appropriate).
        """
        return [self.singlet_states[i-1][0] for i in self.requested_singlets]
    
    @property
    def triplet_energies(self):
        """
        A list of triplet state energies parsed from the Gaussian output files.
        """
        return [self.triplet_states[i-1][1] for i in self.requested_triplets]
    
    @property
    def triplet_levels(self):
        """
        A list of triplet state levels parsed from the Gaussian output files.
        
        These levels take into account both singlets and triplets (if appropriate).
        """
        return [self.triplet_states[i-1][0] for i in self.requested_triplets]
    
    
    def write_molsoc_basis(self, file_name = "molsoc_basis"):
        """
        Write the molsoc basis input file.
        
        :param file_name: The file name to write to.
        """
        with open(file_name, 'w') as basis_file:
            # Write each basis set line.
            for basis_set_line in self.basis_set:
                basis_file.write(basis_set_line)
            
            # Finally write a last END.
            basis_file.write("END")
            
    def write_molsoc_input(self, keywords,  soc_scale = None, file_name = "molsoc.inp",):
        """
        Write the molsoc input file.
        
        :param keywords: List of molsoc keywords (strings).
        :param soc_scale: Scaling factor for Zeff.
        :param file_name: The file name to write to.
        """
        # Use default scaling if none given.
        soc_scale = soc_scale if soc_scale is not None else 1
        
        with open(file_name, 'wt') as inp_file:
            # First, write a header/comment (probably unnecessary).
            inp_file.write('#input for soc in atomic basis\n')
            
            # Now 'keys' (key words?) for molsoc.
            # Write each key separated by whitespace.
            # Might need a final double whitespace after last keyword, not sure...
            inp_file.write("  ".join(keywords))
            inp_file.write("\n")
            
            # Write empty comment (?)
            inp_file.write("#\n")
            
            # Now write geometry.
            for element, x_coord, y_coord, z_coord in self.geometry:
                inp_file.write('{:<4}{:>15.5f}{:>15.5f}{:>15.5f}{:>7.2f}\n'.format(element, x_coord, y_coord, z_coord, soc_scale))
            
            # Finally, end with the 'End' keyword.
            inp_file.write("End\n")
            
    def prepare(self, keywords, soc_scale, output):
        """
        Prepare input files for molsoc.
        
        :param keywords: Molsoc keywords go be written the main input file.
        :param soc_scale: Scaling factor for Zeff.
        :param output: Path to a directory where input files will be written.
        """
        self.output = output
        
        # Now write our molsoc file.
        self.inp_file_name = Path(output, "molsoc.inp")
        self.write_molsoc_input(keywords, soc_scale, self.inp_file_name)
        
        # Write molsoc basis set file.
        self.basis_file_name = Path(output, "molsoc_basis")
        self.write_molsoc_basis(self.basis_file_name)
                    
    
    
    def run(self):
        """
        Run molsoc.
        
        prepare() should be called before this method.
        """
        
        # Run molsoc.
        try:
            subprocess.run(
                (self.MOLSOC_PATH, self.inp_file_name.resolve()),
                check = True,
                universal_newlines = True,
                cwd = self.inp_file_name.parent,
                stdout = subprocess.PIPE,
                stderr = subprocess.STDOUT
            )
        except subprocess.CalledProcessError as e:
            # Error running molsoc.
            getLogger(pysoc.logger_name).error("An error occurred in the molsoc subprogram. Dumping output:\n".format(e.stdout))
            raise e
            
        
        
        
        
        