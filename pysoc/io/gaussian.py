# Methods for reading/writing Gaussian files.

from itertools import islice
from pathlib import Path
import re
import subprocess
import periodictable

from pysoc.io import Molsoc
import cclib

'''read the output from QM calculation for the:
   a. # of basis set, virt orb, occ orb
   b. desired excitation energy 
   c. generate basis set
   d. molecular orbital, orbital energy
   e. CIS coefficients
'''

class RWF_parser():
    """
    Class for parsing data from gaussian binary files (.rwf and .chk).
    """
    
    # Path/command to Gaussian's rwfdump.
    RWFDUMP = "rwfdump"
    
    # A string that we look for in output files that indicates the real data is coming next.
    START_STRING = "Dump of file"
    
    # Some rwfdump codes.
    MO_ENERGY = '522R'
    MOA_COEFFS = '524R'
    MOB_COEFFS = '526R'
    XY_COEFFS = '635R'
    AO_OVERLAP = '514R'
    
    def __init__(self, rwf_file_name):
        """
        Constructor for RWF_parser objects.
        
        :param rwf_file_name: A file to read from. Both .rwf and .chk files are supported.
        """
        self.rwf_file_name = Path(rwf_file_name)
        
    def get_section(self, code):
        """
        Fetch a given section from a rwf file.
        
        The sections to parse are given by a code understood by rwfdump. Each code is up to 4 digits, followed by one of the following letters (taken from rwfdump.hlp):
            - I if the data is to be printed as (decimal) integers
            - H for hexadecimal
            - R for real
            - A for ascii
            
        The data is returned as a list of lines, and will not include the Gaussian preamble.
        
        :param code: A code identifying a section to extract.
        :return: A tuple of the 'Dump of file' line followed by a list of lines.
        """
        try:
            rwfdump_proc =  subprocess.run(
                [self.RWFDUMP, self.rwf_file_name.name, "-", code],
                # Capture both stdout and stderr.
                stdout = subprocess.PIPE,
                stderr = subprocess.STDOUT,
                universal_newlines = True,
                cwd = str(self.rwf_file_name.parent),
                check = True
                )
            
            dumped_data = rwfdump_proc.stdout.split("\n")
            
            #dumped_data =  str(subprocess.check_output([self.RWFDUMP, self.rwf_file_name, "-", code], universal_newlines = True)).split("\n")
        except subprocess.CalledProcessError:
            # rwfdump failed to run, check to see if the given .rwf file actually exists.
            if not self.rwf_file_name.exists():
                raise Exception("Required Gaussian .rwf file '{}' does not exist".format(self.rwf_file_name))
            else:
                raise
                
        # Find start line and remove everything before.
        try:
            start_pos = [index for index, line in enumerate(dumped_data) if self.START_STRING in line][0]
        except IndexError:
            raise Exception("Failed to parse rwfdump output from file '{}', could not find start of dump identified by '{}'".format(self.rwf_file_name, self.START_STRING))
        
        # Next, discard everything before the start line.
        dumped_data = dumped_data[start_pos:]
        
        # The next line is the header, which we'll keep.
        header = dumped_data.pop(0)
        
        # Return header and data.
        return header, dumped_data
        
    def parse(self, code, num_records = -1):
        """
        Fetch and parse given section from a rwf file.
        
        The sections to parse are given by a code understood by rwfdump. Each code is up to 4 digits, followed by one of the following letters (taken from rwfdump.hlp):
            - I if the data is to be printed as (decimal) integers
            - H for hexadecimal
            - R for real
            - A for ascii
            
        The dumped data will additionally be parsed into a number of records (a 1D list), based on num_lines and remaining_fields.
        
        :param code: A code identifying a section to extract.
        :param num_records: The number of records to return from the extracted data. A negative number (the default) will return all available records.
        :return: The processed data (as a 1D list).
        """
        # First, run RWFDUMP to get our data.
        dumped_data = self.get_section(code)[1]
                
        # Convert the remaining data to linear form.
        records = "  ".join(dumped_data).split()
        
        # Only keep the number of records we're interested in.
        if num_records > -1:
            records = records[:num_records]
        
        # And return the data.
        return records
    
    @property
    def num_XY_coefficients(self):
        """
        Parse the number of XY coefficients in the rwf file.
        """
        header = self.get_section(self.XY_COEFFS)[0]
        
        # This header has the format: Dump of file   xxxx length         xxxx (read left to right):
        #                                                                ^^^^
        return int(header.split()[5])


class Gaussian_parser(Molsoc):
    """
    Class for parsing the required output from Gaussian.
    """
    
    # A mapping of subshell labels to number of orbitals (?)
    SHELLS = {"S": 1, "P": 3, "SP": 4, "D": 6, "F": 10}
    
    def __init__(self, log_file_name, rwf_file_name, requested_singlets, requested_triplets):
        """
        Constructor for Gaussian_parser objects.
        
        :param log_file_name: Name/path to the log file to parse.
        :param rwf_file_name: Name/path to the read-write file to parse.
        :param requested_singlets: A list of singlet excited states to calculate SOC for.
        :param requested_triplets: A list of triplet excited states to calculate SOC for.
        """
        super().__init__(requested_singlets, requested_triplets)
        
        # Save our input files.
        self.log_file_name = log_file_name
        self.rwf_file_name = rwf_file_name
        
        # Init attributes.
        # Not totally clear on what this one means. Given by NFC= num in output.
        self.num_frozen_orbitals = 0
        
        # The line number (starting from 0) on which geometry info begins.
        self.geometry_start_line = None
        self.num_atoms = 0
        
        # Like geometry_start_line, but for atomic orbital basis sets.
        self.basis_set_start_line = None
        self.basis_set_end_line = None
        # A list of basis set information. This is cut out from a Gaussian log file without modification.
        self.basis_set = []
        
    
    @classmethod
    def from_output_files(self, log_file_name, *, rwf_file_name = None, **kwargs):
        """
        Create a Gaussian parser from a .log output file.
        
        This constructor is more intelligent than __init__() and will attempt to guess the location of other required output files from the location of the given log file.
        """
        # Convert to path.
        log_file_name = Path(log_file_name)
        
        # Return new object using normal constructor.
        return self(
            log_file_name = log_file_name,
            # Use the specified aux files if given, otherwise guess from the log file and standard file names.
            rwf_file_name = rwf_file_name if rwf_file_name is not None else log_file_name.with_suffix(".rwf"),
            **kwargs
        )
    
    def parse_RWF(self):
        """
        Parse output from the rwf file.
        """
        # Get our parser object.
        rwf_parser = RWF_parser(self.rwf_file_name)
        
        # Now extract each of the sections we require.
        
        # First MO energy.
        #MO_energy = rwf_parser.parse(rwf_parser.MO_ENERGY, self.num_orbitals)[self.num_frozen_orbitals:]
        self.MO_energies = rwf_parser.parse(rwf_parser.MO_ENERGY, self.num_orbitals)[self.num_frozen_orbitals:]
        
        # AO overlap.
        #num_AO = self.num_orbitals * (self.num_orbitals +1) /2
        self.AO_overlaps = rwf_parser.parse(rwf_parser.AO_OVERLAP, self.num_AO_overlaps)
        
        # Alpha molecular orbital coefficients.
        self.MOA_coefficients = rwf_parser.parse(rwf_parser.MOA_COEFFS)
        # Remove frozen orbitals.
        self.MOA_coefficients = self.MOA_coefficients[self.num_frozen_orbitals * self.num_orbitals : self.num_orbitals **2]
        
        # The highest numbered excited state we are interested in.
        max_level = max(self.singlet_levels + self.triplet_levels)
        
        # Don't understand the logic behind this maths; leaving as is for now...
        #ndim = self.num_occupied_orbitals * self.num_virtual_orbitals
        dat_lenth = rwf_parser.num_XY_coefficients
        mseek = int((dat_lenth-12) / (self.ndim*4+1))
        nline = self.ndim * 2 * (max_level +mseek) + 12
        
        # Read ci coefficients (whatever they are).
        self.CI_coefficients = rwf_parser.parse(rwf_parser.XY_COEFFS, nline)
        
        # Don't understand this bit either; leave as is...
        ci_xpy, ci_xmy = [], []
        
        for i in self.singlet_levels + self.triplet_levels:   #singlets come first, then triplets
            np = 12 + self.ndim * 2 * (i-1) 
            
            ci_xpy += self.CI_coefficients[np:np+self.ndim*2] #X+Y(xpy)
            np = 12 + self.ndim * 2 * (mseek+i-1) #X-Y(xmy)
            
            ci_xmy += self.CI_coefficients[np:np+self.ndim*2]
            
        self.CI_coefficients = ci_xpy + ci_xmy
    
    def parse(self):
        """
        Parse required data from the specified files.
        """
        # Open our log file.
        with open(self.log_file_name, 'rt') as log_file:
            
            # Iterate through each line in the file.
            for line_num, line in enumerate(log_file):
                
                if 'NAtoms=' in line:
                    # The number of atoms.
                    self.num_atoms = int(line.split()[1])
                    
                elif 'NFC=' in line:
                    # Orbital info.
                    # Split the line.
                    parts = line.split()
                    
                    # Save the total num of orbitals and also the number frozen.
                    self.num_orbitals = int(parts[1])
                    self.num_frozen_orbitals = int(parts[7])
                    
                elif 'NVA=' in line:
                    # More orbital info.
                    # Split the line.
                    parts = line.split()
                    
                    # Save num occupied and num virtual.
                    self.num_occupied_orbitals = int(parts[3])
                    self.num_virtual_orbitals = int(parts[7])
                    
                elif 'AO basis set' in line:
                    # Start of basis set section.
                    self.basis_set_start_line = line_num +1
                
                elif self.basis_set_start_line is not None and self.basis_set_end_line is None and line == "\n":
                
                #elif 'primitive gaussians' in line:
                    # End of basis set section.
                    self.basis_set_end_line = line_num -1
                    
        # We open our .log file again to start reading from the top.
        with open(self.log_file_name, 'r') as log_file:
            # Parse with cclib.
            ccdata = cclib.io.ccread(log_file)
            
            # Get our geometry.
            try:
                self.geometry = [[str(periodictable.elements[proton_num]), coord[0], coord[1], coord[2]] for proton_num, coord in zip(ccdata.atomnos, ccdata.atomcoords[-1])]
            except Exception as e:
                raise Exception("Failed to parse atom geometry") from e
             
            # Excited states.
            try:
                for es_index, es_symm in enumerate(ccdata.etsyms):
                    # Decide if this is a singlet or triplet.
                    if "Singlet" in es_symm:
                        es_list = self.singlet_states
                    elif "Triplet" in es_symm:
                        es_list = self.triplet_states
                    else:
                        # Unrecognised state.
                        continue
                    
                    # Add data to the identified list.
                    # Each item of this list a two membered list of:
                    # - The order/level of the excited state out of all excite states (not just this mult).
                    # - The energy of the excited state, remembering that etenergies is in wavenumbers.
                    es_list.append([es_index+1, round(self.wavenumbers_to_energy(ccdata.etenergies[es_index]), 4)])
                
            except Exception as e:
                raise Exception("Failed to parse excited states data") from e

        # Check we were able to parse as many excited states as were asked.
        if len(self.requested_singlets) > 0 and max(self.requested_singlets) > len(self.singlet_states):
            raise Exception("Unable to parse requested singlet excited state num {}; only found {} singlets".format(max(self.requested_singlets), len(self.singlet_states)))
        
        if len(self.requested_triplets) > 0 and max(self.requested_triplets) > len(self.triplet_states):
            raise Exception("Unable to parse requested triplet excited state num {}; only found {} triplets".format(max(self.requested_triplets), len(self.triplet_states)))
       
        # Read basis set information.
        with open(self.log_file_name, 'r') as log_file:
            # This is an index specifying the atom we are currently working on (in self.geometry).
            atom_index = 0
            
            # Cut out the basis set section (this can be quite large?).
            for basis_set_line in list(islice(log_file, self.basis_set_start_line, self.basis_set_end_line)):
                
                # Split up the line.
                parts = basis_set_line.split()
                
                if len(parts) == 0:
                    # Empty line.
                    continue
                
                # If the line starts with 'S', 'P', 'SP' etc, then it identifies a new subshell (?).
                # We keep track of all subshells for ao_basis (I think because molsoc needs it for orbital filling, could be wrong).
                if parts[0] in self.SHELLS:
                    # Append the sub_orbital and its mapping to ao_basis...
                    sub_orbital = parts[0]
                    #self.ao_basis.append([sub_orbital, self.SHELLS[sub_orbital]])
                    self.ao_basis.append(self.SHELLS[sub_orbital])
                    
                if len(parts) == 2 and re.search(r'\d \d', basis_set_line):
                    #line = '{}  {}\n'.format(element[i], line.split()[1])
                    try:
                        basis_set_line = '{}  {}\n'.format(self.geometry[atom_index][0], parts[1])
                    except  Exception:
                        pass
                    atom_index += 1
                    
                # Save the line in our attribute.
                self.basis_set.append(basis_set_line)
            
            # Remove the last line (a final '****').
            self.basis_set.pop()
                
        # Get data from the rwf file.
        self.parse_RWF()
    
    @property
    def num_AO_overlaps(self):
        """
        The number of atomic orbital overlaps.
        """
        return int(self.num_orbitals * (self.num_orbitals +1) /2)
    
        
    
    
