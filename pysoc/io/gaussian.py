# Methods for reading/writing Gaussian files.

from itertools import islice
from pathlib import Path
import re
import subprocess

from pysoc.io import Output_parser

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
        self.rwf_file_name = rwf_file_name
        
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
        dumped_data =  str(subprocess.check_output([self.RWFDUMP, self.rwf_file_name, "-", code], universal_newlines = True)).split("\n")
                
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


class Gaussian_parser(Output_parser):
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
        # Save our input files.
        self.log_file_name = log_file_name
        self.rwf_file_name = rwf_file_name
        
        self.requested_singlets = requested_singlets
        self.requested_triplets = requested_triplets
        
        # Init attributes.
        # Orbital information.
        self.num_orbitals = 0
        self.num_occupied_orbitals = 0
        self.num_virtual_orbitals = 0
        # Not totally clear on what this one means. Given by NFC= num in output.
        self.num_frozen_orbitals = 0
        # List of MO energies.
        self.MO_energies = []
        
        # Lists of triplet and singlet energies. Each item is an iterable where the first item is the level (1, 2, 3 etc), and the second the energy (in eV).
        self.singlet_states = []
        self.triplet_states = []
        
        # The line number (starting from 0) on which geometry info begins.
        self.geometry_start_line = None
        self.num_atoms = 0
        # The actual geometry section; atomic numbers and coords.
        # This is a list of iterables of the form [element, x_coord, y_coord, z_coord]
        self.geometry = []
        
        # Like geometry_start_line, but for atomic orbital basis sets.
        self.basis_set_start_line = None
        self.basis_set_end_line = None
        # A list of basis set information. This is cut out from a Gaussian log file without modification.
        self.basis_set = []
        
        # ao_basis is a list of len() == 2 iterables, where the first item is a shell label (S, P, SP etc), and the second is the corresponding occupancy? (1, 3, 4 etc).
        # Not sure what the purpose of ao_basis is.
        self.ao_basis = []
        
        # Alpha and beta MO coefficients.
        self.MOA_coefficients = []
        self.MOB_coefficients = []
        
        # List of atomic orbital overlaps (not really sure what the format of this is).
        self.AO_overlaps = []
        
        # List of CI coefficients (configuration interaction coefficients?)
        self.CI_coefficients = []
    
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
                
                # Look out for specific strings that contain required info.
                if all(s in line for s in ['Charge','Multiplicity']):
                    # Geometry section.
                    self.geometry_start_line = line_num + 1
                    
                elif 'NAtoms=' in line:
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
                    
                elif all(s in line for s in ['Excited State','Triplet']):
                    # Triplet excited state.
                    # Search for all numbers in the line.
                    numbers = self.NUMBER_SEARCH.findall(line)
                    
                    # Add to our list of excited states.
                    self.triplet_states.append([int(numbers[0]), float(numbers[1])])
                    
                elif all(s in line for s in ['Excited State','Singlet']):
                    # Triplet excited state.
                    # Search for all numbers in the line.
                    numbers = self.NUMBER_SEARCH.findall(line)
                    
                    # Add to our list of excited states.
                    self.singlet_states.append([int(numbers[0]), float(numbers[1])])
                    
                elif 'AO basis set' in line:
                    # Start of basis set section.
                    self.basis_set_start_line = line_num +1
                    
                elif 'primitive gaussians' in line:
                    # End of basis set section.
                    self.basis_set_end_line = line_num -1
                    
        # Check we were able to parse as many excited states as were asked.
        if max(self.requested_singlets) > len(self.singlet_states):
            raise Exception("Unable to parse requested singlet excited state num {}; only found {} singlets".format(max(self.requested_singlets), len(self.singlet_states)))
        
        if max(self.requested_triplets) > len(self.triplet_states):
            raise Exception("Unable to parse requested triplet excited state num {}; only found {} triplets".format(max(self.requested_triplets), len(self.triplet_states)))
                    
        # We open our .log file again to start reading from the top.
        with open(self.log_file_name, 'r') as log_file:
            # Cut out the geometry section.
            geometry_section = list(islice(log_file, self.geometry_start_line, self.geometry_start_line + self.num_atoms))
            
            # Split and convert to better types (don't think this is necessary either...)
            for geometry in geometry_section:
                # Split on whitespace (into atom symbol and 3 coords)
                parts = geometry.split()
                
                # Add to our geometry list.
                self.geometry.append([parts[0], float(parts[1]), float(parts[2]), float(parts[3])])    
       
        # Read basis set information.
        with open(self.log_file_name, 'r') as log_file:
            # This is an index specifying the atom we are currently working on (in self.geometry).
            atom_index = 0
            
            # Cut out the basis set section (this can be quite large?).
            for basis_set_line in list(islice(log_file, self.basis_set_start_line, self.basis_set_end_line)):
                
                # Split up the line.
                parts = basis_set_line.split()
                
                # If the line starts with 'S', 'P', 'SP' etc, then it identifies a new subshell (?).
                # We keep track of all subshells for ao_basis (I think because molsoc needs it for orbital filling, could be wrong).
                if parts[0] in self.SHELLS:
                    # Append the sub_orbital and its mapping to ao_basis...
                    sub_orbital = parts[0]
                    #self.ao_basis.append([sub_orbital, self.SHELLS[sub_orbital]])
                    self.ao_basis.append(self.SHELLS[sub_orbital])
                    
                if len(parts) == 2 and re.search(r'\d \d', basis_set_line):
                    #line = '{}  {}\n'.format(element[i], line.split()[1])
                    basis_set_line = '{}  {}\n'.format(self.geometry[atom_index][0], parts[1])
                    atom_index += 1
                    
                # Save the line in our attribute.
                self.basis_set.append(basis_set_line)
            
            # Remove the last line (a final '****').
            self.basis_set.pop()
                
        # Get data from the rwf file.
        self.parse_RWF()

        
    def write_AO_basis(self, file_name = "ao_basis.dat"):
        """
        Write the atomic orbital basis file.
        
        :param file_name: The file name to write to.
        """
        with open(file_name, 'w') as ao_basis_file:
            # First write the number of different subshells and also the total number of atomic orbitals.
            ao_basis_file.write('{}  {}\n'.format(len(self.ao_basis), self.ao_basis_sum))
            
            # Now write the 'number' of each subshell.
            for atomic_orbital in self.ao_basis:
                #ao_basis_file.write('{}  '.format(atomic_orbital[1]))
                ao_basis_file.write('{}  '.format(atomic_orbital))
    
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
                    
    def prepare_molsoc_input(self, keywords, soc_scale, output):
        """
        Prepare input files for molsoc.
        
        :param keywords: List of molsoc keywords (strings).
        :param soc_scale: Scaling factor for Zeff.
        :param output: Path to a directory to write molsoc files.
        :return: A tuple of paths to the files written (molsoc.inp, molsoc_basis, ao_basis.dat)
        """
        # Now write our molsoc file.
        inp_file_name = Path(output, "molsoc.inp")
        self.write_molsoc_input(keywords, soc_scale, inp_file_name)
        
        # Write molsoc basis set file.
        basis_file_name = Path(output, "molsoc_basis")
        self.write_molsoc_basis(basis_file_name)
                    
        # Now write the strange ao_basis.dat.
        AO_basis_file_name = Path(output, "ao_basis.dat")
        self.write_AO_basis(AO_basis_file_name)
        
        return (inp_file_name, basis_file_name, AO_basis_file_name)
    
    @property
    def num_AO_overlaps(self):
        """
        The number of atomic orbital overlaps.
        """
        return int(self.num_orbitals * (self.num_orbitals +1) /2)
    
        
    
    