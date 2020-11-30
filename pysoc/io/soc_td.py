from pysoc.io.file import write_file, read_file
from pathlib import Path
from pysoc.io.gaussian import Gaussian_parser
from pysoc.io.dftb_plus import DFTB_plus_parser
import subprocess
from logging import getLogger
import pysoc
import math

class SOC():
    """
    Class to represent spin orbit coupling values.
    """
    
    def __init__(self, singlet_state, triplet_state, positive_one, zero, negative_one):
        """
        Constructor for SOC class.
        
        :param singlet_state: Symbol for the singlet state that this coupling is between (eg, "S(1)").
        :param triplet_state: Symbol for the triplet state that this coupling is between (eg, "T(1)").
        :param positive_one: SOC with quantum number +1.
        :param zero: SOC with quantum number 0.
        :param negative_one: SOC with quantum number -1.
        """
        self.singlet_state = singlet_state
        self.triplet_state = triplet_state
        self.positive_one = positive_one
        self.zero = zero
        self.negative_one = negative_one
        
    @classmethod
    def from_soc_td(self, soc_td_line):
        """
        Create a SOC class from the output produced by soc_td.
        
        :param soc_td_line: A line of output from the soc_out.dat file.
        """
        # Split on whitespace.
        parts = soc_td_line.split()
        
        # Our state labels are in the second part, which has the format <S*|Hso|T*,1,0,-1>
        identifier = parts[1].split("|")
        
        singlet_state = identifier[0][1:]
        triplet_state = identifier[2].split(",")[0]
        
        # Use our normal constructor.
        return self(singlet_state, triplet_state, float(parts[4]), float(parts[5]), float(parts[6]))
    
    
    @property
    def root_sum_square(self):
        """
        The root sum square of this SOC.
        """
        return math.sqrt(self.positive_one **2 + self.zero **2 + self.negative_one **2)
        
    def __float__(self):
        """
        Floatify this SOC class.
        """
        return self.root_sum_square
    
    def __str__(self):
        """
        Stringify this SOC class.
        """
        return "<{}|Hso|{}> (RSS, +1, 0, -1) (cm-1): {:10.5f} {:10.5f} {:10.5f} {:10.5f}".format(self.singlet_state, self.triplet_state, self.root_sum_square, self.positive_one, self.zero, self.negative_one)
    
    def list(self):
        """
        Listify this SOC.
        """
        return [self.singlet_state, self.triplet_state, self.root_sum_square, self.positive_one, self.zero, self.negative_one]

class Soc_td():
    """
    Class for managing the Soc_td program.
    """
    
    def __init__(self, molsoc):
        """
        Constructor for soc_td.
        
        :param molsoc: A Molsoc wrapper object which contains data used to write input files for soc_td.
        """
        self.molsoc = molsoc
        self.output = None
        self.SOC = []
    
    def write_AO_basis(self, file_name = "ao_basis.dat"):
        """
        Write the atomic orbital basis file.
        
        :param file_name: The file name to write to.
        """
        with open(file_name, 'w') as ao_basis_file:
            # First write the number of different subshells and also the total number of atomic orbitals.
            ao_basis_file.write('{}  {}\n'.format(len(self.molsoc.ao_basis), self.molsoc.ao_basis_sum))
            
            # Now write the 'number' of each subshell.
            for atomic_orbital in self.molsoc.ao_basis:
                #ao_basis_file.write('{}  '.format(atomic_orbital[1]))
                ao_basis_file.write('{}  '.format(atomic_orbital))
                
    
    def prepare(self,  keywords, include_ground, CI_coefficient_threshold):
        """
        Prepare input files for soc_td.
        
        :param output: Path to a directory where input files will be written. 
        """        
        # Now write the strange ao_basis.dat.
        self.AO_basis_file_name = Path(self.molsoc.output, "ao_basis.dat")
        self.write_AO_basis(self.AO_basis_file_name)
        
        # Write MOs and coefficients.
        write_file(self.molsoc.MO_energies, Path(self.molsoc.output, 'mo_ene.dat'), '{}  ')
        write_file(self.molsoc.MOA_coefficients, Path(self.molsoc.output, 'mo_coeff.dat'), '{}  ')
        write_file(self.molsoc.AO_overlaps, Path(self.molsoc.output, 'ao_overlap.dat'), '{}  ')
        write_file(self.molsoc.CI_coefficients, Path(self.molsoc.output, 'ci_coeff.dat'), '{}  ')
        
        soint = []
        Hso = dict(Hso_x='X COMPONENT', Hso_y='Y COMPONENT', Hso_z='Z COMPONENT')
        for i, cont0 in sorted(Hso.items()):
            Hso[i] = read_file(Path(self.molsoc.output, 'soint'), cont0, self.molsoc.num_transitions, 0)
            Hso[i] = [cont1 for k, cont1 in enumerate(Hso[i]) if (k-2)%3==0]
            soint = soint + Hso[i] 
        
        write_file(soint, Path(self.molsoc.output, 'soc_ao.dat'), '{}  ')
        
        
        ########################################################################
        s_matr = read_file(Path(self.molsoc.output, 'molsoc_overlap.dat'), 'AO_overlap', self.molsoc.num_transitions, 0)
        write_file(s_matr, Path(self.molsoc.output, 's_matr.dat'), '{}  ')
        ########################################################################
        #read dipole moment on atomic basis
        dipole_flag = ['False']
        if 'DIP' in keywords:
            dipole_flag = ['True']
            dipint = []
            dip = dict(dip_x='DIM=1', dip_y='DIM=2', dip_z='DIM=3')
            for i, cont0 in sorted(dip.items()):
                dip[i] = read_file(Path(self.molsoc.output, 'molsoc_dipole.dat'), cont0, self.molsoc.num_transitions, 0)
                #dipint = dipint + dip[i] + ['\n']
                dipint = dipint + dip[i] 
        
        write_file(dipint, Path(self.molsoc.output, 'dip_ao.dat'), '{}  ')
        ########################################################################
        #td soc
        #
        
        dat_out = dict(a0_code = ["gauss_tddft" if  type(self.molsoc) == Gaussian_parser else "tddftb"],
                       a0_dipole = dipole_flag,
                       a0_num_g = ["True" if include_ground else "False"], 
                       a0_thresh = [CI_coefficient_threshold],
                       a0_total_s = [len(self.molsoc.requested_singlets)],
                       a1_num_s = self.molsoc.requested_singlets,
                       a1_order_s = self.molsoc.singlet_levels,
                       a2_ene_s = self.molsoc.singlet_energies,
                       a0_total_t = [len(self.molsoc.requested_triplets)], 
                       b1_num_t = self.molsoc.requested_triplets,
                       b1_order_t = self.molsoc.triplet_levels,
                       b2_ener_t = self.molsoc.triplet_energies, 
                       c_num_bov = [self.molsoc.num_orbitals, self.molsoc.num_occupied_orbitals, self.molsoc.num_virtual_orbitals])
        #if QM_program == 'DFTB+':
        if type(self.molsoc) == DFTB_plus_parser:
            dat_out.update({'d_basis_ncart': self.molsoc.ao_ncart}) 
        
        with open(Path(self.molsoc.output, 'soc_td_input.dat'), 'w') as soc_td_input_file:
            for i, line in sorted(dat_out.items()):
                for k in line:
                    soc_td_input_file.write('{:<15}'.format(k))
                soc_td_input_file.write('{:>4}\n'.format(i))
                
                
    def run(self):
        """
        Run soc_td to calculate spin orbit coupling values.
        """
        # Run soc_td.
        try:
            subprocess.run(
                ("soc_td",),
                cwd = self.molsoc.output,
                check = True,
                universal_newlines = True,
                stdout = subprocess.PIPE,
                stderr = subprocess.STDOUT
            )
        except subprocess.CalledProcessError as e:
            # Error running soc_td
            getLogger(pysoc.logger_name).error("An error occurred in the soc_td subprogram. Dumping output:\n".format(e.stdout))
            raise e
        
        # Next read the soc output file and parse.
        with open(Path(self.molsoc.output, "soc_out.dat"), "r") as soc_file:
            # Parse each line.
            self.SOC = [SOC.from_soc_td(soc_td_line) for soc_td_line in soc_file]
                

    @property
    def table(self):
        """
        Get calculated SOC values in a list of list (table format).
        """
        # Get headers
        table = [["Singlet", "Triplet", "RSS (cm-1)", "+1 (cm-1)", "0 (cm-1)", "-1 (cm-1)"]]
         
        # Add data.
        table.extend([SOC.list() for SOC in self.SOC])
        
        # Done.
        return table
        
        
        
        
        