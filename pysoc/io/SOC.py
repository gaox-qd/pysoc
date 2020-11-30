import tempfile
from pathlib import Path

from pysoc.io.dftb_plus import DFTB_plus_parser
from pysoc.io.gaussian import Gaussian_parser
from pysoc.io.soc_td import Soc_td

class Calculator():
    """
    Class for calculating spin-orbit coupling.
    """
    
    def __init__(self,
        calc_file,
        num_singlets = None,
        num_triplets = None,
        QM_program = None,
        **aux_files):
        """
        Main program function for PySOC controller program.
        
        :param calc_file: The main QM output file (.log for Gaussian, .xyz for DFTB+). Other required QM output files will be found automatically based on the location of this file.
        :param num_singlets: The number of singlet excited states to calculate SOC for. This should not exceed the number of singlets calculated by the QM program.
        :param num_triplets: The number of triplet excited states to calculate SOC for. This should not exceed the number of triplets calculated by the QM program.
        :param QM_program: A string identifying the QM program to interface with (currently, one of either 'Gaussian' or 'DFTB+'.
        """
        # TODO: Because of issue #1, we don't allow selecting exact singlets/triplets, we just ask how many.
        requested_singlets = list(range(1, num_singlets +1)) if num_singlets is not None else None
        requested_triplets = list(range(1, num_triplets +1)) if num_triplets is not None else None
                    
                
        # If we weren't told what QM_program to use, guess from the calc_file.
        calc_file = Path(calc_file)
        if QM_program is None:
            # Try and guess from the input file type.
            if calc_file.suffix.lower() == ".log":
                QM_program = "Gaussian"
            elif calc_file.suffix.lower() == ".xyz":
                QM_program = "DFTB+"
            else:
                raise Exception("Could not guess input program type from file '{}'; try specifying explicitly with '--program'".format(calc_file))
            
        # Remove any aux_files that are None.
        aux_files = {key:aux_files[key] for key in aux_files if aux_files[key] is not None}
    
        # Now we need to parse the output from our QM program.
        # Get an appropriate parser.
        if QM_program == 'Gaussian':
            # Get our calculation parser.
            self.molsoc = Gaussian_parser.from_output_files(calc_file, requested_singlets = requested_singlets, requested_triplets = requested_triplets, **aux_files)
            
            # Keywords for molsoc
            self.keywords = ('ANG', 'Zeff', 'DIP')
        
        elif QM_program == 'DFTB+':
            # Get our calculation parser.
            self.molsoc = DFTB_plus_parser.from_output_files(calc_file, requested_singlets = requested_singlets, requested_triplets = requested_triplets, **aux_files)
            
            # Keywords for molsoc
            self.keywords = ('ANG', 'Zeff', 'DIP', 'TDB')
        
        else:
            # We were given something random.
            raise Exception("Unknown or unrecognised program name '{}'".format(QM_program))
        
        # Get our soc_td object.
        self.soc_td = Soc_td(self.molsoc)
            
    def calculate(self,
        output = None,
        SOC_scale = None,
        include_ground = None,
        CI_coefficient_threshold = None,
    ):
        """
        Calculate SOC values.
        
        The calculated SOC can be accessed at self.soc_td.SOC
        
        :param output: Path to a directory where intermediate files will be written. If none is given, a temporary directory will be used (in which case these intermediate files will be unavailable to the user).
        :param SOC_scale: Scaling factor for Zeff.
        :param CI_coefficient_threshold: Threshold for CI (CIS) coefficients.
        :return: The SOC values as a table. The first row contains header information.
        """
        # Set defaults if not given.
        if include_ground is None:
            include_ground = True
            
        if CI_coefficient_threshold is None:
            CI_coefficient_threshold = 1.0e-5
        
        # First, get a temp dir if we need one.
        with tempfile.TemporaryDirectory() as tempdir:
            # Only use if necessary.
            if output is None:
                output = tempdir
                
            # Parse and prepare input for molsoc.
            self.molsoc.parse()
            self.molsoc.prepare(self.keywords, SOC_scale, output)
            
            # Run molsoc.
            self.molsoc.run()
            
            # Prepare input for soc_td.
            self.soc_td.prepare(self.keywords, include_ground, CI_coefficient_threshold)
            
            # Now call soc_td.
            self.soc_td.run()
            
            return self.soc_td.table
                