from .main import Output_parser
from .gaussian import Gaussian_parser
from .dftb_plus import DFTB_plus_parser

def parser_from_program(self, program):
    """
    Get an Output_parser class appropriate to a QM program.
    
    :param program: The name of a supported program (either 'Gaussian' or 'DFTB+'). 
    """
    if program == "Gaussian":
        return Gaussian_parser
    elif program == "DFTB+":
        return DFTB_plus_parser
    else:
        # Unrecognised.
        raise Exception("Unknown or unrecognised program name '{}'".format(program))