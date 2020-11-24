import re

class Output_parser(object):
    """
    Abstract class for output parsers.
    """
    
    # This regex is used to extract numbers from gaussian output.
    # Don't understand the purpose of the boolean or here...
    NUMBER_SEARCH = re.compile(r'[-+]?\d*\.\d+|\d+')
    
            