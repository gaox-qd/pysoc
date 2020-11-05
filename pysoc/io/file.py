import os
def read_file(file_name, start_str, num_lines = 1, remaining_fields = 1):
    """
    Read a given file and return a subsection of the contained data.

    The data is assumed to be in a tabular form, with each record separated by whitespace.
    The start of the tabular data will be the line following the line containing the string start_str.
    num_lines of data will then be read and split()
    Finally, if remaining_fields is not 0, one more row of data will be read in and then split(), with remaining_fields items from the resulting list being appended to the total data.
    The obtained data will then be returned (as a 1D list).
    
    :param file_name: Name/path to the file to read in.
    :param start_str: A string from which reading should start from. The returned data will start from the line after start_str is first encountered.
    :param num_lines: The number of complete lines to read.
    :param remaining_fields: The number of remaining records to read after num_lines.
    :return: The processed data (as a 1D list).
    """
    # The data we will return, a 1D list of extracted fields.
    data_out = []

    # Start reading.
    with open(file_name, 'r') as file:

        # Read through our file until we encounter our start string.
        for line_num, line in enumerate(file):
            if start_str in line:
                break

        # Now read in our file line-by-line
        for line_num, line in enumerate(file):
            data_out.extend(line.split())

            # Stop if we've read enough.
            if line_num == (num_lines -1):
                break

        # Finally, read one more line if we've been asked to.
        if remaining_fields != 0:
            for line in file:
                data_out.extend(line.split()[0:remaining_fields])
                break

    # All done.
    return data_out
   

def write_file(data, file_name, style = "{}"):
    """
    Write to a specified data file.

    :param data: The data to write (a 1D list).
    :param file_name: The file to write to.
    :param style: A formatting string which describes how to write each item in data to the given file.
    """
    # Open our file.
    with open(file_name, "w") as file:
        for line in data:
            # Write each line according to the specified line_style.
            file.write(style.format(line))
            
            
def check_file_exist(file_names):
    """
    Determine whether a (number of) files exists or not.
    
    :param file_name: List of paths of the files to check. If any of the files do not exist, an exception will be raised.
    """
    for file_name in file_names:
        if not os.path.isfile(file_name):
            # TODO: Better exception type?
            raise Exception("Required file '{}' does not exist".format(file_name))



















