import re

class Output_parser(object):
    """
    Abstract class for output parsers.
    """
    
    # This regex is used to extract numbers from gaussian output.
    # Don't understand the purpose of the boolean or here...
    NUMBER_SEARCH = re.compile(r'[-+]?\d*\.\d+|\d+')
    
    @property
    def ndim(self):
        """
        The product of the number of occupied and virtual orbitals.
        """
        return self.num_occupied_orbitals * self.num_virtual_orbitals
    
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