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