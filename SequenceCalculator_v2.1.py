import os
import sys
import pandas as pd
from collections import Counter

class LoadFile:
    '''Loads csv file globally to be used in other functions and classes'''
    @staticmethod
    def resource_path(relative_path):
        """Return the path to an external file located next to the .exe"""
        if getattr(sys, 'frozen', False):
            return os.path.join(os.path.dirname(sys.executable), relative_path)
        return os.path.join(os.path.dirname(__file__), relative_path)
        
    @classmethod
    def get_csv_path(cls):
        return cls.resource_path("amino_acids.csv")

path = LoadFile.get_csv_path()

class CalculatePeptide:
    '''Creates a set for the amino acids for quick look up, validates user input and calculates peptide mass'''
    
    def __init__(self, tokens=None):
        """Load CSV file once when the class is instantiated"""
        self.df = pd.read_csv(path)
        self.valid_amino_acids = set(self.df['AA'].str.strip())
        self.tokens = tokens
    
    def validate_user_sequence(self):
        '''Validates user sequence. Input gives example of how the user should input the sequence'''
        user_sequence = input('Please input your sequence eg: T T Pra C: ')
            
        # Strips whitespace and reverses sequence
        self.tokens = [aa.strip() for aa in user_sequence.split()][::-1]
            
        # Finds invalid amino acids: those that are not in the valid set
        invalid_amino_acids = [aa for aa in self.tokens if aa not in self.valid_amino_acids]
            
        if invalid_amino_acids:
            return f'Invalid amino acid(s): {', '.join(invalid_amino_acids)}. Check sequence is correct and entered as per the example'
        else:
            return f'Your sequence is valid'

    def calculate_sequence_mass(self):
        """Calculate mass of the sequence using the loaded DataFrame"""
        find_masses = dict(zip(self.df['AA'], self.df['MW']))
        validated_sequence_mass = sum(find_masses[aa] for aa in self.tokens)
        return f'Your peptide is {len(self.tokens)} long with a mass of {validated_sequence_mass} g/mol'
    
class VialRack(CalculatePeptide):
    def __init__(self, tokens=None):
        super().__init__(tokens)
        
    def find_occurrences(self):
        return Counter(self.tokens)
    


### Debugging print methods - delete once script works ###

# Create an instance of the class (CSV is loaded once here)
vial_rack_positions = VialRack()

# Now use the instance methods
amino_acids = vial_rack_positions.validate_user_sequence()
sequence_mass = vial_rack_positions.calculate_sequence_mass()
occurrences = vial_rack_positions.find_occurrences()
print(amino_acids)
print(sequence_mass)
print(occurrences)
