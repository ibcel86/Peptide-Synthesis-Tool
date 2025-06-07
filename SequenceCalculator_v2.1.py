import os
import sys
import pandas as pd
from math import *
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
              
        if ' ' not in user_sequence:
            raise ValueError(f"Check peptide sequence has spaces between letters")
        elif invalid_amino_acids:
            raise ValueError(f"Invalid amino acid(s): {', '.join(invalid_amino_acids)}. Check sequence is correct and entered as per the example")
        else:
            return "Your sequence is valid"

    def calculate_sequence_mass(self):
        """Calculate mass of the sequence using the loaded DataFrame"""
        find_masses = dict(zip(self.df['AA'], self.df['MW']))
        validated_sequence_mass = sum(find_masses[aa] for aa in self.tokens)
        return f'Your peptide is {len(self.tokens)} amino acids long with a mass of {validated_sequence_mass} g/mol'
    
class VialRack(CalculatePeptide):
    def __init__(self, tokens=None):
        super().__init__(tokens)
        
    def vial_rack_positions(self, conc=0.4, max_occurrence=6, max_volume=16):
        '''Finds the number of occurrences for each amino acid to find how many vials
        and racks are needed based on vial size and concentration. Builds a map of vials and assigns
        vials to a specific rack'''
        find_mws = dict(zip(self.df['AA'], self.df['MW']))
        amino_acid_occurrences = Counter(self.tokens)
        max_per_vial = floor((max_volume - 1) / 2.5)
        
        output = []
        vial_map = {} 

        rack = 1
        position = 1
        max_positions = 27

        for aa, count in amino_acid_occurrences.items():
            mw = find_mws[aa]

            if count <= max_per_vial:
                splits = [count]
            else:
                splits = []
                remaining = count
                while remaining > 0:
                    chunk = min(remaining, max_per_vial)
                    splits.append(chunk)
                    remaining -= chunk

            for i, split_count in enumerate(splits):
                name = aa if i == 0 else f"{aa}{i+1}"
                mmol_per_occurrence = (max_volume * conc) / max_occurrence
                mmol = split_count * mmol_per_occurrence
                volume = split_count * 2.5 + 1
                mass = mmol * mw / 1000

                output.append({
                    "Amino Acid": name,
                    "Rack": rack,
                    "Position": position,
                    "Occurrences": split_count,
                    "mmol": round(mmol, 2),
                    "Mass (g)": round(mass, 2),
                    "Volume (mL)": round(volume, 2)
                })

                vial_map[name] = (rack, position, split_count)

                position += 1
                if position > max_positions:
                    rack += 1
                    position = 1

        return pd.DataFrame(output), vial_map


### Debugging print methods - delete once script works ###

# Create an instance of the class (CSV is loaded once here)
vial_rack_positions = VialRack()

# Now use the instance methods
amino_acids = vial_rack_positions.validate_user_sequence()
sequence_mass = vial_rack_positions.calculate_sequence_mass()
df, vial_map = vial_rack_positions.vial_rack_positions()
print(amino_acids)
print(sequence_mass)
print(df.to_string(index=False))
