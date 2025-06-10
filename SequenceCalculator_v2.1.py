import os
import sys
import token
import pandas as pd
import openpyxl
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

class DataLoader:
    '''Shared data access for amino acid information'''
    def __init__(self):
        self.df = pd.read_csv(path)
        self.valid_amino_acids = set(self.df['AA'].str.strip())
        self.mw_dict = dict(zip(self.df['AA'], self.df['MW']))

class CalculatePeptide:
    '''Validates user input and calculates peptide mass'''
    
    def __init__(self):
        self.data = DataLoader()
        self.tokens = None
    
    def validate_user_sequence(self):
        '''Validates user sequence. Input gives example of how the user should input the sequence'''
        
        user_sequence = input('Please input your sequence eg: T T Pra C: ')

        # Strips whitespace and reverses sequence
        self.tokens = [aa.strip() for aa in user_sequence.split()][::-1]
            
        # Finds invalid amino acids: those that are not in the valid set
        invalid_amino_acids = [aa for aa in self.tokens if aa not in self.data.valid_amino_acids]
              
        if ' ' not in user_sequence:
            raise ValueError(f"Check peptide sequence has spaces between letters")
        elif invalid_amino_acids:
            raise ValueError(f"Invalid amino acid(s): {', '.join(invalid_amino_acids)}. Check sequence is correct and entered as per the example")
        else:
            return "Your sequence is valid"

    def calculate_sequence_mass(self):
        """Calculate mass of the sequence using the loaded DataFrame"""
        if not self.tokens:
            raise ValueError("No sequence loaded. Run validate_user_sequence() first.")
        
        validated_sequence_mass = sum(self.data.mw_dict[aa] for aa in self.tokens)
        return f'Your peptide is {len(self.tokens)} amino acids long with a mass of {validated_sequence_mass:.2f} g/mol'

class BuildSynthesisPlan():
    '''Calculates vial positions and rack assignments for amino acids'''
    
    def __init__(self, tokens):
        self.data = DataLoader()
        self.tokens = tokens
        
    def vial_rack_positions(self, conc=0.4, max_occurrence=6, max_volume=16):
        '''Finds the number of occurrences for each amino acid to find how many vials
        and racks are needed based on vial size and concentration. Builds a map of vials and assigns
        vials to a specific rack'''
        
        amino_acid_occurrences = Counter(self.tokens)
        max_per_vial = floor((max_volume - 1) / 2.5)
        
        output = []
        vial_map = {} 

        rack = 1
        position = 1
        max_positions = 27

        for aa, count in amino_acid_occurrences.items():
            mw = self.data.mw_dict[aa]

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

        df_output = pd.DataFrame(output)
        return df_output, vial_map
    
    def build_synthesis_plan(self, vial_map):
        
        synthesis_rows = []
        vial_usage_counter = {}  

        for aa in self.tokens:

            related_vials = [v for v in vial_map if v.startswith(aa)]

            # Assign occurrence to first vial that still has capacity
            assigned = False
            for vial_name in related_vials:
                rack, pos, occ = vial_map[vial_name]
                used = vial_usage_counter.get(vial_name, 0)

                if used < occ:
                    vial_usage_counter[vial_name] = used + 1
                    synthesis_rows.append({
                        "Amino Acid": aa,
                        "Vial": vial_name,
                        "Rack": rack,
                        "Position": pos,
                    })
                    assigned = True
                    break
            if not assigned:
                # fallback: no vial found or all full — unlikely if logic correct
                synthesis_rows.append({
                    "Amino Acid": aa,
                    "Vial": "UNKNOWN",
                    "Rack": None,
                    "Position": None,
                })

    
        df_synthesis_plan = pd.DataFrame(synthesis_rows)
        return df_synthesis_plan


### Debugging print methods - delete once script works ###

calc = CalculatePeptide()
amino_acids = calc.validate_user_sequence()  # Gets tokens from user input
sequence_mass = calc.calculate_sequence_mass()

synth_plan = BuildSynthesisPlan(calc.tokens)
df_vial_plan, vial_map = synth_plan.vial_rack_positions()
df_synth_plan = synth_plan.build_synthesis_plan(vial_map)

# Export to Excel

with pd.ExcelWriter("Synthesis Plan.xlsx") as writer:
    df_vial_plan.to_excel(writer, sheet_name="Vial Plan", index=False)
    df_synth_plan.to_excel(writer, sheet_name="Synthesis Plan", index=False)

print(amino_acids)
print(sequence_mass)
print(df_vial_plan.to_string(index=False))
print(df_synth_plan.to_string(index=False))
