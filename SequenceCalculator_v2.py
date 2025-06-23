import os
import sys
import token
import pandas as pd
import openpyxl
import re
from math import *
from collections import Counter

class LoadFile:
    '''Loads csv file globally to be used in other functions and classes'''
    @classmethod
    def resource_path(cls, relative_path):
        """Return the path to an external file located next to the .exe"""
        if getattr(sys, 'frozen', False):
            return os.path.join(os.path.dirname(sys.executable), relative_path)
        return os.path.join(os.path.dirname(__file__), relative_path)
        
    @classmethod
    def get_csv_path(cls):
        return cls.resource_path("amino_acids.csv")
        
class DataLoader:
    '''Shared data access for amino acid information'''
    def __init__(self):
        self.df = pd.read_csv(LoadFile.get_csv_path())
        self.valid_amino_acids = set(self.df['AA'].str.strip())
        self.mw_dict = dict(zip(self.df['AA'], self.df['MW']))
class CalculatePeptide:
    '''Validates user input and calculates peptide mass'''
    
    def __init__(self):
        self.data = DataLoader()
        self.tokens = None
        self.original_tokens = None 
    
    def validate_user_sequence(self, sequence):
        '''Validates user sequence. Input gives example of how the user should input the sequence'''
        
        self.original_tokens = [aa.strip() for aa in sequence.split()]
        # Reverses the sequence input so synthesis plan correctly shows the order of synthesis
        self.tokens = self.original_tokens[::-1]
        # Finds invalid amino acids: those that are not in the .csv file
        invalid_amino_acids = [aa for aa in self.original_tokens if aa not in self.data.valid_amino_acids]
        return self.tokens, self.original_tokens, invalid_amino_acids

    def calculate_sequence_mass(self, sequence):
        """Calculate mass of the sequence using the loaded DataFrame"""
        if not self.tokens:
            raise ValueError("No sequence loaded. Run validate_user_sequence() first.")
        
        validated_sequence_mass = sum(self.data.mw_dict[aa] for aa in self.tokens)
        return validated_sequence_mass
       
class BuildSynthesisPlan():
    '''Calculates vial positions and rack assignments. Exports csv file that is read by the auto-reactor and auto-sampler
    to synthesise the peptide'''
    
    def __init__(self, tokens, original_tokens=None):
        self.data = DataLoader()
        self.tokens = tokens  # This obtains the reversed sequence to ensure correct synthesis order
        self.original_tokens = original_tokens or tokens  # Required original sequence order for some calculations
        
    def vial_rack_positions(self, conc=0.4, max_occurrence=6, max_volume=16):
        '''Finds the number of occurrences for each amino acid to find how many vials
        and racks are needed based on vial size and concentration. Builds a map of vials and assigns
        vials to a specific rack'''
        
        amino_acid_occurrences = Counter(self.tokens)
        max_per_vial = floor(max_volume / 2.5)
        # Vial position logic: tracks vials and racks to ensure vials and racks are in correct order of synthesis,
        # splits vials if the volume exceeds max volume per vial. Split vials are named with a number.
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
                volume = split_count * 2.5
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
    
    def calculate_deprotection_vials_needed(self, max_volume=16, inject_vol=1.5):
        '''Calculate the number of deprotection vials needed based on peptide length and max vial volume'''

        num_deprotection_steps = len(self.tokens)
        samples_per_vial = ceil(max_volume / inject_vol)
        num_vials_needed= ceil(num_deprotection_steps/samples_per_vial)
        
        return num_vials_needed
       
    def build_synthesis_plan(self, vial_map, max_deprotection_volume=16):
        '''Builds the synthesis plan using vial rack positions for amino acids and deprotection vials.
        Outputs the .csv file which is uploaded to the auto-reactor and auto-sampler so the machine knows where and
        when to pick up a vial.'''
        # Calculate required deprotection vials automatically
        num_deprotection_vials = self.calculate_deprotection_vials_needed(max_deprotection_volume)
        
        deprotection_start_pos = 28
        deprotection_rack = 2 # Deprotection vials are always in rack 2
        rack2_start_pos = 28
        rack2_end_pos = 54

        # Check if there are enough rack spaces for deprotection vials
        last_position_needed = deprotection_start_pos + num_deprotection_vials - 1
        if last_position_needed > rack2_end_pos:
            available_positions = rack2_end_pos - deprotection_start_pos + 1
            raise ValueError(
            f"ERROR: Not enough rack space for deprotection vials!\n"
            f"Rack 2 has positions {rack2_start_pos}-{rack2_end_pos}.\n"
            f"Need {num_deprotection_vials} vials but only {available_positions} positions available.\n"
            f"Insufficient rack space for required deprotection vials")

        synthesis_rows = []
        vial_usage_counter = {}
        deprotection_usage_counter = 0

        # Builds a list of deprotection vial positions
        deprotection_positions = [deprotection_start_pos + i for i in range(num_deprotection_vials)]
        uses_per_deprotection_vial = ceil(len(self.tokens) / num_deprotection_vials)

        # Position numbers based on synthesis order
        for synthesis_position, aa in enumerate(self.tokens, 1):
            
            # Gets relevant vial(s)
            related_vials = [v for v in vial_map.keys() if v == aa or (v.startswith(aa) and v[len(aa):].isdigit())]
            related_vials.sort(key=lambda x: (0 if x == aa else int(x[len(aa):])))

            deprotection_vial_index = deprotection_usage_counter // uses_per_deprotection_vial
            if deprotection_vial_index >= len(deprotection_positions):
                deprotection_vial_index = len(deprotection_positions) - 1

            current_deprotection_pos = deprotection_positions[deprotection_vial_index]

            assigned = False
            for vial_name in related_vials:
                rack, pos, occ = vial_map[vial_name]
                used = vial_usage_counter.get(vial_name, 0)

                if used < occ:
                    vial_usage_counter[vial_name] = used + 1

                    # Amino acid step - numbered by synthesis order
                    synthesis_rows.append({
                        "NAME": f"{aa}{synthesis_position}",
                        "FLOW RATE A (ml/min)": 0.889,
                        "FLOW RATE B (ml/min)": 0.444,
                        "FLOW RATE D (ml/min)": 0,
                        "RESIDENCE 2": True,
                        "AUTOSAMPLER SITE A": pos,
                        "REAGENT CONC A (M)": 0.1,
                        "AUTOSAMPLER SITE B": current_deprotection_pos,
                        "REAGENT CONC B (M)": 0.24,
                        "DO NOT FILL": False,
                        "REAGENT USE (ml)": 4,
                        "REACTOR TEMPERATURE 2 (C)": 75,
                        "REACTOR TEMPERATURE 3 (C)": 75,
                        "WHOLE PEAK": False,
                        "DO NOT COLLECT": True,
                        "CLEANING FLOW RATE (ml/min)": 2,
                        "MANUAL CLEAN (ml)": 4
                    })

                    # Get the deprotection vial position from the previous row
                    previous_deprotection_position = synthesis_rows[-1]["AUTOSAMPLER SITE B"]

                    # Deprotection step - numbered by synthesis order
                    synthesis_rows.append({
                        "NAME": f"deprotection {synthesis_position}",
                        "FLOW RATE A (ml/min)": 0,
                        "FLOW RATE B (ml/min)": 0,
                        "FLOW RATE D (ml/min)": 0.8,
                        "RESIDENCE 2": True,
                        "AUTOSAMPLER SITE A": pos,
                        "REAGENT CONC A (M)": 0.1,
                        "AUTOSAMPLER SITE B": previous_deprotection_position,
                        "REAGENT CONC B (M)": 0.1,
                        "DO NOT FILL": False,
                        "REAGENT USE (ml)": 4,
                        "REACTOR TEMPERATURE 2 (C)": 75,
                        "REACTOR TEMPERATURE 3 (C)": 75,
                        "WHOLE PEAK": False,
                        "DO NOT COLLECT": True,
                        "CLEANING FLOW RATE (ml/min)": 2,
                        "MANUAL CLEAN (ml)": 4
                    })
                    
                    deprotection_usage_counter += 1
                    assigned = True
                    break

            if not assigned:
                print(f"Warning: Could not assign vial for amino acid {aa}")
                for name in [f"ERROR_{aa}", f"ERROR_deprotection_{synthesis_position}"]:
                    synthesis_rows.append({
                        "NAME": name,
                        "FLOW RATE A (ml/min)": 0,
                        "FLOW RATE B (ml/min)": 0,
                        "FLOW RATE D (ml/min)": 0,
                        "RESIDENCE 2": False,
                        "AUTOSAMPLER SITE A": 0,
                        "REAGENT CONC A (M)": 0,
                        "AUTOSAMPLER SITE B": 0,
                        "REAGENT CONC B (M)": 0,
                        "DO NOT FILL": True,
                        "REAGENT USE (ml)": 0,
                        "REACTOR TEMPERATURE 2 (C)": 0,
                        "REACTOR TEMPERATURE 3 (C)": 0,
                        "WHOLE PEAK": False,
                        "DO NOT COLLECT": True,
                        "CLEANING FLOW RATE (ml/min)": 0,
                        "MANUAL CLEAN (ml)": 0
                    })

        df_synthesis_plan = pd.DataFrame(synthesis_rows)
        return df_synthesis_plan

class CompareSequences():
    '''Class loads old csv files so scientists can change the peptide sequence, eg: if they substitute
    amino acids, this class runs logic to compare the sequences, finds the different amino acids (eg Pra in place of K)
    and appends the new vial to the end of the rack. The synthesis plan is also modified and saved as a new csv file along
    with the vial map.'''

    def __init__(self):
        self.data = DataLoader()
        self.tokens = None
        self.original_tokens = None

    def extract_old_sequence(self):
        '''Loads the old sequence by extracting it from the synthesis plan csv'''
        old_sequence = os.path.join(os.path.dirname(__file__), "synthesis plan.csv")
        if not os.path.exists(old_sequence):
            raise FileNotFoundError("Synthesis plan not found, please ensure the file is accessible.")
        
        df = pd.read_csv(old_sequence)

        # Filter out only the rows that are amino acid additions only
        aa_rows = df[~df['NAME'].str.contains('deprotection', case=False, na=False)]

        # Extract the amino acid base names from the 'NAME' column by removing trailing digits
        cleaned_tokens = [
            re.sub(r'\d+$', '', name.strip())
            for name in aa_rows['NAME']
        ]

        # Return the original sequence by reversing the extracted sequence from the csv file
        self.original_tokens = cleaned_tokens[::-1]
        return cleaned_tokens

    def load_previous_vial_map(self):

        old_vial_map = os.path.join(os.path.dirname(__file__), "vial plan.csv")
        if not os.path.exists(old_vial_map):
            raise FileNotFoundError("Vial map not found. Please ensure the file is accessible.")
        return f"Previous vial map loaded"

    def compare_sequences(self, cleaned_tokens, modified_sequence):
        """Returns the new tokens in the modified sequence that are different
        from the original sequence. For example:
        old = [T, T, Pra, C, Q, L, I, E]
        new = [T, T, Pra, C, I, L, I, K]
        --> returns ['I', 'K']"""

        for i in range(min(len(cleaned_tokens), len(modified_sequence))):
            if cleaned_tokens[i] != modified_sequence[i]:
                break
        else:
            # If no difference was found in common range
            i = min(len(cleaned_tokens), len(modified_sequence))

        # Return the part of the modified sequence that differs
        return modified_sequence[i:]
    
    def build_new_vial_map(self):
        '''Takes in the modified sequence that extracts the novel amino acids and appends them
        to the new vial map and re-calculates occurrences'''
        pass

    def build_new_synthesis_plan(self):
        '''Builds a new synthesis plan for the modified sequence and appends the vials to the end of 
        the rack, and moves the positions of the deprotection vials if rack space spills over'''
        pass

### Debug Function Tests ### 

    @staticmethod
    def debug_compare_sequences():
        cleaned_tokens = ['T', 'T', 'Pra', 'C', 'Q', 'L', 'I', 'E']
        modified_sequence = ['T', 'T', 'Pra', 'C', 'I', 'L', 'I', 'K']

        # Find index where sequences start to differ
        for i in range(min(len(cleaned_tokens), len(modified_sequence))):
            if cleaned_tokens[i] != modified_sequence[i]:
                break
        else:
            # If no difference was found in common range
            i = min(len(cleaned_tokens), len(modified_sequence))

        # Return the part of the modified sequence that differs
        print(modified_sequence[i:])

        
CompareSequences.debug_compare_sequences()
        







