import os
import sys
import token
import pandas as pd
import openpyxl
import re
from math import *
from math import ceil
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
        
    def vial_rack_positions(self, tokens, conc=0.4, max_occurrence=6, max_volume=16, start_rack=1, start_position=1):
        '''Generates vial map, checks position of previous map when compare sequences is run'''

        amino_acid_occurrences = Counter(tokens)
        max_per_vial = floor(max_volume / 2.5)
        output = []
        vial_map = {} 

        rack = start_rack
        position = start_position
        max_positions = 27 

        for aa, count in amino_acid_occurrences.items():
            mw = self.data.mw_dict[aa]

            splits = []
            while count > 0:
                chunk = min(count, max_per_vial)
                splits.append(chunk)
                count -= chunk

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

        return pd.DataFrame(output), vial_map
    
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

    def __init__(self, builder_instance):
        self.data = DataLoader()
        self.tokens = None
        self.original_tokens = None
        self.builder = builder_instance

    def extract_old_sequence_from_csv(self):
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

    def compare_sequences(self, cleaned_tokens, new_aa):
        """Returns the amino acids in the new sequence that differ from the original.
        Each differing amino acid is marked with '*'.
        Example:
        old = [T, T, Pra, C, Q, L, I, E]
        new = [T, T, Pra, C, I, L, I, K]
        --> returns ['I*', 'K*']
        """

        differences = []

        for old, new in zip(cleaned_tokens, new_aa):
            if old != new:
                differences.append(new + "*")

        # If new_aa is longer than cleaned_tokens, marks additional amino acids
        if len(new_aa) > len(cleaned_tokens):
            add_aa = [aa + "*" for aa in new_aa[len(cleaned_tokens):]]
            differences.extend(add_aa)

        return differences
    
    def build_new_vial_map(self, new_aa):
        '''Takes in new amino acids with * marking changes, recalculates occurrence,
        then generates new vial plan with only the changed amino acids appended to the end.
        '''

        old_vial_path = os.path.join(os.path.dirname(__file__), "vial plan.csv")
        if not os.path.exists(old_vial_path):
            raise FileNotFoundError("Vial map not found. Please ensure the file is accessible.")

        df_old = pd.read_csv(old_vial_path)

        # Find last rack and position
        last_row = df_old.loc[df_old['Rack'].idxmax()]
        last_rack = last_row['Rack']
        last_position = df_old[df_old['Rack'] == last_rack]['Position'].max()

        # Compute starting rack and position
        if last_position >= 27:
            start_rack = last_rack + 1
            start_position = 1
        else:
            start_rack = last_rack
            start_position = last_position + 1

        # Generate new vial plan
        df_new, _ = self.builder.vial_rack_positions(
            start_rack=start_rack,
            start_position=start_position
        )

        # Append and return combined vial map
        df_combined = pd.concat([df_old, df_new], ignore_index=True)
        return df_combined

    def build_new_synthesis_plan():
        pass


