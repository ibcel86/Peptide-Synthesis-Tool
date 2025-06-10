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
    '''Calculates vial positions and rack assignments. Exports csv file that is read by the auto-reactor and auto-sampler
    to synthesise the peptide'''
    
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
        
        # Get user input for number of deprotection vials
        deprotection_start_pos = 28
        deprotection_rack = 2
        rack2_start_pos = 28
        rack2_end_pos = 54  # Assuming rack 2 has 27 positions (28-54)
        
        while True:
            try:
                num_deprotection_vials = int(input("How many deprotection vials do you want to use? "))
                if num_deprotection_vials > 0:
                    # Check if requested vials exceed available rack space
                    last_position_needed = deprotection_start_pos + num_deprotection_vials - 1
                    if last_position_needed > rack2_end_pos:
                        available_positions = rack2_end_pos - deprotection_start_pos + 1
                        print(f"WARNING: Rack 2 has positions {rack2_start_pos}-{rack2_end_pos}.")
                        print(f"Starting from position {deprotection_start_pos}, you can only use {available_positions} vials maximum.")
                        print(f"You requested {num_deprotection_vials} vials, but only {available_positions} positions are available.")
                        
                        continue_choice = input("Do you want to continue with the maximum available vials? (y/n): ").lower()
                        if continue_choice == 'y':
                            num_deprotection_vials = available_positions
                            print(f"Proceeding with {num_deprotection_vials} deprotection vials.")
                            break
                        else:
                            continue
                    else:
                        break
                else:
                    print("Please enter a positive number.")
            except ValueError:
                print("Please enter a valid number.")
        
        synthesis_rows = []
        vial_usage_counter = {}  
        step_counter = 1  # Counter for deprotection steps
        deprotection_usage_counter = 0  # Track how many times we've used deprotection reagent
        
        # Build list of deprotection vial positions
        deprotection_positions = []
        for i in range(num_deprotection_vials):
            deprotection_positions.append(deprotection_start_pos + i)  # Assuming sequential positions
        
        print(f"Using {num_deprotection_vials} deprotection vials at positions: {deprotection_positions}")

        for position, aa in enumerate(self.tokens, 1):  # Start counting from 1
            # Find all related vials for this amino acid (including split vials)
            related_vials = [v for v in vial_map.keys() if v == aa or v.startswith(aa) and v[len(aa):].isdigit()]
    
            # Sort to ensure we use vials in order (aa, aa1, aa2, etc.)
            related_vials.sort(key=lambda x: (0 if x == aa else int(x[len(aa):])))

            # Assign occurrence to first vial that still has capacity
            assigned = False
            for vial_name in related_vials:
                rack, pos, occ = vial_map[vial_name]
                used = vial_usage_counter.get(vial_name, 0)

                if used < occ:
                    vial_usage_counter[vial_name] = used + 1
                    # Add amino acid synthesis step
                    synthesis_rows.append({
                        "NAME": f"{aa}{position}",
                        "FLOW RATE A (ml/min)": 0.889,
                        "FLOW RATE B (ml/min)": 0.444,
                        "FLOW RATE D (ml/min)": 0,
                        "RESIDENCE 2": True,
                        "AUTOSAMPLER SITE A": pos,  
                        "REAGENT CONC A (M)": 0.1,
                        "AUTOSAMPLER SITE B": rack2_start_pos,
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
                    
                    # Determine which deprotection vial to use (every 10 uses, move to next vial)
                    deprotection_vial_index = deprotection_usage_counter // 10
                    if deprotection_vial_index >= len(deprotection_positions):
                        deprotection_vial_index = len(deprotection_positions) - 1  # Use last vial if we exceed
                    
                    current_deprotection_pos = deprotection_positions[deprotection_vial_index]
                    
                    # Add deprotection step immediately after
                    synthesis_rows.append({
                        "NAME": f"deprotection {step_counter}",
                        "FLOW RATE A (ml/min)": 0,
                        "FLOW RATE B (ml/min)": 0,
                        "FLOW RATE D (ml/min)": 0.8,
                        "RESIDENCE 2": True,
                        "AUTOSAMPLER SITE A": 1,
                        "REAGENT CONC A (M)": 0.1,
                        "AUTOSAMPLER SITE B": current_deprotection_pos,
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
                    
                    step_counter += 1
                    deprotection_usage_counter += 1
                    assigned = True
                    break
                    
            if not assigned:
                # fallback: no vial found or all full — unlikely if logic correct
                print(f"Warning: Could not assign vial for amino acid {aa}")
                synthesis_rows.append({
                    "NAME": f"ERROR_{aa}",
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

### Debugging print methods - delete once script works ###

calc = CalculatePeptide()
amino_acids = calc.validate_user_sequence()  # Gets tokens from user input
sequence_mass = calc.calculate_sequence_mass()

synth_plan = BuildSynthesisPlan(calc.tokens)
df_vial_plan, vial_map = synth_plan.vial_rack_positions()
df_synth_plan = synth_plan.build_synthesis_plan(vial_map)

# Export as separate csv files

df_vial_plan.to_csv("Vial_Plan.csv", index=False)
df_synth_plan.to_csv("Synthesis_Plan.csv", index=False)

