import os
import sys
import pandas as pd
import openpyxl
import math
from collections import Counter
import tkinter as tk
from tkinter import messagebox
from tkinter import filedialog

def resource_path(relative_path):
    """Return the path to an external file located next to the .exe"""
    return os.path.join(os.path.dirname(sys.executable), relative_path) if getattr(sys, 'frozen', False) else os.path.join(os.path.dirname(__file__), relative_path)

csv_path = resource_path("amino_acids.csv")

def validate_sequence(sequence, df):
    '''Checks the sequence is valid, strips spaces to ensure three letter amino acids are not confused with single 
    letter amino acids. If the sequence has invalid characters, it is rejected.'''
    valid_aa = set(df['AA'].str.strip().str.upper())
    tokens = [aa.strip().upper() for aa in sequence.split()]
    invalid_aa = [aa for aa in tokens if aa not in valid_aa]

    if invalid_aa:
        messagebox.showerror("Invalid Input", f"Invalid amino acid(s): {', '.join(invalid_aa)}")
        return None
    return tokens

def calculate_sequence_mass(tokens, df):
    '''Calculates the mass of the full peptide'''
    aa_to_mass = dict(zip(df['AA'].str.strip().str.upper(), df['MW']))
    total_mass = sum(aa_to_mass[aa] for aa in tokens)
    return total_mass, len(tokens)

def find_occurrence(tokens):
    '''Treats each amino acid as a token so that they can be individually used in calculations for mass,
    occurrence, mass etc'''
    return Counter(tokens)

def calculate_aa_mass_volume(tokens, df, conc=0.4, max_per_vial=6, max_volume=16):
    '''Calculates the mass and volume required for each amino acid, uses occurrence data
    to calculate number of vials and racks required'''
    aa_to_mw = dict(zip(df['AA'].str.strip().str.upper(), df['MW']))
    counts = Counter(tokens)
    mmol_per_occurrence = max_volume * conc / max_per_vial

    output = []
    vial_map = {} 

    rack = 1
    position = 1
    max_positions = 30

    for aa, count in counts.items():
        aa_upper = aa.upper()
        if aa_upper not in aa_to_mw:
            continue
        mw = aa_to_mw[aa_upper]

        if count > max_per_vial:
            split1 = count // 2
            split2 = count - split1

            for i, split_count in enumerate([split1, split2]):
                name = aa_upper if i == 0 else f"{aa_upper}{i+1}"
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

        else:
            mmol = count * mmol_per_occurrence
            volume = count * 2.5 + 1
            mass = mmol * mw / 1000

            output.append({
                "Amino Acid": aa_upper,
                "Rack": rack,
                "Position": position,
                "Occurrences": count,
                "mmol": round(mmol, 2),
                "Mass (g)": round(mass, 2),
                "Volume (mL)": round(volume, 2)
            })

            vial_map[aa_upper] = (rack, position, count)

            position += 1
            if position > max_positions:
                rack += 1
                position = 1

    return pd.DataFrame(output), vial_map


def build_synthesis_plan(tokens, vial_map, max_per_vial=6):
    
    synthesis_rows = []
    vial_usage_counter = {}  

    for aa in tokens:
        aa_upper = aa.upper()

        # Check if this amino acid has split vials
        # Count total occurrences to determine splits
        # vial names: "A", "A2", ...
        # assign to vial until capacity filled

        # Find all vials that start with aa_upper
        related_vials = [v for v in vial_map if v.startswith(aa_upper)]

        # Assign occurrence to first vial that still has capacity
        assigned = False
        for vial_name in related_vials:
            rack, pos, occ = vial_map[vial_name]
            used = vial_usage_counter.get(vial_name, 0)

            if used < occ:
                vial_usage_counter[vial_name] = used + 1
                synthesis_rows.append({
                    "Amino Acid": aa_upper,
                    "Vial": vial_name,
                    "Rack": rack,
                    "Position": pos,
                })
                assigned = True
                break
        if not assigned:
            # fallback: no vial found or all full — unlikely if logic correct
            synthesis_rows.append({
                "Amino Acid": aa_upper,
                "Vial": "UNKNOWN",
                "Rack": None,
                "Position": None,
            })

    return pd.DataFrame(synthesis_rows)

def process_sequence():
    sequence = entry.get()
    df = pd.read_csv(csv_path)
    tokens = validate_sequence(sequence, df)
    if not tokens:
        return

    total_mass, length = calculate_sequence_mass(tokens, df)
    occurrence = find_occurrence(tokens)
    result_df, vial_map = calculate_aa_mass_volume(tokens, df)

    synthesis_df = build_synthesis_plan(tokens, vial_map)

    output_text.delete("1.0", tk.END)
    output_text.insert(tk.END, f"Sequence length: {length} AA\n")
    output_text.insert(tk.END, f"Total molecular weight: {total_mass:.2f} Da\n\n")
    output_text.insert(tk.END, "Occurrences:\n")
    for aa, count in occurrence.items():
        output_text.insert(tk.END, f"{aa}: {count}\n")

    # Save to Excel with two sheets:
    output_file_path = resource_path("sequence_preparation.xlsx")

    with pd.ExcelWriter(output_file_path) as writer:
        result_df.to_excel(writer, sheet_name="Vial Info", index=False)
        synthesis_df.to_excel(writer, sheet_name="Synthesis Plan", index=False)

    output_text.insert(tk.END, f"\nExcel file '{output_file_path}' generated with synthesis plan.")
    
# --- Tkinter GUI ---
root = tk.Tk()
root.title("Peptide Sequence Tool")

frame = tk.Frame(root, padx=10, pady=10)
frame.pack()

label = tk.Label(frame, text="Enter peptide sequence (e.g. A Pra C D):")
label.pack()

entry = tk.Entry(frame, width=50)
entry.pack()

button = tk.Button(frame, text="Process Sequence", command=process_sequence)
button.pack(pady=5)

output_text = tk.Text(frame, height=20, width=60)
output_text.pack()

root.mainloop()
