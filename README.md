# The Sequence Calculator

When developing novel peptide compounds in drug discovery, scientists often manually calculate the concentration, volume, and mass of each amino acid in a sequence. While this is manageable for short peptides, longer chains (e.g., 100+ mers) make the process error-prone and time-consuming — a single miscalculation can result in weeks of lost work.

The **Sequence Calculator** is a simple but powerful tool that allows scientists to input a peptide sequence and receive an Excel (or CSV) output containing all the data needed to synthesise their peptide on an automated peptide synthesiser.

---

## Key Features

- Calculates amino acid mass, concentration, and volume
- Handles long sequences with repeated residues
- Automatically determines:
  - How many vials are needed
  - How many racks are required (rack size limit = 30)
  - Vial sampling order
- Supports custom amino acids via a user-editable `amino_acid.csv` file
- Exports machine-readable synthesis plans and vial maps

---

## Sampling & Vial Logic

Due to hardware constraints:

- Reactor coil volume = **2.5 mL**
- Vial volume = **16 mL**
- ⇒ A vial can only be sampled **6 times** before a new one is needed

For example, if the amino acid `T` appears **8 times**, the system will split it into two vials: `T` and `T2`. This ensures that the 6-sample limit per vial is not exceeded.

---

## Output Files

- **Vial Plan**: A list of vials required, their positions in the rack, and associated amino acids  
- **Synthesis Plan**: A sampling order used by the peptide synthesiser  
- Both are exported as Excel (v1) or CSV (v2)

---

## Amino Acid Format

- **Single-letter codes** and **custom 3-letter codes** are both supported
- To avoid ambiguity, sequences with 3-letter codes must include spaces:

T T V Q I Pra P R A  ✅
TTVQIPraPRA          ❌

## ⚙Internal Logic (Function Overview)

- resource_path(): Ensures `amino_acid.csv` is found when running as `.exe`
- validate_sequence(): Validates sequence against the amino acid database
- calculate_sequence_mass(): Calculates molecular weight of the sequence
- find_occurrence(): Counts occurrences of each amino acid
- calculate_aa_mass_volume(): Determines vial count, splits, and rack layout
- build_synthesis_plan(): Maps the correct vial sampling order
- process_sequence()`: Ties all logic together through a Tkinter interface

---

## Executable Notes

- Use `pyinstaller` to build the `.exe`
- The `.exe` is ~33 MB; GitHub has a 25 MB limit, so only the Python code is available here
- Ensure the `.exe` and `amino_acid.csv` are placed in the same directory

---

## Version 2: Improved Output + Sequence Comparison

Version 2 introduces significant upgrades:

- **Output Format**: CSV output now adheres to the format required by legacy synthesis machines
- **Compatibility Mode**: Output was reverse-engineered based on exported experiments from the equipment itself
- **Sequence Comparison Tool** 
  - Allows uploading of a previously synthesised sequence
  - Identifies new or changed amino acids
  - Appends new vial positions to the end of existing racks
  - Automates what was previously a manual (and time-consuming) task

---

## Development Status

This project is currently complete! Updates will include bug fixes.

---
## Installation

If you want to create a working .exe, follow these instructions:

1. Place main.py and sequence_calculator.py (or the relevant updated files) into the same directory.
2. Open your terminal or command prompt and navigate to that direction:
    cd path_to_your_directory
3. Make sure you have python installed and added to your system PATH
4. Install pyinstaller if you havent already:
    pip install pyinstaller in your command prompt or terminal
5. Run pyinstaller on main.py script to generate the .exe:

    pyinstaller --onefile main.py
    The --onefile option bundles everything into a single .exe file
   
6. After pyinstaller finishes, find the .exe in the dist folder inside your directory
7. Download amino_acids.csv and keep it in the same direcotry as the .exe file. Not doing so will result in an error. The amino_acid.csv can be edited to add further amino acids 
