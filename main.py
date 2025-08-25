import customtkinter as ctk
from customtkinter import filedialog
from CTkMessagebox import CTkMessagebox
from peptidesequencetool import CalculatePeptide, BuildSynthesisPlan, LoadFile, CompareSequences
import os
import pandas as pd

class TabView(ctk.CTkTabview):
    def __init__(self, master, output_text):
        super().__init__(master)

        self.output_text = output_text
    
        # Tabs
        self.add("Synthesis Planner").grid_columnconfigure(0, weight=1)
        self.add("Modify Synthesis").grid_columnconfigure(0, weight=1)
        self.add("Add Amino Acid").grid_columnconfigure(0, weight=1)
        
        # Synthesis Planner Tab, Entry and Button
        self.title_synthesisplanner = ctk.CTkLabel(self.tab("Synthesis Planner"), text="Synthesis Planner")
        self.title_synthesisplanner.grid(row=0, column=0, padx=10, pady=(10, 0), sticky="w")

        self.entry = ctk.CTkEntry(self.tab("Synthesis Planner"), placeholder_text="Please enter your sequence: ")
        self.entry.grid(row=1, column=0, padx=10, pady=10, sticky="ew")
        self.entry.bind("<Return>", lambda event: self.process_sequence())
        self.submit_button = ctk.CTkButton(self.tab("Synthesis Planner"), text="Submit", command=self.process_sequence)
        self.submit_button.grid(row=2, column=0, padx=10, pady=10)

        # Modify Synthesis Tab, Entry and Button
        self.title_modifysynthesis = ctk.CTkLabel(self.tab("Modify Synthesis"), text="Modify Synthesis")
        self.title_modifysynthesis.grid(row=0, column=0, padx=10, pady=(10, 0), sticky="w")

        self.entry_modify = ctk.CTkEntry(self.tab("Modify Synthesis"), placeholder_text="Please enter your modified sequence: ")
        self.entry_modify.grid(row=1, column=0, padx=10, pady=10, sticky="ew")
        self.entry_modify.bind("<Return>", lambda event: self.process_compared_sequences())
        self.submit_button_modify = ctk.CTkButton(self.tab("Modify Synthesis"), text="Submit", command=self.process_compared_sequences)
        self.submit_button_modify.grid(row=2, column=0, padx=10, pady=10)

        # Add Amino Acid tab
        self.title_add_amino_acid = ctk.CTkLabel(self.tab("Add Amino Acid"), text="Add Amino Acid")
        self.title_add_amino_acid.grid(row=0, column=0, padx=10, pady=(10, 0), sticky="w")

        # Add Amino Acid form
        self.aa_label = ctk.CTkLabel(self.tab("Add Amino Acid"), text="Amino Acid Code:")
        self.aa_label.grid(row=1, column=0, padx=10, pady=(10, 5), sticky="w")
        
        self.entry_aa = ctk.CTkEntry(self.tab("Add Amino Acid"), placeholder_text="e.g., Pra")
        self.entry_aa.grid(row=2, column=0, padx=10, pady=(0, 10), sticky="ew")

        self.mw_label = ctk.CTkLabel(self.tab("Add Amino Acid"), text="Molecular Weight:")
        self.mw_label.grid(row=3, column=0, padx=10, pady=(0, 5), sticky="w")
        
        self.entry_mw = ctk.CTkEntry(self.tab("Add Amino Acid"), placeholder_text="e.g., 335.35")
        self.entry_mw.grid(row=4, column=0, padx=10, pady=(0, 10), sticky="ew")

        self.name_label = ctk.CTkLabel(self.tab("Add Amino Acid"), text="Full Name:")
        self.name_label.grid(row=5, column=0, padx=10, pady=(0, 5), sticky="w")
        
        self.entry_name = ctk.CTkEntry(self.tab("Add Amino Acid"), placeholder_text="e.g., Fmoc-Pra-OH; [0.40M]")
        self.entry_name.grid(row=6, column=0, padx=10, pady=(0, 10), sticky="ew")

        self.submit_button_add_aa = ctk.CTkButton(self.tab("Add Amino Acid"), text="Add Amino Acid", command=self.add_amino_acid)
        self.submit_button_add_aa.grid(row=7, column=0, padx=10, pady=10)

        # View current amino acids button
        self.view_button_aa = ctk.CTkButton(self.tab("Add Amino Acid"), text="View Current Amino Acids", command=self.view_amino_acids)
        self.view_button_aa.grid(row=8, column=0, padx=10, pady=10)

    def process_sequence(self):
        '''Processes user input and outputs data with required outputs messages and error messages'''
        try:
            sequence = self.entry.get()

            calc = CalculatePeptide()
            self.tokens, calc.original_tokens, invalid_amino_acids = calc.validate_user_sequence(sequence)
            
            if not self.tokens:
                CTkMessagebox(title="Error", message="No sequence loaded. Run validate_user_sequence() first.", icon="cancel")
                return
                
            if invalid_amino_acids:
                CTkMessagebox(title="Error", message=f"Invalid amino acid(s): {', '.join(invalid_amino_acids)} Check sequence is correct and entered as per the example", icon="cancel")
                return
                
            validated_sequence_mass = calc.calculate_sequence_mass(sequence)

            synthesis_plan = BuildSynthesisPlan(calc.tokens, calc.original_tokens)
            df_vial_plan, vial_map = synthesis_plan.vial_rack_positions(self.tokens)
            df_synth_plan = synthesis_plan.build_synthesis_plan(vial_map)

            # Get the correct output path and save files with full path
            CTkMessagebox(title="Save File", message="Click OK to save vial plan", icon="info")
            output_path = LoadFile.resource_path("")
            initial_path = LoadFile.resource_path("vial plan.csv")  # gives full path to default file
            vial_plan_path = filedialog.asksaveasfilename(
                initialdir = os.path.dirname(initial_path),  # start in folder of the default file
                initialfile = os.path.basename(initial_path),  # default filename
                title="Save Vial Plan CSV",
                defaultextension=".csv",
                filetypes=(("CSV files", "*.csv"), ("All files", "*.*"))
            )

            if not vial_plan_path:
                return
            CTkMessagebox(title="Save File", message="Click OK to save synthesis plan", icon="info")
            initial_path = LoadFile.resource_path("synthesis plan.csv")  # gives full path to default file
            synthesis_plan_path = filedialog.asksaveasfilename(
                initialdir=os.path.dirname(initial_path),  # start in folder of the default file
                initialfile=os.path.basename(initial_path),  # default filename
                title="Save Synthesis Plan CSV",
                defaultextension=".csv",
                filetypes=(("CSV files", "*.csv"), ("All files", "*.*"))
            )

            if not synthesis_plan_path:  # user canceled
                return
            
            df_vial_plan.to_csv(vial_plan_path, index=False)
            df_synth_plan.to_csv(synthesis_plan_path, index=False)

            # Output success message
            self.output_text.delete("1.0", "end")
            self.output_text.insert("end", f"Your peptide contains {len(self.tokens)} amino acids\n")
            self.output_text.insert("end", f"Your peptide has a mass of: {validated_sequence_mass:.2f} g/mol\n\n")
            self.output_text.insert("end", f"Success! Your CSV files were saved in:\n{output_path}\n")
            
        except ValueError as e:
            CTkMessagebox(title="Error", message=str(e), icon="cancel")
        except Exception as e:
            CTkMessagebox(title="Error", message=f"An unexpected error occurred: {str(e)}", icon="cancel")

    def process_compared_sequences(self):
        '''Executes functions to compare sequences and output modified vial map and synthesis plan'''
        try:
            # Get user input
            new_sequence = self.entry_modify.get()

            # Early validation
            if ' ' not in new_sequence:
                CTkMessagebox(title="Error", message="Check peptide sequence has spaces between letters", icon="cancel")
                return

            # Validate sequence
            calc = CalculatePeptide()
            self.tokens, calc.original_tokens, invalid_amino_acids = calc.validate_user_sequence(new_sequence)

            if not self.tokens:
                CTkMessagebox(title="Error", message="No sequence loaded. Run validate_user_sequence() first.", icon="cancel")
                return

            if invalid_amino_acids:
                CTkMessagebox(
                    title="Error",
                    message=f"Invalid amino acid(s): {', '.join(invalid_amino_acids)}. Check sequence is correct.",
                    icon="cancel"
                )
                return

            validated_sequence_mass = calc.calculate_sequence_mass(new_sequence)

            # Build synthesis plan instance
            builder_instance = BuildSynthesisPlan(self.tokens, calc.original_tokens)

            # Ask user to select old synthesis plan CSV
            CTkMessagebox(title="Load Synthesis Plan", message="Load prior Synthesis Plan", icon="info")
            old_synthesis_path = filedialog.askopenfilename(
                title="Select Old Synthesis Plan CSV",
                filetypes=(("CSV files", "*.csv"), ("All files", "*.*"))
            )

            if not old_synthesis_path:  # user canceled
                return

            # Ask user to select old vial plan CSV
            CTkMessagebox(title="Load Vial Plan", message="Load prior Vial Plan", icon="info")
            old_vial_path = filedialog.askopenfilename(
                title="Select Old Vial Plan CSV",
                filetypes=(("CSV files", "*.csv"), ("All files", "*.*"))
            )

            if not old_vial_path:  # user canceled
                return

            # Compare sequences using user-selected files
            comparison = CompareSequences(builder_instance, old_synthesis_path, old_vial_path)
            old_sequence = comparison.extract_old_sequence_from_csv()
            new_only = comparison.compare_sequences(old_sequence, self.tokens)
            df_combined = comparison.build_new_vial_map(new_only)
            comparison.tokens = self.tokens
            new_synthesis_plan = comparison.build_new_synthesis_plan(df_combined)

            # Save new vial plan CSV
            CTkMessagebox(title="Save File", message="Click OK to save vial plan", icon="info")
            initial_path = LoadFile.resource_path("new vial plan.csv")
            vial_plan_path = filedialog.asksaveasfilename(
                initialdir=os.path.dirname(initial_path),
                initialfile=os.path.basename(initial_path),
                title="Save New Vial Plan CSV",
                defaultextension=".csv",
                filetypes=(("CSV files", "*.csv"), ("All files", "*.*"))
            )

            if not vial_plan_path:
                return

            # Save new synthesis plan CSV
            CTkMessagebox(title="Save File", message="Click OK to save synthesis plan", icon="info")
            initial_path = LoadFile.resource_path("new synthesis plan.csv")
            synthesis_plan_path = filedialog.asksaveasfilename(
                initialdir=os.path.dirname(initial_path),
                initialfile=os.path.basename(initial_path),
                title="Save New Synthesis Plan CSV",
                defaultextension=".csv",
                filetypes=(("CSV files", "*.csv"), ("All files", "*.*"))
            )
            
            if not synthesis_plan_path:
                return

            # Write CSVs
            df_combined.to_csv(vial_plan_path, index=False)
            new_synthesis_plan.to_csv(synthesis_plan_path, index=False)

            # Output success message
            self.output_text.delete("1.0", "end")
            self.output_text.insert("end", f"Your peptide contains {len(self.tokens)} amino acids\n")
            self.output_text.insert("end", f"Your peptide has a mass of: {validated_sequence_mass:.2f} g/mol\n\n")
            self.output_text.insert("end", f"Success! Your CSV files were saved:\n")
            self.output_text.insert("end", f"- {vial_plan_path}\n- {synthesis_plan_path}")

        except FileNotFoundError as e:
            CTkMessagebox(title="Error", message=f"File not found: {str(e)}", icon="cancel")
        except ValueError as e:
            CTkMessagebox(title="Error", message=str(e), icon="cancel")
        except Exception as e:
            CTkMessagebox(title="Error", message=f"An unexpected error occurred: {str(e)}", icon="cancel")

    def add_amino_acid(self):
        """Adds a new amino acid to amino_acids.csv (AA, MW, Name)."""
        try:
            aa_code = self.entry_aa.get().strip()
            mw_text = self.entry_mw.get().strip()
            full_name = self.entry_name.get().strip()

            if not aa_code:
                CTkMessagebox(title="Error", message="Please enter an amino acid code.", icon="cancel")
                return
            try:
                mw = float(mw_text)
            except ValueError:
                CTkMessagebox(title="Error", message="Molecular weight must be a valid number.", icon="cancel")
                return
            if not full_name:
                full_name = f"Fmoc-{aa_code}-OH; [0.40M]"

            # Ensure CSV exists and has schema
            csv_path = LoadFile.ensure_csv_schema()

            # Load, normalise columns/order
            df = pd.read_csv(csv_path)
            for col in ['AA','MW','Name']:
                if col not in df.columns:
                    df[col] = pd.Series(dtype='object')
            df = df[['AA','MW','Name']]
            df['AA'] = df['AA'].astype(str).str.strip()

            # Case-insensitive match on AA
            mask = df['AA'].str.casefold() == aa_code.casefold()
            if mask.any():
                df.loc[mask, ['MW','Name']] = [mw, full_name]
                action = "updated"
            else:
                df = pd.concat([df, pd.DataFrame([{'AA': aa_code, 'MW': mw, 'Name': full_name}])],
                            ignore_index=True)
                action = "added"

            # Atomic write to avoid partial/truncated files
            tmp_path = csv_path + ".tmp"
            df.to_csv(tmp_path, index=False)
            os.replace(tmp_path, csv_path)

            # Clear fields + UI feedback
            self.entry_aa.delete(0, 'end')
            self.entry_mw.delete(0, 'end')
            self.entry_name.delete(0, 'end')

            self.output_text.delete("1.0", "end")
            self.output_text.insert("end", f"Success! Amino acid '{aa_code}' {action}.\n")
            self.output_text.insert("end", f"CSV: {csv_path}")
            CTkMessagebox(title="Success", message=f"Amino acid '{aa_code}' {action}.", icon="check")

        except PermissionError:
            CTkMessagebox(
                title="Error",
                message="Can't write to amino_acids.csv (permission denied). "
                        "Move the .exe and CSV to a writable folder (e.g., Desktop/Documents).",
                icon="cancel"
            )
        except Exception as e:
            CTkMessagebox(title="Error", message=f"Unexpected error: {e}", icon="cancel")


    def view_amino_acids(self):
        try:
            csv_path = LoadFile.ensure_csv_schema()
            df = pd.read_csv(csv_path)
            for col in ['AA','MW','Name']:
                if col not in df.columns:
                    df[col] = ""

            self.output_text.delete("1.0", "end")
            self.output_text.insert("end", "Current Amino Acids:\n" + "="*50 + "\n\n")

            for _, row in df[['AA','MW','Name']].iterrows():
                self.output_text.insert("end", f"Code: {row['AA']}\n")
                self.output_text.insert("end", f"MW: {row['MW']}\n")
                self.output_text.insert("end", f"Name: {row['Name']}\n")
                self.output_text.insert("end", "-"*30 + "\n")

            self.output_text.insert("end", f"\nTotal amino acids: {len(df)}")
            self.output_text.update_idletasks()
            self.output_text.see("end")

        except Exception as e:
            CTkMessagebox(title="Error", message=f"Error loading amino acids: {e}", icon="cancel")

class App(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("Peptide Sequence Tool")
        ctk.set_appearance_mode("dark")

        # Get native screen resolution
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()

        # Set window to full screen size
        self.geometry(f"{screen_width}x{screen_height}")

        self.grid_columnconfigure(0, weight=4)
        self.grid_rowconfigure(1, weight=4)

        self.output_text = ctk.CTkTextbox(self, height=360, width=480)
        self.output_text.grid(row=3, column=0, padx=10, pady=10, sticky="ew")

        self.tabview = TabView(self, self.output_text)
        self.tabview.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")

if __name__ == "__main__":
    app = App()
    app.mainloop()