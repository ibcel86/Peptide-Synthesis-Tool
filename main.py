import customtkinter as ctk
from customtkinter import filedialog
from CTkMessagebox import CTkMessagebox
from SequenceCalculator_v2 import CalculatePeptide, BuildSynthesisPlan, LoadFile, CompareSequences

class TabView(customtkinter.CTkTabview):
    def __init__(self, master, output_text):
        super().__init__(master)

        self.output_text = output_text
    
        # Tabs
        self.add("Synthesis Planner").grid_columnconfigure(0, weight=1)
        self.add("Modify Synthesis").grid_columnconfigure(0, weight=1)
        
        # Synthesis Planner Tab, Entry and Button
        self.title_synthesisplanner = ctk.CTkLabel(self.tab("Synthesis Planner"), text="Synthesis Planner")
        self.title_synthesisplanner.grid(row=0, column=0, padx=10, pady=(10, 0), sticky="w")

        self.entry = ctk.CTkEntry(self.tab("Synthesis Planner"), placeholder_text="Please enter your sequence eg T T Pra C: ")
        self.entry.grid(row=1, column=0, padx=10, pady=10, sticky="ew")
        self.entry.bind("<Return>", lambda event: self.process_sequence())
        self.submit_button = ctk.CTkButton(self.tab("Synthesis Planner"), text="Submit", command=self.process_sequence)
        self.submit_button.grid(row=2, column=0, padx=10, pady=10)

        # Modify Synthesis Tab, Entry and Button
        self.title_modifysynthesis = ctk.CTkLabel(self.tab("Modify Synthesis"), text="Modify Synthesis")
        self.title_modifysynthesis.grid(row=0, column=0, padx=10, pady=(10, 0), sticky="w")

        self.entry_modify = ctk.CTkEntry(self.tab("Modify Synthesis"), placeholder_text="Please enter your modified sequence eg T T Pra C: ")
        self.entry_modify.grid(row=1, column=0, padx=10, pady=10, sticky="ew")
        self.entry_modify.bind("<Return>", lambda event: self.process_compared_sequences())
        self.submit_button_modify = ctk.CTkButton(self.tab("Modify Synthesis"), text="Submit", command=self.process_compared_sequences)
        self.submit_button_modify.grid(row=2, column=0, padx=10, pady=10)

    def process_sequence(self):
        '''Processes user input and outputs data with required outputs messages and error messages'''
        try:
            sequence = self.entry.get()
            
            # Early validation checks
            if ' ' not in sequence:
                CTkMessagebox(title="Error", message="Check peptide sequence has spaces between letters", icon="cancel")
                return
                
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
            output_path = LoadFile.resource_path("")
            vial_plan_path = LoadFile.resource_path("vial plan.csv")
            synthesis_plan_path = LoadFile.resource_path("synthesis plan.csv")
            
            df_vial_plan.to_csv(vial_plan_path, index=False)
            df_synth_plan.to_csv(synthesis_plan_path, index=False)

            # Output success message
            self.output_text.delete("1.0", "end")
            self.output_text.insert("end", f"Your peptide contains {len(self.tokens)} amino acids\n")
            self.output_text.insert("end", f"Your peptide has a mass of: {validated_sequence_mass:.2f} g/mol\n\n")
            self.output_text.insert("end", f"Success! Your CSV files were saved in:\n{output_path}\n")
            self.output_text.insert("end", f"Files created:\n- vial plan.csv\n- synthesis plan.csv")
            
        except ValueError as e:
            CTkMessagebox(title="Error", message=str(e), icon="cancel")
        except Exception as e:
            CTkMessagebox(title="Error", message=f"An unexpected error occurred: {str(e)}", icon="cancel")

    def process_compared_sequences(self):
        '''Executes functions to compare sequences and output modified vial map and synthesis plan'''
        try:
            new_sequence = self.entry_modify.get()

            # Early validation checks
            if ' ' not in new_sequence:
                CTkMessagebox(title="Error", message="Check peptide sequence has spaces between letters", icon="cancel")
                return
                
            calc = CalculatePeptide()
            self.tokens, calc.original_tokens, invalid_amino_acids = calc.validate_user_sequence(new_sequence)
            
            if not self.tokens:
                CTkMessagebox(title="Error", message="No sequence loaded. Run validate_user_sequence() first.", icon="cancel")
                return
                
            if invalid_amino_acids:
                CTkMessagebox(title="Error", message=f"Invalid amino acid(s): {', '.join(invalid_amino_acids)} Check sequence is correct and entered as per the example", icon="cancel")
                return
                
            validated_sequence_mass = calc.calculate_sequence_mass(new_sequence)

            builder_instance = BuildSynthesisPlan(self.tokens, calc.original_tokens)
            comparison = CompareSequences(builder_instance)

            old_sequence = comparison.extract_old_sequence_from_csv()
            new_only = comparison.compare_sequences(old_sequence, self.tokens)
            df_combined = comparison.build_new_vial_map(new_only)
            comparison.tokens = self.tokens
            new_synthesis_plan = comparison.build_new_synthesis_plan(df_combined)

            # Get the correct output path
            output_path = LoadFile.resource_path("")
            
            # Save files with full path
            vial_plan_path = LoadFile.resource_path("new vial plan.csv")
            synthesis_plan_path = LoadFile.resource_path("new synthesis plan.csv")
            
            df_combined.to_csv(vial_plan_path, index=False)
            new_synthesis_plan.to_csv(synthesis_plan_path, index=False)

            # Output success message
            self.output_text.delete("1.0", "end")
            self.output_text.insert("end", f"Your peptide contains {len(self.tokens)} amino acids\n")
            self.output_text.insert("end", f"Your peptide has a mass of: {validated_sequence_mass:.2f} g/mol\n\n")
            self.output_text.insert("end", f"Success! Your CSV files were saved in:\n{output_path}\n")
            self.output_text.insert("end", f"Files created:\n- new vial plan.csv\n- new synthesis plan.csv")
            
        except FileNotFoundError as e:
            CTkMessagebox(title="Error", message=f"File not found: {str(e)}", icon="cancel")
        except ValueError as e:
            CTkMessagebox(title="Error", message=str(e), icon="cancel")
        except Exception as e:
            CTkMessagebox(title="Error", message=f"An unexpected error occurred: {str(e)}", icon="cancel")

class App(customtkinter.CTk):
    '''Activates the GUI etc'''
    def __init__(self):
        super().__init__()

        self.title("Peptide Sequence Tool")
        ctk.set_appearance_mode("dark")
        self.geometry("720x480")
        self.grid_columnconfigure(0, weight=4)
        self.grid_rowconfigure(1, weight=4)

        self.output_text = ctk.CTkTextbox(self, height=100, width=600)
        self.output_text.grid(row=3, column=0, padx=10, pady=10, sticky="ew")

        self.tabview = TabView(self, self.output_text)
        self.tabview.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")

if __name__ == "__main__":
    app = App()
    app.mainloop()