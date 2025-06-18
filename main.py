import customtkinter
from CTkMessagebox import CTkMessagebox
from SequenceCalculator_v2 import CalculatePeptide, BuildSynthesisPlan, LoadFile


class TabView(customtkinter.CTkTabview):
    def __init__(self, master, output_text):
        super().__init__(master)

        self.output_text = output_text
    
        # Tabs
        self.add("Synthesis Planner").grid_columnconfigure(0, weight=1)
        self.add("Modify Synthesis").grid_columnconfigure(0, weight=1)
        
        # Synthesis Planner Tab, Entry and Button
        self.title_synthesisplanner = customtkinter.CTkLabel(self.tab("Synthesis Planner"), text="Synthesis Planner")
        self.title_synthesisplanner.grid(row=0, column=0, padx=10, pady=(10, 0), sticky="w")

        self.entry = customtkinter.CTkEntry(self.tab("Synthesis Planner"), placeholder_text="Please enter your sequence eg T T Pra C: ")
        self.entry.grid(row=1, column=0, padx=10, pady=10, sticky="ew")
        self.entry.bind("<Return>", lambda event: self.process_sequence())
        self.submit_button = customtkinter.CTkButton(self.tab("Synthesis Planner"), text="Submit", command=self.process_sequence)
        self.submit_button.grid(row=2, column=0, padx=10, pady=10)

        # Modify Synthesis Tab, Entry and Button
        self.title_modifysynthesis = customtkinter.CTkLabel(self.tab("Modify Synthesis"), text="Modify Synthesis")
        self.title_modifysynthesis.grid(row=0, column=0, padx=10, pady=(10, 0), sticky="w")

        self.entry_modify = customtkinter.CTkEntry(self.tab("Modify Synthesis"), placeholder_text="Please enter your sequence eg T T Pra C: ")
        self.entry_modify.grid(row=1, column=0, padx=10, pady=10, sticky="ew")
        self.entry_modify.bind("<Return>", lambda event: self.process_sequence())
        self.submit_button_modify = customtkinter.CTkButton(self.tab("Modify Synthesis"), text="Submit", command=self.process_sequence)
        self.submit_button_modify.grid(row=2, column=0, padx=10, pady=10)

    def process_sequence(self):

        try:
            sequence = self.entry.get()
            calc = CalculatePeptide()
            self.tokens, calc.original_tokens, invalid_amino_acids = calc.validate_user_sequence(sequence)
            validated_sequence_mass = calc.calculate_sequence_mass(sequence)

            synthesis_plan = BuildSynthesisPlan(calc.tokens, calc.original_tokens)
            df_vial_plan, vial_map = synthesis_plan.vial_rack_positions()
            df_synth_plan = synthesis_plan.build_synthesis_plan(vial_map)

            output_path = LoadFile.resource_path("")

            df_vial_plan.to_csv("vial plan.csv", index=False)
            df_synth_plan.to_csv("synthesis plan.csv", index=False)


            # Output text
            self.output_text.delete("1.0", "end")
            self.output_text.insert("end", f"Your peptide contains {len(self.tokens)} amino acids\n")
            self.output_text.insert("end", f"Your peptide has a mass of: {validated_sequence_mass:.2f} g/mol\n\n")

            if ' ' not in sequence:
                CTkMessagebox(title="Error", message="Check peptide sequence has spaces between letters", icon="cancel")
            elif not self.tokens:
                CTkMessagebox(title="Error", message="No sequence loaded. Run validate_user_sequence() first.", icon="cancel")
            elif invalid_amino_acids:
                CTkMessagebox(title="Error", message=f"Invalid amino acid(s): {', '.join(invalid_amino_acids)} Check sequence is correct and entered as per the example"
                              , icon="cancel")
            else:
                return self.output_text.insert("end", f"Success!, your csv files saved in {output_path}")


        except ValueError as e:
            CTkMessagebox(title="Error", message=str(e), icon="cancel")
            
class App(customtkinter.CTk):
    '''Activates the GUI etc'''
    def __init__(self):
        super().__init__()

        self.title("Peptide Sequence Tool")
        customtkinter.set_appearance_mode("dark")
        self.geometry("720x480")
        self.grid_columnconfigure(0, weight=4)
        self.grid_rowconfigure(1, weight=4)

        self.output_text = customtkinter.CTkTextbox(self, height=100, width=600)
        self.output_text.grid(row=3, column=0, padx=10, pady=10, sticky="ew")

        self.tabview = TabView(self, self.output_text)
        self.tabview.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")

if __name__ == "__main__":
    app = App()
    app.mainloop()