import customtkinter
from CTkMessagebox import CTkMessagebox
from SequenceCalculator_v2 import CalculatePeptide, BuildSynthesisPlan


class TabView(customtkinter.CTkTabview):
    def __init__(self, master):
        super().__init__(master)
    
        # Tabs
        self.add("Synthesis Planner").grid_columnconfigure(0, weight=1)
        self.add("Modify Synthesis").grid_columnconfigure(0, weight=1)
        
        # Synthesis Planner Tab, Entry and Button
        self.title_synthesisplanner = customtkinter.CTkLabel(self.tab("Synthesis Planner"), text="Synthesis Planner")
        self.title_synthesisplanner.grid(row=0, column=0, padx=10, pady=(10, 0), sticky="w")

        self.entry = customtkinter.CTkEntry(self.tab("Synthesis Planner"), placeholder_text="Please enter your sequence eg T T Pra C: ")
        self.entry.grid(row=0, column=0, padx=10, pady=10, sticky="w")
        self.entry.bind("<Return>", lambda event: self.process_sequence())
        self.submit_button = customtkinter.CTkButton(self.tab("Synthesis Planner"), text="Submit", command=self.process_sequence)
        self.submit_button.grid(row=2, column=0, padx=10, pady=10, sticky="w")

        # Modify Synthesis Tab, Entry and Button
        self.title_modifysynthesis = customtkinter.CTkLabel(self.tab("Modify Synthesis"), text="Modify Synthesis")
        self.title_modifysynthesis.grid(row=0, column=0, padx=10, pady=(10, 0), sticky="e")

    def process_sequence(self):
        sequence = self.entry.get()
        self.calc = CalculatePeptide()

        try:
            valid_sequence = self.calc.validate_user_sequence(sequence)
            tokens = self.calc.tokens
            original_tokens = self.calc.original_tokens
            planner = BuildSynthesisPlan(tokens, original_tokens)
            df_vial_plan, vial_map = planner.vial_rack_positions()
            df_synth_plan = planner.build_synthesis_plan(vial_map)

            df_vial_plan.to_csv("vial plan.csv", index=False)
            df_synth_plan.to_csv("synthesis plan.csv", index=False)

            CTkMessagebox(title="Success", message="CSV files saved successfully.", icon="check")
        except ValueError as e:
            CTkMessagebox(title="Error", message=str(e), icon="cancel")



class App(customtkinter.CTk):
    '''Activates the GUI etc'''
    def __init__(self):
        super().__init__()

        self.title("Peptide Sequence Tool")
        customtkinter.set_appearance_mode("dark")
        self.geometry("720x480")
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)

        self.tabview = TabView(self)
        self.tabview.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")


