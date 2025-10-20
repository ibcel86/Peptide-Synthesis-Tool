"""
main.py

CustomTkinter GUI for the Peptide Sequence Tool.

This module provides the main graphical user interface for peptide sequence
validation, vial plan generation, synthesis plan building, and amino acid
management. It integrates with the backend logic defined in
`peptidesequencetool.py`.
"""

from __future__ import annotations
import os
import pandas as pd
import customtkinter as ctk
from customtkinter import filedialog
from CTkMessagebox import CTkMessagebox

from peptidesequencetool import CalculatePeptide, BuildSynthesisPlan, LoadFile, CompareSequences


class TabView(ctk.CTkTabview):
    """Main application tab view for peptide sequence operations."""

    def __init__(self, master: ctk.CTk, output_text: ctk.CTkTextbox) -> None:
        """Initialize the TabView containing all main tool tabs.

        Args:
            master (ctk.CTk): Parent window.
            output_text (ctk.CTkTextbox): Text box for displaying messages and results.
        """
        super().__init__(master)
        self.output_text = output_text

        # Tabs
        self.add("Synthesis Planner").grid_columnconfigure(0, weight=1)
        self.add("Modify Synthesis").grid_columnconfigure(0, weight=1)
        self.add("Add Amino Acid").grid_columnconfigure(0, weight=1)

        # --- Synthesis Planner Tab ---
        self.title_synthesisplanner = ctk.CTkLabel(
            self.tab("Synthesis Planner"), text="Synthesis Planner"
        )
        self.title_synthesisplanner.grid(row=0, column=0, padx=10, pady=(10, 0), sticky="w")

        self.entry = ctk.CTkEntry(
            self.tab("Synthesis Planner"), placeholder_text="Please enter your sequence: "
        )
        self.entry.grid(row=1, column=0, padx=10, pady=10, sticky="ew")
        self.entry.bind("<Return>", lambda event: self.process_sequence())
        self.submit_button = ctk.CTkButton(
            self.tab("Synthesis Planner"), text="Submit", command=self.process_sequence
        )
        self.submit_button.grid(row=2, column=0, padx=10, pady=10)

        # --- Modify Synthesis Tab ---
        self.title_modifysynthesis = ctk.CTkLabel(
            self.tab("Modify Synthesis"), text="Modify Synthesis"
        )
        self.title_modifysynthesis.grid(row=0, column=0, padx=10, pady=(10, 0), sticky="w")

        self.entry_modify = ctk.CTkEntry(
            self.tab("Modify Synthesis"),
            placeholder_text="Please enter your modified sequence: ",
        )
        self.entry_modify.grid(row=1, column=0, padx=10, pady=10, sticky="ew")
        self.entry_modify.bind("<Return>", lambda event: self.process_compared_sequences())
        self.submit_button_modify = ctk.CTkButton(
            self.tab("Modify Synthesis"), text="Submit", command=self.process_compared_sequences
        )
        self.submit_button_modify.grid(row=2, column=0, padx=10, pady=10)

        # --- Add Amino Acid Tab ---
        self.title_add_amino_acid = ctk.CTkLabel(
            self.tab("Add Amino Acid"), text="Add Amino Acid"
        )
        self.title_add_amino_acid.grid(row=0, column=0, padx=10, pady=(10, 0), sticky="w")

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

        self.entry_name = ctk.CTkEntry(
            self.tab("Add Amino Acid"), placeholder_text="e.g., Fmoc-Pra-OH; [0.40M]"
        )
        self.entry_name.grid(row=6, column=0, padx=10, pady=(0, 10), sticky="ew")

        self.submit_button_add_aa = ctk.CTkButton(
            self.tab("Add Amino Acid"), text="Add Amino Acid", command=self.add_amino_acid
        )
        self.submit_button_add_aa.grid(row=7, column=0, padx=10, pady=10)

        self.view_button_aa = ctk.CTkButton(
            self.tab("Add Amino Acid"),
            text="View Current Amino Acids",
            command=self.view_amino_acids,
        )
        self.view_button_aa.grid(row=8, column=0, padx=10, pady=10)

    # -------------------------------------------------------------------------
    # Sequence processing and synthesis plan generation
    # -------------------------------------------------------------------------

    def process_sequence(self) -> None:
        """Validate a peptide sequence and generate vial and synthesis plans."""
        try:
            sequence = self.entry.get()
            calc = CalculatePeptide()
            tokens, calc.original_tokens, invalid_amino_acids = calc.validate_user_sequence(sequence)

            if not tokens:
                CTkMessagebox(title="Error", message="No sequence loaded.", icon="cancel")
                return
            if invalid_amino_acids:
                CTkMessagebox(
                    title="Error",
                    message=f"Invalid amino acid(s): {', '.join(invalid_amino_acids)}",
                    icon="cancel",
                )
                return

            validated_mass = calc.calculate_sequence_mass(sequence)
            synthesis_plan = BuildSynthesisPlan(calc.tokens, calc.original_tokens)
            df_vial_plan, vial_map = synthesis_plan.vial_rack_positions(tokens)
            df_synth_plan = synthesis_plan.build_synthesis_plan(vial_map)

            CTkMessagebox(title="Save File", message="Click OK to save vial plan", icon="info")
            output_path = LoadFile.resource_path("")
            initial_path = LoadFile.resource_path("vial plan.csv")
            vial_plan_path = filedialog.asksaveasfilename(
                initialdir=os.path.dirname(initial_path),
                initialfile=os.path.basename(initial_path),
                title="Save Vial Plan CSV",
                defaultextension=".csv",
                filetypes=(("CSV files", "*.csv"), ("All files", "*.*")),
            )
            if not vial_plan_path:
                return

            CTkMessagebox(title="Save File", message="Click OK to save synthesis plan", icon="info")
            initial_path = LoadFile.resource_path("synthesis plan.csv")
            synthesis_plan_path = filedialog.asksaveasfilename(
                initialdir=os.path.dirname(initial_path),
                initialfile=os.path.basename(initial_path),
                title="Save Synthesis Plan CSV",
                defaultextension=".csv",
                filetypes=(("CSV files", "*.csv"), ("All files", "*.*")),
            )
            if not synthesis_plan_path:
                return

            df_vial_plan.to_csv(vial_plan_path, index=False)
            df_synth_plan.to_csv(synthesis_plan_path, index=False)

            self.output_text.delete("1.0", "end")
            self.output_text.insert(
                "end",
                f"Your peptide contains {len(tokens)} amino acids\n"
                f"Mass: {validated_mass:.2f} g/mol\n\n"
                f"Vial and synthesis plans saved successfully.\n{output_path}\n",
            )

        except Exception as e:
            CTkMessagebox(title="Error", message=f"An error occurred: {e}", icon="cancel")

    def process_compared_sequences(self) -> None:
        """Compare modified and previous peptide sequences, updating vial/synthesis plans."""
        try:
            new_sequence = self.entry_modify.get()
            if " " not in new_sequence:
                CTkMessagebox(title="Error", message="Add spaces between amino acids.", icon="cancel")
                return

            calc = CalculatePeptide()
            tokens, calc.original_tokens, invalid_amino_acids = calc.validate_user_sequence(
                new_sequence
            )
            if not tokens or invalid_amino_acids:
                CTkMessagebox(title="Error", message="Invalid sequence.", icon="cancel")
                return

            validated_mass = calc.calculate_sequence_mass(new_sequence)
            builder_instance = BuildSynthesisPlan(tokens, calc.original_tokens)

            CTkMessagebox(title="Load", message="Load prior Synthesis Plan", icon="info")
            old_synthesis_path = filedialog.askopenfilename(
                title="Select Old Synthesis Plan CSV",
                filetypes=(("CSV files", "*.csv"), ("All files", "*.*")),
            )
            if not old_synthesis_path:
                return

            CTkMessagebox(title="Load", message="Load prior Vial Plan", icon="info")
            old_vial_path = filedialog.askopenfilename(
                title="Select Old Vial Plan CSV",
                filetypes=(("CSV files", "*.csv"), ("All files", "*.*")),
            )
            if not old_vial_path:
                return

            comparison = CompareSequences(builder_instance, old_synthesis_path, old_vial_path)
            old_sequence = comparison.extract_old_sequence_from_csv()
            new_only = comparison.compare_sequences(old_sequence, tokens)
            df_combined = comparison.build_new_vial_map(new_only)
            comparison.tokens = tokens
            new_synth_plan = comparison.build_new_synthesis_plan(df_combined)

            CTkMessagebox(title="Save", message="Save updated vial plan", icon="info")
            initial_path = LoadFile.resource_path("new vial plan.csv")
            vial_plan_path = filedialog.asksaveasfilename(
                initialdir=os.path.dirname(initial_path),
                initialfile=os.path.basename(initial_path),
                title="Save New Vial Plan CSV",
                defaultextension=".csv",
                filetypes=(("CSV files", "*.csv"), ("All files", "*.*")),
            )
            if not vial_plan_path:
                return

            CTkMessagebox(title="Save", message="Save updated synthesis plan", icon="info")
            initial_path = LoadFile.resource_path("new synthesis plan.csv")
            synthesis_plan_path = filedialog.asksaveasfilename(
                initialdir=os.path.dirname(initial_path),
                initialfile=os.path.basename(initial_path),
                title="Save New Synthesis Plan CSV",
                defaultextension=".csv",
                filetypes=(("CSV files", "*.csv"), ("All files", "*.*")),
            )
            if not synthesis_plan_path:
                return

            df_combined.to_csv(vial_plan_path, index=False)
            new_synth_plan.to_csv(synthesis_plan_path, index=False)

            self.output_text.delete("1.0", "end")
            self.output_text.insert(
                "end",
                f"Peptide contains {len(tokens)} amino acids\n"
                f"Mass: {validated_mass:.2f} g/mol\n\n"
                f"Updated plans saved:\n{vial_plan_path}\n{synthesis_plan_path}",
            )

        except Exception as e:
            CTkMessagebox(title="Error", message=f"An error occurred: {e}", icon="cancel")

    # -------------------------------------------------------------------------
    # Amino acid csv file management
    # -------------------------------------------------------------------------

    def add_amino_acid(self) -> None:
        """Add a new amino acid entry to the amino_acids.csv file."""
        try:
            aa_code = self.entry_aa.get().strip()
            mw_text = self.entry_mw.get().strip()
            full_name = self.entry_name.get().strip() or f"Fmoc-{aa_code}-OH; [0.40M]"

            if not aa_code:
                CTkMessagebox(title="Error", message="Enter an amino acid code.", icon="cancel")
                return
            try:
                mw = float(mw_text)
            except ValueError:
                CTkMessagebox(title="Error", message="Molecular weight must be numeric.", icon="cancel")
                return

            csv_path = LoadFile.ensure_csv_schema()
            df = pd.read_csv(csv_path)
            for col in ["AA", "MW", "Name"]:
                if col not in df.columns:
                    df[col] = pd.Series(dtype="object")
            df = df[["AA", "MW", "Name"]]
            df["AA"] = df["AA"].astype(str).str.strip()

            mask = df["AA"].str.casefold() == aa_code.casefold()
            if mask.any():
                df.loc[mask, ["MW", "Name"]] = [mw, full_name]
                action = "updated"
            else:
                df = pd.concat(
                    [df, pd.DataFrame([{"AA": aa_code, "MW": mw, "Name": full_name}])],
                    ignore_index=True,
                )
                action = "added"

            tmp_path = csv_path + ".tmp"
            df.to_csv(tmp_path, index=False)
            os.replace(tmp_path, csv_path)

            self.entry_aa.delete(0, "end")
            self.entry_mw.delete(0, "end")
            self.entry_name.delete(0, "end")

            self.output_text.delete("1.0", "end")
            self.output_text.insert(
                "end", f"Success: Amino acid '{aa_code}' {action}.\nFile: {csv_path}"
            )
            CTkMessagebox(title="Success", message=f"Amino acid '{aa_code}' {action}.", icon="check")

        except Exception as e:
            CTkMessagebox(title="Error", message=f"Unexpected error: {e}", icon="cancel")

    def view_amino_acids(self) -> None:
        """Display the current amino acid table in the output text box."""
        try:
            csv_path = LoadFile.ensure_csv_schema()
            df = pd.read_csv(csv_path)
            for col in ["AA", "MW", "Name"]:
                if col not in df.columns:
                    df[col] = ""

            self.output_text.delete("1.0", "end")
            self.output_text.insert("end", "Current Amino Acids:\n" + "=" * 50 + "\n\n")
            for _, row in df[["AA", "MW", "Name"]].iterrows():
                self.output_text.insert(
                    "end",
                    f"Code: {row['AA']}\nMW: {row['MW']}\nName: {row['Name']}\n{'-'*30}\n",
                )
            self.output_text.insert("end", f"\nTotal amino acids: {len(df)}")
        except Exception as e:
            CTkMessagebox(title="Error", message=f"Error loading amino acids: {e}", icon="cancel")


class App(ctk.CTk):
    """Main application window for the Peptide Sequence Tool."""

    def __init__(self) -> None:
        """Initialize the main CustomTkinter application."""
        super().__init__()
        self.title("Peptide Sequence Tool")
        ctk.set_appearance_mode("dark")

        screen_width, screen_height = self.winfo_screenwidth(), self.winfo_screenheight()
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
