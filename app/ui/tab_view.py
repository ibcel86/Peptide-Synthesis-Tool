from __future__ import annotations
import os
import pandas as pd
import customtkinter as ctk
from customtkinter import filedialog
from CTkMessagebox import CTkMessagebox
from app.core.sequence_processor import ValidatePeptide
from app.core.synthesis_builder import BuildSynthesisPlan
from app.core.sequence_comparator import CompareSequences
from app.io.csv_loader import LoadFile

class PeptideTabView(ctk.CTkTabview):
    """Main application tab view for peptide sequence operations."""

    def __init__(self, master: ctk.CTk, output_text: ctk.CTkTextbox) -> None:
        super().__init__(master)
        self.output_text = output_text

        self.add("Synthesis Planner").grid_columnconfigure(0, weight=1)
        self.add("Modify Synthesis").grid_columnconfigure(0, weight=1)
        self.add("Add Amino Acid").grid_columnconfigure(0, weight=1)

        self._build_synthesis_tab()
        self._build_modify_tab()
        self._build_add_amino_acid_tab()

    def _build_synthesis_tab(self) -> None:
        tab = self.tab("Synthesis Planner")

        self.title_synthesisplanner = ctk.CTkLabel(tab, text="Synthesis Planner")
        self.title_synthesisplanner.grid(row=0, column=0, padx=10, pady=(10, 0), sticky="w")

        self.entry = ctk.CTkEntry(tab, placeholder_text="Please enter your sequence:")
        self.entry.grid(row=1, column=0, padx=10, pady=10, sticky="ew")
        self.entry.bind("<Return>", lambda event: self.process_sequence())

        self.submit_button = ctk.CTkButton(
            tab, text="Submit", command=self.process_sequence
        )
        self.submit_button.grid(row=2, column=0, padx=10, pady=10)

    def _build_modify_tab(self) -> None:
        tab = self.tab("Modify Synthesis")

        self.title_modifysynthesis = ctk.CTkLabel(tab, text="Modify Synthesis")
        self.title_modifysynthesis.grid(row=0, column=0, padx=10, pady=(10, 0), sticky="w")

        self.entry_modify = ctk.CTkEntry(
            tab,
            placeholder_text="Please enter your modified sequence:",
        )
        self.entry_modify.grid(row=1, column=0, padx=10, pady=10, sticky="ew")
        self.entry_modify.bind("<Return>", lambda event: self.process_compared_sequences())

        self.submit_button_modify = ctk.CTkButton(
            tab, text="Submit", command=self.process_compared_sequences
        )
        self.submit_button_modify.grid(row=2, column=0, padx=10, pady=10)

    def _build_add_amino_acid_tab(self) -> None:
        tab = self.tab("Add Amino Acid")

        self.title_add_amino_acid = ctk.CTkLabel(tab, text="Add Amino Acid")
        self.title_add_amino_acid.grid(row=0, column=0, padx=10, pady=(10, 0), sticky="w")

        self.aa_label = ctk.CTkLabel(tab, text="Amino Acid Code:")
        self.aa_label.grid(row=1, column=0, padx=10, pady=(10, 5), sticky="w")

        self.entry_aa = ctk.CTkEntry(tab, placeholder_text="e.g., Pra")
        self.entry_aa.grid(row=2, column=0, padx=10, pady=(0, 10), sticky="ew")

        self.mw_label = ctk.CTkLabel(tab, text="Molecular Weight:")
        self.mw_label.grid(row=3, column=0, padx=10, pady=(0, 5), sticky="w")

        self.entry_mw = ctk.CTkEntry(tab, placeholder_text="e.g., 335.35")
        self.entry_mw.grid(row=4, column=0, padx=10, pady=(0, 10), sticky="ew")

        self.name_label = ctk.CTkLabel(tab, text="Full Name:")
        self.name_label.grid(row=5, column=0, padx=10, pady=(0, 5), sticky="w")

        self.entry_name = ctk.CTkEntry(
            tab, placeholder_text="e.g., Fmoc-Pra-OH; [0.40M]"
        )
        self.entry_name.grid(row=6, column=0, padx=10, pady=(0, 10), sticky="ew")

        self.submit_button_add_aa = ctk.CTkButton(
            tab, text="Add Amino Acid", command=self.add_amino_acid
        )
        self.submit_button_add_aa.grid(row=7, column=0, padx=10, pady=10)

        self.view_button_aa = ctk.CTkButton(
            tab,
            text="View Current Amino Acids",
            command=self.view_amino_acids,
        )
        self.view_button_aa.grid(row=8, column=0, padx=10, pady=10)

    def process_sequence(self) -> None:
        """Validate a peptide sequence and generate vial and synthesis plans."""
        try:
            sequence = self.entry.get().strip()
            calc = ValidatePeptide()
            tokens, original_tokens, invalid_amino_acids = calc.validate_user_sequence(sequence)

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

            validated_mass = calc.calculate_sequence_mass()
            synthesis_plan = BuildSynthesisPlan(tokens, original_tokens)
            df_vial_plan, vial_map = synthesis_plan.vial_rack_positions(tokens)
            df_synth_plan = synthesis_plan.build_synthesis_plan(vial_map)

            CTkMessagebox(title="Save File", message="Click OK to save vial plan", icon="info")
            initial_vial_path = LoadFile.resource_path("vial plan.csv")
            vial_plan_path = filedialog.asksaveasfilename(
                initialdir=os.path.dirname(initial_vial_path),
                initialfile=os.path.basename(initial_vial_path),
                title="Save Vial Plan CSV",
                defaultextension=".csv",
                filetypes=(("CSV files", "*.csv"), ("All files", "*.*")),
            )
            if not vial_plan_path:
                return

            CTkMessagebox(title="Save File", message="Click OK to save synthesis plan", icon="info")
            initial_synthesis_path = LoadFile.resource_path("synthesis plan.csv")
            synthesis_plan_path = filedialog.asksaveasfilename(
                initialdir=os.path.dirname(initial_synthesis_path),
                initialfile=os.path.basename(initial_synthesis_path),
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
                "Vial and synthesis plans saved successfully.\n"
                f"{vial_plan_path}\n{synthesis_plan_path}\n",
            )

        except Exception as e:
            CTkMessagebox(title="Error", message=f"An error occurred: {e}", icon="cancel")

    def process_compared_sequences(self) -> None:
        """Compare modified and previous peptide sequences, updating vial/synthesis plans."""
        try:
            new_sequence = self.entry_modify.get().strip()

            calc = ValidatePeptide()
            tokens, original_tokens, invalid_amino_acids = calc.validate_user_sequence(new_sequence)

            if not tokens or invalid_amino_acids:
                CTkMessagebox(title="Error", message="Invalid sequence.", icon="cancel")
                return

            validated_mass = calc.calculate_sequence_mass()
            builder_instance = BuildSynthesisPlan(tokens, original_tokens)

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
            initial_vial_path = LoadFile.resource_path("new vial plan.csv")
            vial_plan_path = filedialog.asksaveasfilename(
                initialdir=os.path.dirname(initial_vial_path),
                initialfile=os.path.basename(initial_vial_path),
                title="Save New Vial Plan CSV",
                defaultextension=".csv",
                filetypes=(("CSV files", "*.csv"), ("All files", "*.*")),
            )
            if not vial_plan_path:
                return

            CTkMessagebox(title="Save", message="Save updated synthesis plan", icon="info")
            initial_synthesis_path = LoadFile.resource_path("new synthesis plan.csv")
            synthesis_plan_path = filedialog.asksaveasfilename(
                initialdir=os.path.dirname(initial_synthesis_path),
                initialfile=os.path.basename(initial_synthesis_path),
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
                CTkMessagebox(
                    title="Error",
                    message="Molecular weight must be numeric.",
                    icon="cancel",
                )
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
                    f"Code: {row['AA']}\nMW: {row['MW']}\nName: {row['Name']}\n{'-' * 30}\n",
                )

            self.output_text.insert("end", f"\nTotal amino acids: {len(df)}")

        except Exception as e:
            CTkMessagebox(title="Error", message=f"Error loading amino acids: {e}", icon="cancel")