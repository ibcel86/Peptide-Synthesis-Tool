"""
peptidesequencetool.py

Core logic for the Peptide Sequence Tool.

This module provides classes for peptide sequence validation, molecular mass
calculation, vial mapping, and synthesis plan generation. It also includes
functionality for comparing old and modified sequences and updating vial plans
accordingly.
"""

from __future__ import annotations
import os
import sys
import re
from math import ceil, floor
from collections import Counter
from typing import Any, Dict, List, Tuple

import pandas as pd


class LoadFile:
    """Utility class for handling resource file paths and ensuring CSV schema."""

    @classmethod
    def resource_path(cls, relative_path: str) -> str:
        """Return the absolute path to a resource, compatible with PyInstaller builds.

        Args:
            relative_path (str): Path relative to the script or executable.

        Returns:
            str: Absolute path to the resource.
        """
        if getattr(sys, "frozen", False):
            return os.path.join(os.path.dirname(sys.executable), relative_path)
        return os.path.join(os.path.dirname(__file__), relative_path)

    @classmethod
    def get_csv_path(cls) -> str:
        """Return the absolute path to the amino acid CSV file.

        Returns:
            str: Full path to amino_acids.csv.
        """
        return cls.resource_path("amino_acids.csv")

    @classmethod
    def ensure_csv_schema(cls) -> str:
        """Ensure `amino_acids.csv` exists and follows the expected schema.

        If the file does not exist it is created with columns ['AA','MW','Name'].
        Missing columns are added and the order is normalized when the file exists.

        Returns:
            str: Path to the validated or newly created CSV file.
        """
        path = cls.get_csv_path()
        if not os.path.exists(path):
            pd.DataFrame(columns=["AA", "MW", "Name"]).to_csv(path, index=False)
            return path

        df = pd.read_csv(path)
        for col in ["AA", "MW", "Name"]:
            if col not in df.columns:
                df[col] = pd.Series(dtype="object")
        df = df[["AA", "MW", "Name"]]
        df.to_csv(path, index=False)
        return path


class DataLoader:
    """Load amino acid data from CSV into memory."""

    def __init__(self) -> None:
        path = LoadFile.ensure_csv_schema()
        self.df: pd.DataFrame = pd.read_csv(path)
        self.df["AA"] = self.df["AA"].astype(str).str.strip()
        self.valid_amino_acids: set[str] = set(self.df["AA"])
        self.mw_dict: Dict[str, float] = dict(zip(self.df["AA"], self.df["MW"]))


class CalculatePeptide:
    """Validate peptide sequences and calculate molecular mass."""

    def __init__(self) -> None:
        self.data = DataLoader()
        self.tokens: List[str] | None = None
        self.original_tokens: List[str] | None = None

    def _tokenize_sequence(self, sequence: str) -> List[str]:
        """Tokenize a sequence using known amino acid codes.

        Args:
            sequence (str): Raw peptide sequence (no spaces).

        Returns:
            List[str]: Tokenized list of amino acids.
        """
        valid_aas = sorted(self.data.valid_amino_acids, key=len, reverse=True)
        tokens: List[str] = []
        i = 0
        while i < len(sequence):
            match = None
            for aa in valid_aas:
                if sequence[i:].startswith(aa):
                    match = aa
                    tokens.append(match)
                    i += len(match)
                    break
            if not match:
                tokens.append(sequence[i])
                i += 1
        return tokens

    def validate_user_sequence(self, sequence: str) -> Tuple[List[str], List[str], List[str]]:
        """Validate and tokenize a user-provided sequence.

        Sequence should contain no spaces and use one-letter or multi-letter
        amino acid codes present in `amino_acids.csv`.

        Args:
            sequence (str): Peptide sequence without spaces.

        Returns:
            Tuple[List[str], List[str], List[str]]:
                - tokens: sequence reversed for synthesis order
                - original_tokens: original order tokens
                - invalid_amino_acids: list of invalid tokens (if any)

        Raises:
            ValueError: If an invalid amino acid is encountered.
        """
        sequence = sequence.strip()
        valid_aas = sorted(self.data.valid_amino_acids, key=len, reverse=True)
        tokens: List[str] = []
        i = 0

        while i < len(sequence):
            match = None
            for aa in valid_aas:
                if sequence.startswith(aa, i):
                    match = aa
                    tokens.append(match)
                    i += len(aa)
                    break
            if not match:
                raise ValueError(f"Invalid amino acid at position {i+1}: '{sequence[i:]}'")

        self.original_tokens = tokens
        self.tokens = tokens[::-1]
        invalid_amino_acids = [aa for aa in tokens if aa not in self.data.valid_amino_acids]
        return self.tokens, self.original_tokens, invalid_amino_acids

    def calculate_sequence_mass(self, sequence: str) -> float:
        """Calculate total molecular mass of a validated peptide sequence.

        Args:
            sequence (str): The validated peptide sequence (for clarity).

        Returns:
            float: Molecular mass in g/mol.

        Raises:
            ValueError: If no validated sequence has been loaded.
        """
        if not self.tokens:
            raise ValueError("No sequence loaded. Run validate_user_sequence() first.")
        return sum(self.data.mw_dict[aa] for aa in self.tokens)


class BuildSynthesisPlan:
    """Generate vial mappings and synthesis plans for automated peptide synthesis."""

    def __init__(self, tokens: List[str], original_tokens: List[str] | None = None) -> None:
        self.data = DataLoader()
        self.tokens = tokens
        self.original_tokens = original_tokens or tokens

    def vial_rack_positions(
        self,
        tokens: List[str],
        conc: float = 0.4,
        max_occurrence: int = 6,
        max_volume: int = 16,
        start_rack: int = 1,
        start_position: int = 1,
    ) -> Tuple[pd.DataFrame, Dict[str, Tuple[int, int, int]]]:
        """Generate vial map and rack positions.

        Args:
            tokens (List[str]): Amino acid tokens (synthesis order or otherwise).
            conc (float): Stock concentration (M).
            max_occurrence (int): Max occurrences per vial.
            max_volume (int): Max vial volume (mL).
            start_rack (int): Starting rack index.
            start_position (int): Starting position index.

        Returns:
            Tuple[pd.DataFrame, dict]: DataFrame of vial map and mapping dict {name: (rack, pos, occ)}.
        """
        amino_acid_occurrences = Counter(tokens)
        max_per_vial = floor(max_volume / 2.5)
        output: List[Dict[str, Any]] = []
        vial_map: Dict[str, Tuple[int, int, int]] = {}

        rack = start_rack
        position = start_position
        max_positions = 27

        for aa, count in amino_acid_occurrences.items():
            mw = self.data.mw_dict[aa]
            splits: List[int] = []
            while count > 0:
                chunk = min(count, max_per_vial)
                splits.append(chunk)
                count -= chunk

            for i, split_count in enumerate(splits):
                name = aa if i == 0 else f"{aa}{i+1}"
                mmol = split_count * ((max_volume * conc) / max_occurrence)
                volume = split_count * 2.5
                mass = mmol * mw / 1000

                output.append(
                    {
                        "Amino Acid": name,
                        "Rack": rack,
                        "Position": position,
                        "Occurrences": split_count,
                        "mmol": round(mmol, 2),
                        "Mass (g)": round(mass, 2),
                        "Volume (mL)": round(volume, 2),
                    }
                )
                vial_map[name] = (rack, position, split_count)
                position += 1
                if position > max_positions:
                    rack += 1
                    position = 1

        return pd.DataFrame(output), vial_map

    def calculate_deprotection_vials_needed(self, max_volume: int = 16, inject_vol: float = 1.5) -> int:
        """Calculate number of deprotection vials required.

        Args:
            max_volume (int): Maximum vial volume (mL).
            inject_vol (float): Injection volume (mL).

        Returns:
            int: Number of deprotection vials.
        """
        num_deprotection_steps = len(self.tokens)
        samples_per_vial = ceil(max_volume / inject_vol)
        return ceil(num_deprotection_steps / samples_per_vial)

    def build_synthesis_plan(self, vial_map: Dict[str, Tuple[int, int, int]], max_deprotection_volume: int = 16) -> pd.DataFrame:
        """Build a synthesis plan DataFrame based on vial mapping.

        Args:
            vial_map (dict): Mapping of amino acid name -> (rack, position, occurrences).
            max_deprotection_volume (int): Max vial volume (mL).

        Returns:
            pd.DataFrame: Synthesis plan suitable for export.

        Raises:
            ValueError: If insufficient rack space exists for deprotection vials.
        """
        num_deprotection_vials = self.calculate_deprotection_vials_needed(max_deprotection_volume)
        deprotection_start_pos = 28
        rack2_end_pos = 54

        last_position_needed = deprotection_start_pos + num_deprotection_vials - 1
        if last_position_needed > rack2_end_pos:
            available_positions = rack2_end_pos - deprotection_start_pos + 1
            raise ValueError(
                f"Not enough rack space for deprotection vials. Need {num_deprotection_vials}, available {available_positions}."
            )

        synthesis_rows: List[Dict[str, Any]] = []
        vial_usage_counter: Dict[str, int] = {}
        deprotection_usage_counter = 0
        deprotection_positions = [deprotection_start_pos + i for i in range(num_deprotection_vials)]
        uses_per_deprotection_vial = ceil(len(self.tokens) / num_deprotection_vials)

        for synthesis_position, aa in enumerate(self.tokens, 1):
            related_vials = [
                v for v in vial_map.keys() if v == aa or (v.startswith(aa) and v[len(aa):].isdigit())
            ]
            related_vials.sort(key=lambda x: (0 if x == aa else int(x[len(aa):])))

            deprotection_vial_index = min(
                deprotection_usage_counter // uses_per_deprotection_vial,
                len(deprotection_positions) - 1,
            )
            current_deprotection_pos = deprotection_positions[deprotection_vial_index]

            assigned = False
            for vial_name in related_vials:
                rack, pos, occ = vial_map[vial_name]
                used = vial_usage_counter.get(vial_name, 0)

                if used < occ:
                    vial_usage_counter[vial_name] = used + 1
                    synthesis_rows.append(
                        {
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
                            "MANUAL CLEAN (ml)": 4,
                        }
                    )
                    synthesis_rows.append(
                        {
                            "NAME": f"deprotection {synthesis_position}",
                            "FLOW RATE A (ml/min)": 0,
                            "FLOW RATE B (ml/min)": 0,
                            "FLOW RATE D (ml/min)": 0.8,
                            "RESIDENCE 2": True,
                            "AUTOSAMPLER SITE A": pos,
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
                            "MANUAL CLEAN (ml)": 4,
                        }
                    )
                    deprotection_usage_counter += 1
                    assigned = True
                    break

            if not assigned:
                synthesis_rows.append(
                    {
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
                        "MANUAL CLEAN (ml)": 0,
                    }
                )

        return pd.DataFrame(synthesis_rows)


class CompareSequences:
    """Compare and update vial maps / synthesis plans after sequence modifications."""

    def __init__(self, builder_instance: BuildSynthesisPlan, old_synthesis_path: str, old_vial_path: str) -> None:
        self.builder = builder_instance
        self.old_synthesis_path = old_synthesis_path
        self.old_vial_path = old_vial_path
        self.tokens: List[str] | None = None
        self.original_tokens: List[str] | None = None
        self.data = DataLoader()

    def extract_old_sequence_from_csv(self) -> List[str]:
        """Extract old peptide sequence tokens from an existing synthesis plan CSV.

        Returns:
            List[str]: Tokens from the old sequence in forward order.

        Raises:
            FileNotFoundError: If the synthesis plan file cannot be found.
        """
        if not os.path.exists(self.old_synthesis_path):
            raise FileNotFoundError("Synthesis plan not found, please ensure the file is accessible.")

        df = pd.read_csv(self.old_synthesis_path)
        df.columns = df.columns.str.strip()
        aa_rows = df[~df["NAME"].str.contains("deprotection", case=False, na=False)]
        cleaned_tokens = [re.sub(r"\d+$", "", name.strip()) for name in aa_rows["NAME"]]
        self.original_tokens = cleaned_tokens[::-1]
        return cleaned_tokens

    def compare_sequences(self, cleaned_tokens: List[str], new_aa: List[str]) -> List[str]:
        """Return amino acids present in the new sequence that differ from the old.

        Args:
            cleaned_tokens (List[str]): Tokens from the old sequence (forward order).
            new_aa (List[str]): Tokens from the new sequence (forward order).

        Returns:
            List[str]: Amino acids that need to be added to the vial map.
        """
        differences: List[str] = [new for old, new in zip(cleaned_tokens, new_aa) if old != new]
        if len(new_aa) > len(cleaned_tokens):
            differences.extend(new_aa[len(cleaned_tokens):])
        return differences

    def build_new_vial_map(self, new_aa: List[str]) -> pd.DataFrame:
        """Build an updated vial map by appending new amino acids to the existing vial map CSV.

        Args:
            new_aa (List[str]): New amino acids that should be added.

        Returns:
            pd.DataFrame: Combined DataFrame of the old and new vial mappings.

        Raises:
            FileNotFoundError: If the old vial map CSV cannot be found.
        """
        if not os.path.exists(self.old_vial_path):
            raise FileNotFoundError("Vial map not found. Please ensure the file is accessible.")

        df_old = pd.read_csv(self.old_vial_path)
        df_old.columns = df_old.columns.str.strip()

        last_row = df_old.loc[df_old["Rack"].idxmax()]
        last_rack = int(last_row["Rack"])
        last_position = int(df_old[df_old["Rack"] == last_rack]["Position"].max())

        max_positions = 27
        max_per_vial = 6
        if last_position >= max_positions:
            start_rack = last_rack + 1
            start_position = 1
        else:
            start_rack = last_rack
            start_position = last_position + 1

        pattern = re.compile(r"^([A-Za-z]+)(\d+)?$")
        aa_max_index: Dict[str, int] = {}
        for name in df_old["Amino Acid"]:
            match = pattern.match(str(name))
            if match:
                base = match.group(1)
                idx = int(match.group(2)) if match.group(2) else 1
                aa_max_index[base] = max(aa_max_index.get(base, 0), idx)

        cleaned_new_aa = [aa.replace("*", "") for aa in new_aa]
        new_occurrences = Counter(cleaned_new_aa)

        output: List[Dict[str, Any]] = []
        rack = start_rack
        position = start_position

        for aa in cleaned_new_aa:
            if new_occurrences[aa] == 0:
                continue

            total_count = new_occurrences[aa]
            new_occurrences[aa] = 0
            splits: List[int] = []
            while total_count > 0:
                chunk = min(total_count, max_per_vial)
                splits.append(chunk)
                total_count -= chunk

            start_index = aa_max_index.get(aa, 0)

            for i, split_count in enumerate(splits):
                suffix = "" if start_index == 0 and i == 0 else str(start_index + i + 1)
                name = f"{aa}{suffix}"
                mmol = split_count * (16 * 0.4) / 6
                mass = mmol * self.data.mw_dict.get(aa, 0) / 1000
                volume = split_count * 2.5

                output.append(
                    {
                        "Amino Acid": name,
                        "Rack": rack,
                        "Position": position,
                        "Occurrences": split_count,
                        "mmol": round(mmol, 2),
                        "Mass (g)": round(mass, 2),
                        "Volume (mL)": round(volume, 2),
                    }
                )

                position += 1
                if position > max_positions:
                    rack += 1
                    position = 1

        df_new = pd.DataFrame(output)
        df_combined = pd.concat([df_old, df_new], ignore_index=True)
        return df_combined

    def build_new_synthesis_plan(self, df_combined: pd.DataFrame) -> pd.DataFrame:
        """Build a new synthesis plan DataFrame using the updated combined vial map.

        Args:
            df_combined (pd.DataFrame): Combined vial map (old + appended new).

        Returns:
            pd.DataFrame: Updated synthesis plan.
        """
        vial_map: Dict[str, Tuple[int, int, int]] = {
            row["Amino Acid"]: (int(row["Rack"]), int(row["Position"]), int(row["Occurrences"]))
            for _, row in df_combined.iterrows()
        }
        builder = BuildSynthesisPlan(self.tokens or [], self.original_tokens or [])
        return builder.build_synthesis_plan(vial_map)
