from __future__ import annotations
import os
import re
from collections import Counter
from typing import Any, Dict, List, Tuple
import pandas as pd
from app.io.csv_loader import DataLoader
from app.core.synthesis_builder import BuildSynthesisPlan
class CompareSequences:
    """Compare and update vial maps / synthesis plans after sequence modifications."""

    def __init__(
        self,
        builder_instance: BuildSynthesisPlan,
        old_synthesis_path: str,
        old_vial_path: str,
    ) -> None:
        self.builder = builder_instance
        self.old_synthesis_path = old_synthesis_path
        self.old_vial_path = old_vial_path
        self.tokens: List[str] | None = None
        self.original_tokens: List[str] | None = None
        self.data = DataLoader()

    def extract_old_sequence_from_csv(self) -> List[str]:
        """Extract old peptide sequence tokens from an existing synthesis plan CSV."""
        if not os.path.exists(self.old_synthesis_path):
            raise FileNotFoundError(
                "Synthesis plan not found, please ensure the file is accessible."
            )

        df = pd.read_csv(self.old_synthesis_path)
        df.columns = df.columns.str.strip()
        aa_rows = df[~df["NAME"].str.contains("deprotection", case=False, na=False)]
        cleaned_tokens = [re.sub(r"\d+$", "", name.strip()) for name in aa_rows["NAME"]]
        self.original_tokens = cleaned_tokens[::-1]
        return cleaned_tokens

    def compare_sequences(self, cleaned_tokens: List[str], new_aa: List[str]) -> List[str]:
        """Return amino acids present in the new sequence that differ from the old."""
        differences: List[str] = [new for old, new in zip(cleaned_tokens, new_aa) if old != new]
        if len(new_aa) > len(cleaned_tokens):
            differences.extend(new_aa[len(cleaned_tokens):])
        return differences

    def build_new_vial_map(self, new_aa: List[str]) -> pd.DataFrame:
        """Build an updated vial map by appending new amino acids to the existing vial map CSV."""
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
                mass = mmol * self.data.amino_acids[aa].molecular_weight / 1000
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
        """Build a new synthesis plan DataFrame using the updated combined vial map."""
        vial_map: Dict[str, Tuple[int, int, int]] = {
            row["Amino Acid"]: (int(row["Rack"]), int(row["Position"]), int(row["Occurrences"]))
            for _, row in df_combined.iterrows()
        }
        builder = BuildSynthesisPlan(self.tokens or [], self.original_tokens or [])
        return builder.build_synthesis_plan(vial_map)