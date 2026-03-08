from __future__ import annotations
from math import ceil, floor
from collections import Counter
from typing import Any, Dict, List, Tuple
import pandas as pd
from app.io.csv_loader import DataLoader

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