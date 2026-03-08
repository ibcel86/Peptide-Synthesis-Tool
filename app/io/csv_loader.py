from __future__ import annotations
import csv
import os
import sys
from typing import Dict
import pandas as pd
from app.models.amino_acids import AminoAcid


class LoadFile:
    """Utility class for handling resource file paths and ensuring CSV schema."""

    @classmethod
    def resource_path(cls, relative_path: str) -> str:
        """Return absolute path to a resource, compatible with PyInstaller builds."""
        if getattr(sys, "frozen", False):
            return os.path.join(os.path.dirname(sys.executable), relative_path)
        return os.path.join(os.path.dirname(__file__), relative_path)

    @classmethod
    def get_csv_path(cls) -> str:
        """Return absolute path to amino_acids.csv."""
        return cls.resource_path("amino_acids.csv")

    @classmethod
    def ensure_csv_schema(cls) -> str:
        """Ensure amino_acids.csv exists and has the expected columns."""
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


def load_amino_acids(filepath: str) -> dict[str, AminoAcid]:
    """Load amino acids from CSV into a dictionary keyed by amino acid code."""
    amino_acids: dict[str, AminoAcid] = {}

    with open(filepath, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)

        for row in reader:
            code = row["AA"].strip()
            if not code:
                continue

            aa = AminoAcid(
                code=code,
                molecular_weight=float(row["MW"]),
                name=row["Name"].strip(),
            )
            amino_acids[aa.code] = aa

    return amino_acids


class DataLoader:
    """Load amino acid data from CSV into memory."""

    def __init__(self) -> None:
        path = LoadFile.ensure_csv_schema()

        self.df: pd.DataFrame = pd.read_csv(path)
        self.df["AA"] = self.df["AA"].astype(str).str.strip()

        self.amino_acids: dict[str, AminoAcid] = load_amino_acids(path)
        self.valid_amino_acids: set[str] = set(self.amino_acids.keys())
        self.mw_dict: Dict[str, float] = {
            code: aa.molecular_weight for code, aa in self.amino_acids.items()
        }