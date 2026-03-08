from __future__ import annotations
import os
import sys
from typing import Dict
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