from __future__ import annotations
from typing import List, Tuple
from app.io.csv_loader import DataLoader

class ValidatePeptide:
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