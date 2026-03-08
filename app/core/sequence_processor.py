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
        """Tokenize a sequence using known amino acid codes."""
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
                tokens.append(sequence[i])
                i += 1

        return tokens

    def validate_user_sequence(self, sequence: str) -> Tuple[List[str], List[str], List[str]]:
        """Validate and tokenize a user-provided sequence."""
        sequence = sequence.strip()
        tokens = self._tokenize_sequence(sequence)

        invalid_amino_acids = [aa for aa in tokens if aa not in self.data.valid_amino_acids]
        if invalid_amino_acids:
            raise ValueError(
                f"Invalid amino acid(s) found: {', '.join(invalid_amino_acids)}"
            )

        self.original_tokens = tokens
        self.tokens = tokens[::-1]
        return self.tokens, self.original_tokens, invalid_amino_acids

    def calculate_sequence_mass(self) -> float:
        """Calculate total molecular mass of the validated peptide sequence."""
        if not self.tokens:
            raise ValueError("No sequence loaded. Run validate_user_sequence() first.")

        return sum(
            self.data.amino_acids[aa].molecular_weight
            for aa in self.tokens
        )