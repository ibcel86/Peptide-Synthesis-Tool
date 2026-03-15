import pytest
from app.core.sequence_processor import ValidatePeptide

def test_validate_user_sequence_returns_forward_and_reverse_tokens():
    validator = ValidatePeptide()
    reversed_tokens, original_tokens, invalid = validator.validate_user_sequence("CPK")

    assert reversed_tokens == ["K", "P", "C"]
    assert original_tokens == ["C", "P", "K"]
    assert invalid == []

def test_validate_user_sequence_raises_for_invalid_amino_acids():
    validator = ValidatePeptide()

    with pytest.raises(ValueError, match="Invalid amino acid"):
        validator.validate_user_sequence("CPZ")

def test_calculate_sequence_mass_returns_sum_of_token_masses():
    validator = ValidatePeptide()
    validator.validate_user_sequence("CP")

    expected = (
        validator.data.amino_acids["P"].molecular_weight
        + validator.data.amino_acids["C"].molecular_weight
    )

    assert validator.calculate_sequence_mass() == expected

def test_calculate_sequence_mass_raises_if_sequence_not_validated():
    validator = ValidatePeptide()

    with pytest.raises(ValueError, match="No sequence loaded"):
        validator.calculate_sequence_mass()