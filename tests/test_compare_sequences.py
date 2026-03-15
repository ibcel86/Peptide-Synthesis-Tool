import pytest
from app.core.sequence_comparator import CompareSequences
from app.core.synthesis_builder import BuildSynthesisPlan

def test_compare_sequences_returns_changed_amino_acids():
    builder = BuildSynthesisPlan(tokens=["C"])
    comparer = CompareSequences(builder, "old_plan.csv", "old_vial.csv")

    result = comparer.compare_sequences(
        cleaned_tokens=["C", "P", "K"],
        new_aa=["C", "V", "K"],
    )

    assert result == ["V"]

def test_compare_sequences_appends_extra_amino_acids_when_new_sequence_is_longer():
    builder = BuildSynthesisPlan(tokens=["C"])
    comparer = CompareSequences(builder, "old_plan.csv", "old_vial.csv")

    result = comparer.compare_sequences(
        cleaned_tokens=["C", "P"],
        new_aa=["C", "P", "K", "V"],
    )

    assert result == ["K", "V"]

def test_extract_old_sequence_from_csv_raises_if_file_missing():
    builder = BuildSynthesisPlan(tokens=["C"])
    comparer = CompareSequences(builder, "does_not_exist.csv", "old_vial.csv")

    with pytest.raises(FileNotFoundError, match="Synthesis plan not found"):
        comparer.extract_old_sequence_from_csv()


def test_build_new_vial_map_raises_if_vial_map_missing():
    builder = BuildSynthesisPlan(tokens=["C"])
    comparer = CompareSequences(builder, "old_plan.csv", "missing_vial.csv")

    with pytest.raises(FileNotFoundError, match="Vial map not found"):
        comparer.build_new_vial_map(["C"])

