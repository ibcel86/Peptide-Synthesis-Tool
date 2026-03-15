import pytest
from app.core.synthesis_builder import BuildSynthesisPlan

def test_vial_rack_positions_counts_occurrences_correctly():
    builder = BuildSynthesisPlan(tokens=["C", "C", "P"])
    df, vial_map = builder.vial_rack_positions(builder.tokens)

    c_row = df[df["Amino Acid"] == "C"].iloc[0]
    p_row = df[df["Amino Acid"] == "P"].iloc[0]

    assert c_row["Occurrences"] == 2
    assert p_row["Occurrences"] == 1

def test_vial_rack_positions_splits_when_occurrences_exceed_vial_capacity():
    builder = BuildSynthesisPlan(tokens=["C"] * 8)
    df, vial_map = builder.vial_rack_positions(builder.tokens)

    c_rows = df[df["Amino Acid"].str.startswith("C")]

    assert len(c_rows) == 2
    assert list(c_rows["Occurrences"]) == [6, 2]

def test_calculate_deprotection_vials_needed_returns_expected_value():
    builder = BuildSynthesisPlan(tokens=["C"] * 12)
    result = builder.calculate_deprotection_vials_needed(max_volume=16, inject_vol=1.5)

    assert result == 2

def test_build_synthesis_plan_creates_two_rows_per_token_when_mapping_exists():
    builder = BuildSynthesisPlan(tokens=["C", "P"])
    vial_map = {
        "C": (1, 1, 1),
        "P": (1, 2, 1),
    }

    df = builder.build_synthesis_plan(vial_map)

    assert len(df) == 4
    assert list(df["NAME"]) == [
        "C1",
        "deprotection 1",
        "P2",
        "deprotection 2",
    ]


def test_build_synthesis_plan_creates_error_row_if_no_vial_available():
    builder = BuildSynthesisPlan(tokens=["C"])
    vial_map = {}

    df = builder.build_synthesis_plan(vial_map)

    assert len(df) == 1
    assert df.iloc[0]["NAME"] == "ERROR_C"