import pandas as pd

from app.io.csv_loader import LoadFile, DataLoader, load_amino_acids


def test_ensure_csv_schema_creates_file_if_missing(tmp_path, monkeypatch):
    csv_path = tmp_path / "amino_acids.csv"

    monkeypatch.setattr(LoadFile, "get_csv_path", classmethod(lambda cls: str(csv_path)))

    result = LoadFile.ensure_csv_schema()

    assert result == str(csv_path)
    assert csv_path.exists()

    df = pd.read_csv(csv_path)
    assert list(df.columns) == ["AA", "MW", "Name"]
    assert df.empty


def test_ensure_csv_schema_adds_missing_columns_and_orders_them(tmp_path, monkeypatch):
    csv_path = tmp_path / "amino_acids.csv"

    pd.DataFrame({
        "AA": ["C", "P"],
        "Name": ["Cys", "Pro"],
    }).to_csv(csv_path, index=False)

    monkeypatch.setattr(LoadFile, "get_csv_path", classmethod(lambda cls: str(csv_path)))

    LoadFile.ensure_csv_schema()

    df = pd.read_csv(csv_path)
    assert list(df.columns) == ["AA", "MW", "Name"]
    assert df.loc[0, "AA"] == "C"
    assert df.loc[1, "AA"] == "P"


def test_load_amino_acids_returns_dictionary_keyed_by_code(tmp_path):
    csv_path = tmp_path / "amino_acids.csv"

    pd.DataFrame({
        "AA": ["C", "P"],
        "MW": [121.16, 115.13],
        "Name": ["Cys", "Pro"],
    }).to_csv(csv_path, index=False)

    amino_acids = load_amino_acids(str(csv_path))

    assert set(amino_acids.keys()) == {"C", "P"}
    assert amino_acids["C"].code == "C"
    assert amino_acids["C"].molecular_weight == 121.16
    assert amino_acids["C"].name == "Cys"


def test_load_amino_acids_skips_blank_codes(tmp_path):
    csv_path = tmp_path / "amino_acids.csv"

    pd.DataFrame({
        "AA": ["C", "", "P"],
        "MW": [121.16, 999.99, 115.13],
        "Name": ["Cys", "Ignore", "Pro"],
    }).to_csv(csv_path, index=False)

    amino_acids = load_amino_acids(str(csv_path))

    assert set(amino_acids.keys()) == {"C", "P"}


def test_load_amino_acids_strips_whitespace_from_fields(tmp_path):
    csv_path = tmp_path / "amino_acids.csv"

    pd.DataFrame({
        "AA": [" C ", " P "],
        "MW": [121.16, 115.13],
        "Name": [" Cys ", " Pro "],
    }).to_csv(csv_path, index=False)

    amino_acids = load_amino_acids(str(csv_path))

    assert set(amino_acids.keys()) == {"C", "P"}
    assert amino_acids["C"].name == "Cys"
    assert amino_acids["P"].name == "Pro"


def test_dataloader_loads_dataframe_and_lookup_structures(tmp_path, monkeypatch):
    csv_path = tmp_path / "amino_acids.csv"

    pd.DataFrame({
        "AA": ["C", "P", "TTAC"],
        "MW": [121.16, 115.13, 300.50],
        "Name": ["Cys", "Pro", "Special"],
    }).to_csv(csv_path, index=False)

    monkeypatch.setattr(LoadFile, "get_csv_path", classmethod(lambda cls: str(csv_path)))

    loader = DataLoader()

    assert not loader.df.empty
    assert set(loader.valid_amino_acids) == {"C", "P", "TTAC"}
    assert set(loader.amino_acids.keys()) == {"C", "P", "TTAC"}
    assert loader.mw_dict["C"] == 121.16
    assert loader.mw_dict["TTAC"] == 300.50


def test_dataloader_strips_whitespace_from_aa_column(tmp_path, monkeypatch):
    csv_path = tmp_path / "amino_acids.csv"

    pd.DataFrame({
        "AA": [" C ", " P "],
        "MW": [121.16, 115.13],
        "Name": ["Cys", "Pro"],
    }).to_csv(csv_path, index=False)

    monkeypatch.setattr(LoadFile, "get_csv_path", classmethod(lambda cls: str(csv_path)))

    loader = DataLoader()

    assert list(loader.df["AA"]) == ["C", "P"]
    assert "C" in loader.valid_amino_acids
    assert "P" in loader.valid_amino_acids


def test_ensure_csv_schema_preserves_existing_rows(tmp_path, monkeypatch):
    csv_path = tmp_path / "amino_acids.csv"

    pd.DataFrame({
        "AA": ["C"],
        "MW": [121.16],
        "Name": ["Cys"],
    }).to_csv(csv_path, index=False)

    monkeypatch.setattr(LoadFile, "get_csv_path", classmethod(lambda cls: str(csv_path)))

    LoadFile.ensure_csv_schema()

    df = pd.read_csv(csv_path)
    assert len(df) == 1
    assert df.loc[0, "AA"] == "C"
    assert df.loc[0, "MW"] == 121.16
    assert df.loc[0, "Name"] == "Cys"