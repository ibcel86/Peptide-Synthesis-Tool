# Peptide Synthesis Tool

A Python application for generating automated peptide synthesis plans and vial maps for flow-based peptide synthesizers.

This tool validates peptide sequences, calculates molecular mass, assigns reagent vials based on usage constraints, and generates synthesis plans compatible with automated synthesis systems.

The project demonstrates modular software architecture, domain-specific scientific computation, and automated testing.

---

## Features

- Validate peptide sequences against a defined amino acid database
- Tokenize peptide sequences including multi-character amino acid codes
- Calculate total peptide molecular mass
- Automatically allocate reagent vials based on occurrence limits
- Split amino acids across multiple vials when capacity is exceeded
- Generate full synthesis plans for automated peptide synthesizers
- Compare modified peptide sequences against previous synthesis plans
- Update vial maps when sequences change
- Load and validate amino acid data from CSV

---

## Example Workflow

1. User provides peptide sequence

CPKTTAC


2. Sequence is validated against the amino acid database

3. Molecular mass is calculated

4. A vial map is generated

Example:

| Amino Acid | Rack | Position | Occurrences |
|------------|------|----------|-------------|
| C | 1 | 1 | 2 |
| P | 1 | 2 | 1 |
| K | 1 | 3 | 1 |

5. A synthesis plan is produced for automated execution

---

## Project Structure
```text
app/
├── core/
│ ├── sequence_processor.py
│ ├── sequence_comparator.py
│ └── synthesis_builder.py
│
├── io/
│ └── csv_loader.py
│
├── models/
│ └── amino_acids.py
│
├── ui/
│ └── gui components
│
tests/
└── unit tests for core logic
```

## Core Modules

| Module | Purpose |
|------|------|
| `sequence_processor` | peptide validation and molecular mass calculation |
| `synthesis_builder` | vial allocation and synthesis plan generation |
| `sequence_comparator` | detect sequence differences and update vial maps |
| `csv_loader` | load amino acid data and enforce CSV schema |

---

## Technologies

- Python
- Pandas
- Pytest
- Modular project architecture
- CSV-based data modelling

## Testing

The project includes automated unit tests covering:

- peptide validation
- molecular mass calculation
- vial allocation logic
- synthesis plan generation
- sequence comparison
- CSV data loading

Tests are written using **pytest**.

Run all tests:

```bash
python -m pytest -v
```
Example output:

============================= test session starts =============================

tests/test_compare_sequences.py .....
tests/test_data_loader.py ........
tests/test_synthesis_builder.py .....
tests/test_validate_peptide.py ....

============================== 21 passed =====================================

## All core scientific logic is validated by automated tests.

Why this project exists

This tool was developed to support peptide synthesis workflows by automating repetitive planning tasks such as vial allocation and synthesis plan generation.

The project demonstrates how domain expertise in chemistry can be combined with software engineering practices to build practical scientific tools.