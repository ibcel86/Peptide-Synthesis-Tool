# Peptide Synthesis Tool

## Automated Workflow Engine for Laboratory Synthesis Planning

A Python desktop application that replaces a multi-day manual calculation 
process with a GUI-driven workflow, producing validated, machine-readable 
output files in seconds.

Built to solve a real production problem: scientists manually calculating 
amino acid quantities for peptide synthesis were spending days on error-prone 
spreadsheet work. A single miscalculation could invalidate weeks of 
downstream lab work. This tool eliminates that entirely.

---

## Tech Stack

- **Language**: Python
- **GUI**: Tkinter
- **Output**: CSV / Excel via automated file generation
- **Packaging**: PyInstaller (.exe deployment)
- **Data**: User-configurable amino acid database (CSV)

---

## What It Does

- Accepts a peptide sequence as input via GUI
- Calculates reagent quantities, vial counts, and rack layouts automatically
- Handles hardware constraints (vial volume limits, rack sizing, sampling order)
- Exports two machine-readable files: a vial plan and a synthesis plan
- v2: reverse-engineered legacy machine output format for compatibility with 
  older synthesis hardware
- v2: sequence comparison tool that diffs two synthesis runs and generates 
  an incremental plan — automating another previously manual workflow

---

## Impact

- Reduced synthesis planning time from **days to seconds**
- Eliminated manual calculation errors that previously caused weeks of lost work
- Deployed and used in active drug discovery research
