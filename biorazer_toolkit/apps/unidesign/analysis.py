from pathlib import Path
import pandas as pd


def parse_summary(summary_file: str | Path) -> pd.DataFrame:
    """
    Parse the UniDesign_bestseqs.txt summary file into a pandas DataFrame.
    """
    col_names = [
        "sequence",
        "beststruct_pdb",
        "seq_recovery",
        "total_energy",
        "evolutionary_energy",
        "physical_energy_without_binding",
        "binding_energy",
        "unsatisfied_catalytic_constraints",
    ]
    all_data = []
    with open(summary_file, "r") as file:
        for line in file:
            if line.startswith("#"):
                continue
            cols = line.strip().split()
            row_data = dict()
            for i, col_name in enumerate(col_names):
                row_data[col_name] = cols[i]
            all_data.append(row_data)
    df = pd.DataFrame(all_data)
    df["unsatisfied_catalytic_constraints"] = df[
        "unsatisfied_catalytic_constraints"
    ].astype(int)
    return df


def fix_unidesign_pdb(pdb_file: str | Path):
    """
    PDBs from UniDesign are not in the correct format for biotite.
    This function replaces ATOM lines[60: 66] with '  0.00'
    """
    pdb_file = Path(pdb_file)
    with open(pdb_file, "r") as file:
        lines = file.readlines()

    for i, line in enumerate(lines):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            lines[i] = line[:60] + "  0.00" + line[66:]

    fixed_pdb = pdb_file.with_stem(pdb_file.stem + "_fixed")
    with open(fixed_pdb, "w") as file:
        file.writelines(lines)

    return fixed_pdb
