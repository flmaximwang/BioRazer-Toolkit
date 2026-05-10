from biorazer.structure.io.protein import PDB2SEQ
from pathlib import Path
import shutil, re


def reduce_structures(structure_dir: str | Path):
    """
    Remove redundant structures from the given directory.
    Redundant structures are those that have identical sequences.
    The first occurrence of each unique sequence is kept.
    """

    structure_dir = Path(structure_dir)
    output_dir = structure_dir.with_name(structure_dir.name + "_reduced")
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
    seq_dict = {}
    pdb_files = sorted(
        list(structure_dir.glob("*.pdb")),
        key=lambda x: int(re.match(r"match(\d+)", x.stem).group(1)),
    )
    for pdb_file in pdb_files:
        seq = list(PDB2SEQ(input_file=pdb_file).read().values())[0]
        if seq not in seq_dict.values():
            seq_dict[pdb_file.stem] = seq
            output_path = output_dir / pdb_file.name
            shutil.copyfile(pdb_file, output_path)
