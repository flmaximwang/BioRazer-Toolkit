"""
对于一个模型进行基本的评估
- 是否有空间冲突, 冲突出现的位置
- 整体的能量水平, 能量较高的局部位置
"""

from pathlib import Path
from biorazer.structure.io.protein import PDB2STRUCT
from biorazer.structure.analysis.static.report import report_intra_steric_clashes
from biorazer.structure.analysis.static.check import is_steric_clashed
from biorazer_ex.apps.rosetta.execution import RosettaApp


structure = PDB2STRUCT(input_file="fixbb-2025.10.27/pdbs/1.pdb").read()
print(f"Is the structure steric clashed? {is_steric_clashed(structure)}")

# Report steric clashes in the structure, selecting all atoms, and outputting in PyMOL format.
report_intra_steric_clashes(structure, selection="all", fmt="pymol")

# Report overall energy of the structure, and output in a simple text format.
total_energy_app = RosettaApp(app_dir="/Users/maxim/Applications/Rosetta")
total_energy_app.use_app(["default", "score_jd2"])
total_energy_app.run_with_structure(
    structure, "-score:weights", "ref2015", "-out:file:scorefile", "total_energy.sc"
)

# Report per-residue energy breakdown, using the default energy function, and output in a silent file format.
energy_breakdown_app = RosettaApp(app_dir="/Users/maxim/Applications/Rosetta")
energy_breakdown_app.use_app(["default", "residue_energy_breakdown"])
energy_breakdown_app.run_with_structure(
    structure, "-score:weights", "ref2015", "-out:file:silent", "energy_breakdown.sc"
)
