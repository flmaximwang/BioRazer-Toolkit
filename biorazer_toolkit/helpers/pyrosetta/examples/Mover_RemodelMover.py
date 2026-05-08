from pathlib import Path
import pyrosetta
from pyrosetta.io import pose_from_pdb
from pyrosetta.rosetta.basic import options
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.protocols.moves import PyMOLMover
from pyrosetta.rosetta.protocols.forge.remodel import RemodelMover
from biorazer_toolkit.apps.rosetta.blueprint import rosetta_config, Blueprint

pyrosetta.init()

pdb_init = "data/EXP/6C51_E4.pdb"
pose_init = pose_from_pdb(pdb_init)
pymol = PyMOLMover()
pymol.apply(pose_init)


rosetta_config.set_rosetta_dir(
    "/Users/maxim/Applications/Rosetta/rosetta.source.release-362"
)
bp = Blueprint.from_pdb(pdb_init, chain="E")
bp_data = bp.get_data()
loop_seqs = ["SSGG", "SSGSGSGG", "SGSG"]
loop_extensions = [(1, 1), (2, 1), (1, 1)]
break_points = [25, 50, 75]
for loop_seq, break_point, extend in zip(loop_seqs, break_points, loop_extensions):
    bp.insert_seq(loop_seq, break_point, ss="L", extend=extend)
bp.set_ss(49, "H")
bp_file = "1.3_3loops.bp"
bp.to_file(bp_file)

#
rm_mover = RemodelMover()
# options.set_real_option('run:show_simulation_in_pymol', 1.0)
options.set_file_option("remodel:blueprint", bp_file)
options.set_integer_option("remodel:dr_cycles", 3)  # Default
options.set_real_option(
    "remodel:RemodelLoopMover:max_linear_chainbreak", 0.07
)  # Default
options.set_boolean_option("remodel:quick_and_dirty", False)  # Default
options.set_boolean_option("remodel:design:find_neighbors", False)  # Default
options.set_boolean_option("remodel:design:design_neighbors", False)  # Default
rm_mover.register_options()
# rm_mover.max_linear_chainbreak(0.07) # Default
rm_mover.apply(pose_init)

pose_init.dump_pdb(str(Path(pdb_init).with_stem(Path(pdb_init).stem + "_with_loops_2")))
pymol.apply(pose_init)
