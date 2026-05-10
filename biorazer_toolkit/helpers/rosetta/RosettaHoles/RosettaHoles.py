"""

At step 4: Packing filters - Low threshold, designs with Rosetta holes score > 1.75 were filtered out. This is a very lenient threshold, and it is recommended to use a more stringent threshold for future designs.
At step 7: Packing filters - High threshold, designs with Rosetta holes score > 0.5 were filtered out. This is a more stringent threshold, and it is recommended to use an even more stringent threshold for future designs.

How to interpret the holes score output: 
A value of 0.0 means on par with native structures observed in the PDB; 
positive is worse (more voids), negative is better (less voids); 
this is for the default HolesParams used here, which is dec15; 
it is not recommended to change this unless you know the inner workings of what this code is doing.

References:
- /home/wangfanlin1/home/Repositories/PDE002/pipelines/pipeline-3.1.0/2_StructurePrediction/scripts/RosettaHoles
- Brunette, et al., Exploring the repeat protein universe through computational protein design, Nature, 2015
"""


from pathlib import Path
import pandas as pd
from biorazer_toolkit.apps.rosetta import RosettaApp
from biorazer_toolkit.apps.rosetta.scorefile import Scorefile

ROSETTA_SCRIPT_BIN = "/opt/rosetta.source.release-408/main/source/bin/rosetta_scripts.default.linuxgccrelease"
DALPHABALL_BIN = "/opt/rosetta.source.release-408/main/source/external/DAlpahBall/DAlphaBall.gcc"
script_dir = Path(__file__).parent
rosetta_script_xml = script_dir / "RosettaHoles.xml"
input_output_map = {
    "data/fastdesign_runs_1-comparison_107_83_0049_sample-0_model.pdb": "test/merged.sc",
}
log_path = "test/rosetta_holes.log"

app = RosettaApp(
    app_bin=ROSETTA_SCRIPT_BIN
)
app.set_default_handler(
    handler_types=["FileHandler"],
    file_path=log_path,
    file_mode="w"
)

for k, v in input_output_map.items():
    app.run(
        "-holes:dalphaball", f"{DALPHABALL_BIN}",
        "-s", f"{k}",
        "-parser:protocol", f"{rosetta_script_xml}",
        "-nstruct", "1",
        "-out:file:scorefile", f"{v}",
        "-out:file:pdb", f"{str(Path(v).parent)}",
        mode="subprocess.run"
    )

scorefiles = [Scorefile.from_sc(v) for v in input_output_map.values()]
result_df = pd.concat([scorefile.data for scorefile in scorefiles], ignore_index=True)
