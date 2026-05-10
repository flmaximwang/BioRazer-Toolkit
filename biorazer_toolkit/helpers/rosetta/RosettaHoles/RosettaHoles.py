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
        mode="subprocess.run"
    )

scorefiles = [Scorefile.from_sc(v) for v in input_output_map.values()]
result_df = pd.concat([scorefile.data for scorefile in scorefiles], ignore_index=True)
