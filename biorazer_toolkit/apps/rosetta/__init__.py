from pathlib import Path


class RosettaConfig:
    def __init__(self):
        self.ROSETTA_DIR = ""
        self.ROSETTA_BIN = ""
        self.ROSETTA_TOOLS = ""

    def set_rosetta_dir(self, path):
        self.ROSETTA_DIR = str(Path(path))
        self.ROSETTA_BIN = f"{self.ROSETTA_DIR}/main/source/bin"
        self.ROSETTA_TOOLS = f"{self.ROSETTA_DIR}/main/tools"
        self.check()

    def check(self):
        try:
            list(Path(self.ROSETTA_BIN).glob("rosetta_scripts*"))[0]
            # print(f"Rosetta directory {self.ROSETTA_DIR} is set correctly.")
        except IndexError:
            raise FileNotFoundError(
                f"Rosetta is not installed at {self.ROSETTA_BIN}.\n"
                "Please set the Rosetta directory first.\n"
                "You can use <module>.rosetta_config.set_rosetta_dir(path).\n"
                "A correct directory looks like /Users/maxim/Applications/Rosetta/rosetta.source.release-362"
            )


config = RosettaConfig()
rosetta_config = config

from .scorefile import Scorefile
from .execution import RosettaApp
