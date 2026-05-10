from pathlib import Path
from ..basic import App


class UniDesignConfig(App):

    def run_with_structure(
        self,
        atom_array,
        *args,
        input_file_format="pdb",
        **kwargs,
    ):
        raise NotImplementedError(
            "UniDesignConfig is a configuration holder and does not execute structure-based runs."
        )

    def set(self, app_dir):
        self.DIR = app_dir
        self.BIN = f"{self.DIR}/UniDesign"
        self.LIBRARY = f"{self.DIR}/library"
        self.WREAD = f"{self.DIR}/wread"

    def check(self):
        assert Path(self.BIN).exists(), f"UniDesign binary not found at {self.BIN}"
        assert Path(
            self.LIBRARY
        ).exists(), f"UniDesign library not found at {self.LIBRARY}"
        assert Path(self.WREAD).exists(), f"UniDesign wread not found at {self.WREAD}"

    def get(self):
        """

        Return
        ------
        dict
            Dictionary with the config parameters.
            - DIR: UniDesign directory
            - BIN: UniDesign binary path
            - LIBRARY: UniDesign library path
            - WREAD: UniDesign wread path
        """
        return {
            "DIR": self.DIR,
            "BIN": self.BIN,
            "LIBRARY": self.LIBRARY,
            "WREAD": self.WREAD,
        }
