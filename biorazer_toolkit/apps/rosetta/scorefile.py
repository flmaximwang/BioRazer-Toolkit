from dataclasses import dataclass
from pathlib import Path
from io import StringIO
import pandas as pd


@dataclass
class Scorefile:
    data: pd.DataFrame

    @classmethod
    def from_sc(cls, sc_file: str | Path):
        sc_file_str = ""
        with open(sc_file, "r") as f:
            for line in f:
                if line.startswith("SCORE:"):
                    sc_file_str += line
        df = pd.read_csv(StringIO(sc_file_str), comment="#", sep=r"\s+")
        df = df.iloc[:, 1:]  # drop the first column which is just 'SCORE:'
        return cls(data=df)

    def to_csv(self):
        self.data.to_csv(index=False)
