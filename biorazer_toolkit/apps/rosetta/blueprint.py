from dataclasses import dataclass, field
import os
import subprocess
from pathlib import Path
import pandas as pd
from . import config as rosetta_config


@dataclass
class Blueprint:
    data: pd.DataFrame = field(default_factory=lambda: Blueprint.default_data())

    @staticmethod
    def default_data() -> pd.DataFrame:
        data = pd.DataFrame(columns=["res_id_pos", "aa", "ss", "command"])
        col_type_dict = dict(
            res_id_pos=int,
            aa=str,
            ss=str,
            command=str,
        )
        for col, col_type in col_type_dict.items():
            data[col] = data[col].astype(col_type)
        return data

    @classmethod
    def from_pdb(cls, pdb, chain, script_dir: str | Path):
        bp_str = cls.pdb2str(pdb, chain, script_dir)
        return cls.from_str(bp_str)

    @staticmethod
    def pdb2str(pdb, chain, script_dir: str | Path):
        """
        Generate a Blueprint string from a PDB file with Rosetta's getBluePrintFromCoords.pl script.
        """
        try:
            bp_str = subprocess.check_output(
                [
                    f"{script_dir}/getBluePrintFromCoords.pl",
                    "-pdbfile",
                    pdb,
                    "--chain",
                    chain,
                ],
                text=True,
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Failed to convert PDB {pdb} to Blueprint string: {e}")

        return bp_str.strip()

    @staticmethod
    def pdb2bp(pdb, bp, chain):

        bp_path = Path(bp)
        if not bp_path.parent.exists():
            bp_path.parent.mkdir(parents=True)

        bp_str = Blueprint.pdb2str(pdb, chain)
        with open(bp, "w") as f:
            f.write(bp_str)
        return bp

    @staticmethod
    def from_str(bp_str):
        blueprint = Blueprint()
        data = blueprint.data

        for line in bp_str.strip().split("\n"):
            parts = line.strip().split()
            new_row = dict(
                res_id_pos=int(parts[0]),
                aa=parts[1],
                ss=parts[2],
                command=" ".join(parts[3:]) if len(parts) > 3 else "",
            )
            data.loc[len(data)] = new_row

        return blueprint

    @staticmethod
    def from_bp(bp):

        with open(bp, "r") as f:
            bp_str = "".join(f.readlines())
        return Blueprint.from_str(bp_str)

    def to_bp(self, bp_file):
        data = self.data
        with open(bp_file, "w") as f:
            for i in data.index:
                line = f"{data.loc[i, 'res_id_pos']} {data.loc[i, 'aa']} {data.loc[i, 'ss']} {data.loc[i, 'command']}\n"
                f.write(line)

    def get_res_is_pos_index(self, res_id_pos):
        if not isinstance(res_id_pos, int):
            raise ValueError(f"res_id_pos must be an integer, got {type(res_id_pos)}")
        if res_id_pos < 0:
            raise ValueError(f"res_id_pos must be non-negative, got {res_id_pos}")
        index_res_id_pos = (
            self.get_data().index[self.get_data()["res_id_pos"] == res_id_pos].tolist()
        )
        if len(index_res_id_pos) == 0:
            raise ValueError(f"No entry found for res_id_pos {res_id_pos}.")
        if len(index_res_id_pos) > 1:
            raise ValueError(f"Multiple entries found for res_id_pos {res_id_pos}.")
        return index_res_id_pos[0]

    def set_ss(self, res_id_pos, ss):
        index_res_id_pos = self.get_res_is_pos_index(res_id_pos)
        bp_data = self.get_data()
        bp_data.loc[index_res_id_pos, "ss"] = ss

    def set_command(self, res_id_pos, command):
        index_res_id_pos = self.get_res_is_pos_index(res_id_pos)
        bp_data = self.get_data()
        bp_data.loc[index_res_id_pos, "command"] = command

    def set_aa(self, res_id_pos, aa):
        index_res_id_pos = self.get_res_is_pos_index(res_id_pos)
        bp_data = self.get_data()
        if len(aa) != 1:
            raise ValueError(f"AA must be a single character, got {aa}.")
        bp_data.loc[index_res_id_pos, "aa"] = aa

    def insert_seq(self, seq, res_id_pos, ss="L", extend=[1, 1]):
        bp_data = self.get_data()
        index_res_id_pos = self.get_res_is_pos_index(res_id_pos)

        for index_tmp in range(
            index_res_id_pos + 1 - extend[0], index_res_id_pos + 1 + extend[1]
        ):
            aa_tmp = bp_data.loc[index_tmp, "aa"]
            bp_data.loc[index_tmp, "ss"] = ss
            bp_data.loc[index_tmp, "command"] = f"PIKAA {aa_tmp}"

        index_left = list(bp_data.index[: index_res_id_pos + 1])
        index_right = list(bp_data.index[index_res_id_pos + 1 :] + len(seq))
        bp_data.index = index_left + index_right
        for i, aa in enumerate(seq, start=1):
            new_row = dict(
                res_id_pos=0,
                aa="X",
                ss=ss,
                command=f"PIKAA {aa}",
            )
            bp_data.loc[index_res_id_pos + i] = new_row
        bp_data.sort_index(inplace=True)
