import subprocess, tempfile
from pathlib import Path
import numpy as np
import biotite.structure as bio_struct
from biorazer.structure.io.protein import PDB2PDB
from ..basic import App


class PulchraFixBB(App):

    def __init__(self, app_dir=None, app_bin="pulchra", logger=None):
        super().__init__(app_dir, app_bin, logger=logger)

    def run(self, atom_array: bio_struct.AtomArray, cwd=None) -> bio_struct.AtomArray:
        pulchra_exe = self.bin
        cache_dir = tempfile.gettempdir()
        with tempfile.TemporaryDirectory(dir=cache_dir) as temp_dir:
            chain_ids = list(set(atom_array.chain_id))
            chain_ids.sort()
            fixed_arrays = []
            for chain_id in chain_ids:
                mask = atom_array.chain_id == chain_id
                masked_array = atom_array[mask]
                input_file = Path(temp_dir) / "input.pdb"
                output_file = Path(temp_dir) / "input.rebuilt.pdb"
                fixed_file = Path(temp_dir) / "fixed.pdb"
                PDB2PDB(output_file=input_file).write(masked_array)
                cmd = [str(pulchra_exe), "input.pdb"]
                subprocess.run(cmd, check=True, cwd=temp_dir, stdout=subprocess.DEVNULL)
                with open(output_file) as f:
                    with open(fixed_file, "w") as fixed_f:
                        for line in f:
                            if line.startswith("ATOM"):
                                fixed_f.write(line[:-1] + "  1.00  0.00\n")
                            else:
                                fixed_f.write(line)
                fixed_array = PDB2PDB(input_file=fixed_file).read()
                fixed_array.chain_id = np.array([chain_id] * len(fixed_array))
                fixed_arrays.append(fixed_array)
        fixed_array = bio_struct.concatenate(fixed_arrays)
        return fixed_array

    def run_with_structure(
        self,
        atom_array: bio_struct.AtomArray,
        *args,
        input_file_format="pdb",
        cwd=None,
        **kwargs,
    ) -> bio_struct.AtomArray:
        normalized_format = input_file_format.lower().lstrip(".")
        if normalized_format != "pdb":
            raise ValueError(
                "PulchraFixBB currently only supports PDB temporary input."
            )
        return self.run(atom_array=atom_array, cwd=cwd)
