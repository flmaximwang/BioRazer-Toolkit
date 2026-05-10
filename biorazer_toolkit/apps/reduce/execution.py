import tempfile
from pathlib import Path
from biotite.structure import AtomArray
from biorazer.structure.io.protein import PDB2PDB
from ...utils.structure_file import atom_array_as_temp_file
from ..basic import App


class Reduce(App):

    def __init__(self, app_dir=None, app_bin="reduce", log_name="reduce"):
        super().__init__(app_dir, app_bin, log_name=log_name)

    def run_with_structure(
        self,
        atom_array: AtomArray,
        *args,
        input_file_format="pdb",
        cwd=".",
        get_output=True,
        verbose=True,
        mode="pty",
        **kwargs,
    ):
        with atom_array_as_temp_file(
            atom_array,
            temp_file_format=input_file_format,
        ) as input_file:
            return self.run(
                *args,
                str(input_file),
                cwd=cwd,
                get_output=get_output,
                verbose=verbose,
                mode=mode,
                **kwargs,
            )


class ReduceFile(Reduce):
    """
    Run reduce in terms of input and output files.
    """

    def run(self, *args, cwd=".", output_file=None, **kwargs):
        result = super().run(*args, cwd=cwd, get_output=True, **kwargs)
        if output_file:
            with open(output_file, "w") as f:
                f.write(result)


class ReduceArray(Reduce):
    """
    Run reduce in terms of Biotite AtomArray.
    """

    def run(self, atom_array: AtomArray, *args, **kwargs):
        with tempfile.TemporaryDirectory(dir=tempfile.gettempdir()) as temp_dir:
            output_file = Path(temp_dir) / "reduce_output.pdb"
            with atom_array_as_temp_file(atom_array, dir=temp_dir) as input_file:
                result = super().run(*args, str(input_file), get_output=True, **kwargs)
            with open(output_file, "w") as f:
                f.write(result)
            reduced_array = PDB2PDB(input_file=output_file).read()
            return reduced_array

    def run_with_structure(
        self,
        atom_array: AtomArray,
        *args,
        input_file_format="pdb",
        **kwargs,
    ):
        normalized_format = input_file_format.lower().lstrip(".")
        if normalized_format != "pdb":
            raise ValueError(
                "ReduceArray only supports PDB temporary input because reduce output parsing expects PDB text."
            )
        return self.run(atom_array, *args, **kwargs)
