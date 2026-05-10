from ..basic import App
from pathlib import Path
from ...utils.structure_file import atom_array_as_temp_file


class MASTERApp(App):
    """
    Do not provide app_bin
    """

    def run_createPDS(self, *args, get_output=False, verbose=True, **kwargs):
        self.bin = self.dir / "bin" / "createPDS"
        return self.run(*args, get_output=get_output, verbose=verbose, **kwargs)

    def run_parsePDS(self, *args, get_output=False, verbose=True, **kwargs):
        self.bin = self.dir / "bin" / "parsePDS"
        return self.run(*args, get_output=get_output, verbose=verbose, **kwargs)

    def run_master(self, *args, get_output=False, verbose=True, **kwargs):
        self.bin = self.dir / "bin" / "master"
        return self.run(*args, get_output=get_output, verbose=verbose, **kwargs)

    def run_with_structure(
        self,
        atom_array,
        *args,
        input_file_format="pdb",
        runner="run_master",
        input_file_flag=None,
        get_output=False,
        verbose=True,
        **kwargs,
    ):
        run_method = getattr(self, runner, None)
        if run_method is None or not callable(run_method):
            raise ValueError(f"Unknown MASTER runner: {runner}")

        with atom_array_as_temp_file(
            atom_array,
            temp_file_format=input_file_format,
        ) as input_file:
            run_args = list(args)
            if input_file_flag is None:
                run_args.append(str(input_file))
            else:
                run_args.extend([str(input_file_flag), str(input_file)])
            return run_method(
                *run_args,
                get_output=get_output,
                verbose=verbose,
                **kwargs,
            )
