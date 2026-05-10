import re, sys
from ..basic import App
import os
from ...utils.structure_file import atom_array_as_temp_file


class RosettaApp(App):

    def _check_dir(self):
        # self.app_dir 必须如 */rosetta.source.*/main
        if not self.app_dir:
            raise ValueError("Rosetta directory is not set.")
        if not self.app_dir.is_dir():
            raise ValueError(f"Rosetta directory {self.app_dir} does not exist.")
        if not re.match(r".*rosetta\.source\..*[/\\]main", str(self.app_dir)):
            raise ValueError(
                f"Rosetta directory {self.app_dir} is not correct. It should be like */rosetta.source.*/main"
            )

    def run(self, *args, cwd=".", get_output=False, verbose=True, mode="pty", **kwargs):
        return super().run(
            *args, cwd=cwd, get_output=get_output, verbose=verbose, mode=mode, **kwargs
        )

    def run_with_structure(
        self,
        atom_array,
        *args,
        input_file_format="pdb",
        input_file_flag="-s",
        cwd=".",
        get_output=False,
        verbose=True,
        mode="pty",
        **kwargs,
    ):
        with atom_array_as_temp_file(
            atom_array,
            temp_file_format=input_file_format,
        ) as input_file:
            run_args = list(args)
            if input_file_flag is None:
                run_args.append(str(input_file))
            else:
                run_args.extend([str(input_file_flag), str(input_file)])
            return self.run(
                *run_args,
                cwd=cwd,
                get_output=get_output,
                verbose=verbose,
                mode=mode,
                **kwargs,
            )

    def find_app(self, include_keywords: list[str], exclude_keywords: list[str] = []):

        results = []
        for identifier in (self.app_dir / "source" / "bin").glob("*default*"):
            if all(keyword in identifier.name for keyword in include_keywords):
                if not any(keyword in identifier.name for keyword in exclude_keywords):
                    results.append(identifier)
        return results

    def use_app(self, app_keywords: list[str]):
        results = self.find_app(app_keywords)
        if not results:
            raise RuntimeError(
                f"Failed to find app with keywords {app_keywords} in {self.app_dir}"
            )
        if len(results) > 1:
            raise RuntimeError(
                f"Found multiple apps with keywords {app_keywords} in {self.app_dir}: {results}"
            )
        self.app_bin = results[0]
        self.logger.info(f"Using app {self.app_bin} for keywords {app_keywords}")

    def use_script(self):
        """使用  Rosetta Script 运行"""
        self._check_dir()
        if sys.platform == "darwin":
            self.app_bin = (
                self.app_dir / "source" / "bin" / "rosetta_scripts.macosclangrelease"
            )
        elif sys.platform == "linux":
            self.app_bin = (
                self.app_dir / "source" / "bin" / "rosetta_scripts.linuxgccrelease"
            )
        elif sys.platform == "win32":
            self.app_bin = (
                self.app_dir
                / "source"
                / "bin"
                / "rosetta_scripts.windowsgccrelease.exe"
            )
        else:
            raise RuntimeError(f"Unsupported platform: {sys.platform}")

    def find_tool(self, tool_keyword):

        results = []
        for parent, dirs, files in ((self.app_dir) / "tools").walk():
            if not ("tools" in str(parent).split(os.path.sep)):
                continue
            for file in files:
                if tool_keyword in file:
                    results.append(parent / file)
        return results
