import os, subprocess, shutil
from pathlib import Path
from biotite.structure import AtomArray
from ...utils.structure_file import call_with_structure_file
from .config import UniDesignConfig

config = None


def init(app_dir):
    global config
    config = UniDesignConfig(None, None)
    config.set(app_dir)
    config.check()


def run(input_pdb, output_dir, cmd_kwargs, resfile="RESFILE.txt", log_name="UniDesign"):
    global config
    assert (
        config is not None
    ), "UniDesignConfig is not initialized. Call init(app_dir) first."
    config_dict = config.get()
    unidesign_bin = config_dict["BIN"]
    unidesign_library = config_dict["LIBRARY"]
    unidesign_wread = config_dict["WREAD"]

    input_pdb_path = Path(input_pdb).resolve()
    output_dir_path = Path(output_dir)
    if not output_dir_path.exists():
        output_dir_path.mkdir(parents=True, exist_ok=True)
    for file_name, file in zip(
        ["UniDesign", "library", "wread"],
        [unidesign_bin, unidesign_library, unidesign_wread],
    ):
        if (output_dir_path / file_name).is_symlink() and (
            output_dir_path / file_name
        ).exists():
            continue
        else:
            (output_dir_path / file_name).symlink_to(Path(file).resolve())

    cwd_ori = os.getcwd()
    os.chdir(output_dir)
    shutil.copyfile(input_pdb_path, input_pdb_path.name)
    unidesign_command_list = ["./UniDesign"]
    unidesign_command_list.append(f"--pdb={input_pdb_path.name}")
    unidesign_command_list.append(f"--resfile={resfile}")
    for key, value in cmd_kwargs.items():
        if isinstance(value, bool):
            if value:
                unidesign_command_list.append(f"--{key}")
        else:
            unidesign_command_list.append(f"--{key}={value}")
    with open(f"{log_name}.log", "w") as log_file:
        result = subprocess.run(
            unidesign_command_list, stdout=log_file, stderr=log_file
        )
    if result.returncode != 0:
        raise RuntimeError(f"UniDesign failed. See {log_name}.log for details.")
    os.chdir(cwd_ori)


def run_with_structure(
    atom_array: AtomArray,
    output_dir,
    cmd_kwargs,
    resfile="RESFILE.txt",
    log_name="UniDesign",
    input_file_format: str = "pdb",
):
    """
    Run UniDesign from an in-memory structure by materializing a temporary file.

    Parameters
    ----------
    atom_array:
        Input structure in biotite AtomArray format.
    input_file_format:
        Temporary file format used as app input. Default is "pdb".
    """
    return call_with_structure_file(
        atom_array,
        run,
        output_dir,
        cmd_kwargs,
        resfile=resfile,
        log_name=log_name,
        temp_file_format=input_file_format,
    )


# Backward compatibility alias
run_structure = run_with_structure
