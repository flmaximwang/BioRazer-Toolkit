import logging, site, re, sys, os
from pathlib import Path

import pyrosetta
from pyrosetta.rosetta.core.select import get_residue_set_from_subset
from pyrosetta.rosetta.core.pose import PDBInfo, Pose


def print_pyrosetta_dir():
    site_packages_dir = site.getsitepackages()
    pyrosetta_path = Path(site_packages_dir[0]) / "pyrosetta"
    return str(pyrosetta_path)


def redict_to_file(
    target_file,
    logging_level=logging.INFO,
    format="[%(asctime)s] %(message)s",
    mode="a",
):
    """
    Execute this function before pyrosetta.init() to redirect the logger to a file.
    """

    logger = logging.getLogger("rosetta")  # Make the logger named rosetta
    logger.setLevel(logging_level)

    ###################### Constructing FileHandler ######################
    target_file = Path(target_file)
    file_handler = logging.FileHandler(target_file, mode=mode)

    class EscapeFormatter(logging.Formatter):
        def format(self, record):
            s = super().format(record)
            return re.sub(r"\x1b\[[0-9;]*m", "", s)  # 移除 ANSI 转义序列)

    file_handler.setFormatter(EscapeFormatter(format))
    file_handler.setLevel(logging_level)

    if logger.hasHandlers():
        logger.handlers.clear()
    logger.addHandler(file_handler)
    ######################################################################

    pyrosetta.logging_support.set_logging_sink()  # Make sure Pyrosetta uses logger "rosetta"
    return logger


def quiet_init(**kwargs):
    with open(os.devnull, "w") as devnull:
        sys.stdout = devnull
        sys.stderr = devnull
        pyrosetta.init(**kwargs)
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__


def generate_pose_from_residue_mask(pose, residue_mask):
    """
    Generate a new pose from the residue mask.
    All water molecules will be lost.
    PDBInfo of the new pose has not been correctlyly updated yet.
    """
    residue_set = get_residue_set_from_subset(residue_mask)
    new_pose = Pose()
    for residue_id in residue_set:
        new_pose.append_residue_by_jump(
            pose.residue(residue_id), len(new_pose.residues)
        )

    pdb_info = pose.pdb_info()
    new_pdb_info = PDBInfo(new_pose.size())

    for new_residue_index, residue_id in enumerate(residue_set, start=1):
        new_pdb_info.set_resinfo(
            new_residue_index, pdb_info.chain(residue_id), pdb_info.number(residue_id)
        )

    new_pose.pdb_info(new_pdb_info)

    return new_pose
