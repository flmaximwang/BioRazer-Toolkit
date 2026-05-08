import numpy as np
from biotite.structure import AtomArray
from biotite.structure import hbond
import hydride
from pymol import cmd
from biorazer.structure.analysis.static.select import (
    mask_interface_atoms,
    mask_buried_unsat_hbond_atoms,
)
from biorazer.structure.analysis.static.check import is_hydrided
from biorazer.structure.analysis.static.report import report_hbonds as br_report_hbonds
from biorazer.structure.io import STRUCT2PDB
from .io import load_atom_array

__all__ = [
    "load_atom_array",
    "report_interface_residues",
    "report_hbond",
    "report_interface_buried_unsat_hbond",
]


def report_interface_residues(
    atom_array: AtomArray,
    selection1,
    selection2,
    distance_cutoff=3.5,
    model_name="default",
):
    interface_atom_mask_1, interface_atom_mask_2, _ = mask_interface_atoms(
        atom_array,
        selection1=selection1,
        selection2=selection2,
        distance_cutoff=distance_cutoff,
    )

    interfaces = [atom_array[interface_atom_mask_1], atom_array[interface_atom_mask_2]]
    interface_residues = [[], []]
    for i in range(2):
        for atom in interfaces[i]:
            identifier = (atom.chain_id, atom.res_id)
            if not identifier in interface_residues[i]:
                interface_residues[i].append(identifier)

    for i in range(2):
        cmd.select(f"{model_name}_interface_tmp_{i+1}", "not all")
        for chain_id, res_id in interface_residues[i]:
            cmd.select(
                f"{model_name}_interface_tmp_{i+1}",
                f"/{model_name}//{chain_id}/{res_id}/ or {model_name}_interface_tmp_{i+1}",
            )
        cmd.create(
            f"{model_name}_interface_sphere_{i+1}", f"{model_name}_interface_tmp_{i+1}"
        )
        cmd.hide("everything", f"{model_name}_interface_sphere_{i+1}")
        cmd.show("sphere", f"{model_name}_interface_sphere_{i+1} and not c. C+N+O+CA")
        cmd.create(
            f"{model_name}_interface_sticks_{i+1}", f"{model_name}_interface_tmp_{i+1}"
        )
        cmd.hide("everything", f"{model_name}_interface_sticks_{i+1} and not n. C+N+O")
        cmd.show("sticks", f"{model_name}_interface_sticks_{i+1}")
        cmd.delete(f"{model_name}_interface_tmp_{i+1}")


def report_hbonds(
    atom_array: AtomArray,
    selection1,
    selection2,
    model_name="default",
    cutoff_dist=2.5,
    cutoff_angle=120,
    donor_elements=("O", "N", "S"),
    acceptor_elements=("O", "N", "S"),
    periodic=False,
):
    """

    Parameters
    ----------
    format: str
        - pymol: print PyMOL commands to visualize the hydrogen bonds
        - list: return a list of tuples (donor, hydrogen, acceptor)
        - text: print hydrogen bonds in text format
    """

    hbonds = br_report_hbonds(
        atom_array,
        selection1=selection1,
        selection2=selection2,
        cutoff_dist=cutoff_dist,
        cutoff_angle=cutoff_angle,
        donor_elements=donor_elements,
        acceptor_elements=acceptor_elements,
        periodic=periodic,
        fmt="list",
    )

    for i, (donor, hydrogen, acceptor) in enumerate(hbonds, start=1):
        cmd.distance(
            f"hbond_{i}",
            f"{model_name}//{donor.chain_id}/{donor.res_id}/{donor.atom_name}",
            f"{model_name}//{acceptor.chain_id}/{acceptor.res_id}/{acceptor.atom_name}",
        )

    return hbonds


def report_buried_unsat_hbonds(
    atom_array: AtomArray,
    selection1,
    selection2,
    model_name="default",
    sasa_kwargs=dict(sasa_cutoff=0.5, probe_radius=1.4),
    interface_kwargs=dict(distance_cutoff=3.5),
    hbond_kwargs=dict(),
):
    """
    Report unsatisfied hydrogen bonds in the interface between two selections of atoms.

    Parameters
    ----------
    atom_array : bio_struct.AtomArray
        The biotite structure array containing the atoms.
    selection1 : np.ndarray
        The atom mask for the 1st selection.
    selection2 : np.ndarray
        The atom mask for the 2nd selection.
    format : str, optional
        The output format, by default "pymol".
    """

    unsat_hbond_mask_1, unsat_hbond_mask_2, _ = mask_buried_unsat_hbond_atoms(
        atom_array,
        selection1=selection1,
        selection2=selection2,
        sasa_kwargs=sasa_kwargs,
        interface_kwargs=interface_kwargs,
        hbond_kwargs=hbond_kwargs,
    )

    buried_unsat_hbonds = [
        atom_array[unsat_hbond_mask_1],
        atom_array[unsat_hbond_mask_2],
    ]
    for i in range(2):
        cmd.select(f"{model_name}_unsat_hbond_tmp_{i+1}", "not all")
        for atom in buried_unsat_hbonds[i]:
            chain_id = atom.chain_id
            res_id = atom.res_id
            atom_name = atom.atom_name
            cmd.select(
                f"{model_name}_unsat_hbond_tmp_{i+1}",
                f"/{model_name}//{chain_id}/{res_id}/{atom_name} or {model_name}_unsat_hbond_tmp_{i+1}",
            )
        cmd.create(
            f"{model_name}_unsat_hbond_sphere_{i+1}",
            f"{model_name}_unsat_hbond_tmp_{i+1}",
        )
        cmd.hide("everything", f"{model_name}_unsat_hbond_sphere_{i+1}")
        cmd.show("sphere", f"{model_name}_unsat_hbond_sphere_{i+1}")
        cmd.delete(f"{model_name}_unsat_hbond_tmp_{i+1}")
    return buried_unsat_hbonds
