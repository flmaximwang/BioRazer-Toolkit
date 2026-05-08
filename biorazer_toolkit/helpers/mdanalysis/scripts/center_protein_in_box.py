import numpy as np
from tqdm import tqdm
import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.transformations import wrap, center_in_box, unwrap


def center_monomer_protein_in_box_in_memory(u: mda.Universe):

    # Check whether the universe is 'in_memory'
    if not isinstance(u.trajectory, MemoryReader):
        raise ValueError("The universe must not be 'in_memory'.")

    protein = u.select_atoms("protein")
    not_protein = u.select_atoms("not protein")
    for ts in tqdm(u.trajectory):
        protein.unwrap(compound="fragments")
        protein_center = protein.center_of_mass(wrap=True)
        dim = ts.triclinic_dimensions
        box_center = np.sum(dim, axis=0) / 2
        u.atoms.translate(box_center - protein_center)
        not_protein.wrap(compound="residues")


def center_monomer_protein_in_box(u: mda.Universe):

    protein = u.select_atoms("protein")
    not_protein = u.select_atoms("not protein")
    transforms = [
        unwrap(protein),
        center_in_box(protein, wrap=True),
        wrap(not_protein),
    ]
    u.trajectory.add_transformations(*transforms)


def center_multimer_protein_in_box(u: mda.Universe, center_chainid: str):

    protein = u.select_atoms("protein")
    center_chain = u.select_atoms(f"protein and chainid {center_chainid}")
    other_chains = u.select_atoms(f"protein and not chainid {center_chainid}")
    not_protein = u.select_atoms("not protein")
    transforms = [
        unwrap(protein),