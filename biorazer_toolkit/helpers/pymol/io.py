import shutil
from biotite.structure import AtomArray
from pymol import cmd
from biorazer.structure.io import STRUCT2PDB


def load_atom_array(atom_array: AtomArray, model_name):
    tmp_pdb = "tmp.pdb"
    STRUCT2PDB("", tmp_pdb).write(atom_array)
    cmd.load(tmp_pdb, model_name)
    shutil.rmtree(tmp_pdb)
