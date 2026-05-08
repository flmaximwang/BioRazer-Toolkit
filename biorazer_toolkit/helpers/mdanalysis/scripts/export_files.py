import MDAnalysis as mda


def export_xtc(u: mda.Universe, output_xtc):
    with mda.Writer(output_xtc, n_atoms=u.atoms.n_atoms) as W:
        for ts in u.trajectory:
            W.write(u.atoms)
