import MDAnalysis as mda


def slice_trajectory(universe_args, slice, sliced_xtc):

    u = mda.Universe(*universe_args)
    all = u.select_atoms("all")
    with mda.Writer(sliced_xtc, all.n_atoms) as W:
        for ts in u.trajectory[slice]:
            W.write(all)
            u.trajectory.next()
