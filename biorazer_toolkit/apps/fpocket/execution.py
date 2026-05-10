import subprocess, re, shutil, warnings, time
from pathlib import Path
import pandas as pd
from sklearn.cluster import DBSCAN
from biotite.structure import AtomArray
from biorazer.structure.io.protein import PDB2STRUCT, CIF2STRUCT
from biorazer.structure.util.geometry.surface import fibonacci_surface_grid
from biorazer.structure.util.dictionary.radius import vdw_radii
from biorazer.structure.analysis.static.select import mask_atoms_within_distance
from ...utils.structure_file import call_with_structure_file

warnings.filterwarnings("ignore", module="biotite.structure.io.pdb")


def collect(
    output_dir: str | Path, input_suffix=".pdb", probe_radius=1.5, n_points=512
):
    output_dir_path = Path(output_dir)
    pdb_cif_path = output_dir_path.parent / (
        output_dir_path.stem.removesuffix("_out") + input_suffix
    )
    pocket_info_txt = output_dir_path / f"{pdb_cif_path.stem}_info.txt"
    keys = [
        "Pocket ID",
        "Score",
        "Druggability Score",
        "Number of Alpha Spheres",
        "Total SASA",
        "Polar SASA",
        "Apolar SASA",
        "Volume",
        "Mean local hydrophobic density",
        "Mean alpha sphere radius",
        "Mean alp. sph. solvent access",
        "Apolar alpha sphere proportion",
        "Hydrophobicity score",
        "Volume score",
        "Polarity score",
        "Charge score",
        "Proportion of polar atoms",
        "Alpha sphere density",
        "Cent. of mass - Alpha Sphere max dist",
        "Flexibility",
        "Type",
    ]

    all_data = []
    row_data = {}
    with open(pocket_info_txt, "r") as f:
        for line in f:
            if line.startswith("Pocket"):
                if len(row_data) > 0:
                    all_data.append(row_data)
                row_data = {}
                matched = re.match(r"Pocket (\d+).*", line)
                pocket_id = int(matched.group(1))
                row_data["Pocket ID"] = pocket_id
                continue
            if len(line.strip()) == 0:
                continue
            col_name = line.split(":")[0].strip()
            value = float(line.split(":")[1].strip())
            row_data[col_name] = value
    pocket_df = pd.DataFrame(data=all_data, columns=keys)

    pocket_pdbs = output_dir_path.glob("pockets/pocket*_vert.pqr")
    pocket_pdbs = list(pocket_pdbs)
    pocket_pdbs.sort(
        key=lambda x: int(re.match(r"pocket(\d+)_vert\.pqr", x.name).group(1))
    )
    if pdb_cif_path.suffix.lower() == ".cif":
        structure = CIF2STRUCT(pdb_cif_path, "").read()
    elif pdb_cif_path.suffix.lower() == ".pdb":
        structure = PDB2STRUCT(pdb_cif_path, "").read()
    else:
        raise ValueError(f"Unsupported file format: {pdb_cif_path.suffix}")
    for i, pdb in enumerate(pocket_pdbs):
        pocket = PDB2STRUCT(pdb, "").read()
        neighbor_mask = mask_atoms_within_distance(pocket, structure, distance=5.0)
        neighbors = structure[neighbor_mask]
        fibonacci_grid = fibonacci_surface_grid(
            neighbors.coord,
            vdw_radii=vdw_radii(neighbors.element),
            probe_radius=probe_radius,
            n_points=n_points,
        )
        cluster_res = DBSCAN(eps=1.0, min_samples=10).fit(fibonacci_grid)
        surface_indices = set(
            cluster_res.labels_[cluster_res.labels_ >= 0]
        )  # minus indices represent samples that are less than min_samples, which should be noise
        if len(surface_indices) == 0:
            raise RuntimeError(
                "No surface points found. Try increasing the distance threshold in mask_atom_within_distance."
            )
        elif len(surface_indices) == 1:
            pocket_df.loc[i, "Type"] = "Pocket"
        else:
            pocket_df.loc[i, "Type"] = "Cavity"

    return pocket_df


def run(pdb_cif, *args, n_points=512, probe_radius=1.5, **kwargs):
    """
    Python wrapper for fpocket command line tool.

    Parameters
    ----------
    pdb_cif: str
        Path to the PDB or CIF file to analyze.
    n_points: int
        Number of points on the Fibonacci sphere for surface sampling. Default is 512.
    probe_radius: float
        Radius of the probe used for surface sampling. Default is 1.5 Å.
    args: list[str]
        arguments without input values for command line
    kwargs: dict
        arguments with input values for command line

    CLI Usage
    ---------
    Mandatory parameters",
        fpocket -f --file pdb or cif file
        [ fpocket -F --fileList fileList ]


    Optional output parameters
        -x --calculate_interaction_grids      : Specify this flag if you want fpocket to
                                                calculate VdW and Coulomb grids for each pocket
        -d --pocket_descr_stdout              : Put this flag if you want to write fpocket
                                                descriptors to the standard output


    Optional input parameters
        -l --model_number (int)               : Number of Model to analyze.
        -y --topology_file (string)           : File name of a topology file (Amber prmtop).
        -r --custom_ligand (string)           : String specifying a ligand like:
                                                residuenumber:residuename:chain_code (ie. 1224:PU8:A).
        -P --custom_pocket (string)           : String specifying a pocket like:
                                                residuenumber1:insertion_code1('-' if empty):chain_code1.residuenumber2:insertion_code2:chain_code2 (ie. 138:-:A.139:-:A).
        -u --min_n_explicit_pocket (int)      : If explicit pocket provided, minimum number
                                                of atoms of an alpha sphere that have to be in the selected pocket.
        -a --chain_as_ligand (char)           : Character specifying a chain as a ligand


    Optional pocket detection parameters (default parameters)
        -m --min_alpha_size (float)           : Minimum radius of an alpha-sphere.    (3.4)
        -M --max_alpha_size (float)           : Maximum radius of an alpha-sphere.    (6.2)
        -D --clustering_distance (float)      : Distance threshold for clustering algorithm   (2.4)
        -C --clustering_method (char)         : Specify the clustering method wanted for
                                                grouping voronoi vertices together (s)
                                                s : single linkage clustering
                                                m : complete linkage clustering
                                                a : average linkage clustering
                                                c : centroid linkage clustering
        -e --clustering_measure (char)        : Specify the distance measure for clustering   (e)
                                                e : euclidean distance
                                                b : Manhattan distance
        -i --min_spheres_per_pocket (int)     : Minimum number of a-sphere per pocket.        (15)
        -p --ratio_apol_spheres_pocket (float): Minimum proportion of apolar sphere in
                                                a pocket (remove otherwise) (0.0)
        -A --number_apol_asph_pocket (int)    : Minimum number of apolar neighbor for
                                                an a-sphere to be considered as apolar.   (3)
        -v --iterations_volume_mc (integer)   : Number of Monte-Carlo iteration for the
                                                calculation of each pocket volume.(300)
        -c --drop_chains (char)               : Name of the chains to be deleted before pocket detection,
                                                able to delete up to (20) chains (ie : -c A,B,E)
        -k --keep_chains (char)               : Name of the chains to be kept before pocket detection,
                                                able to keep up to (20) chains (ie : -k A,B,C,E)
        -a --chain_as_ligand (char)           : consider this chain as a ligand explicitly (i.e. -a D)
        -w --write_mode (char)                : Writing mode to be used after pocket detection,
                                                d -> default (same format outpout as input)
                                                b or both -> both pdb and mmcif | p or pdb ->pdb | m or cif or mmcif-> mmcif
    For more information: http://fpocket.sourceforge.net
    """

    pdb_cif_path = Path(pdb_cif)
    output_dir_path = pdb_cif_path.parent / (pdb_cif_path.stem + "_out")
    if output_dir_path.exists():
        shutil.rmtree(output_dir_path)

    # fpocket_start = time.time()
    cmd_args = ["fpocket", "-f", pdb_cif]
    cmd_args.extend(args)
    for key, value in kwargs.items():
        cmd_args.append(str(key))
        cmd_args.append(str(value))
    p = subprocess.run(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.stdout, p.stderr
    if p.returncode != 0:
        raise RuntimeError(f"fpocket failed with error: {stderr.decode('utf-8')}")
    # fpocket_end = time.time()
    # print(f"fpocket finished in {fpocket_end - fpocket_start:.2f} seconds")

    # collection_start = time.time()
    pocket_df = collect(
        output_dir_path,
        input_suffix=pdb_cif_path.suffix,
        probe_radius=probe_radius,
        n_points=n_points,
    )
    pocket_df.to_csv(output_dir_path / f"{pdb_cif_path.stem}_pockets.csv", index=False)
    return pocket_df


def run_with_structure(
    atom_array: AtomArray,
    *args,
    input_file_format: str = "pdb",
    n_points=512,
    probe_radius=1.5,
    **kwargs,
):
    return call_with_structure_file(
        atom_array,
        run,
        *args,
        temp_file_format=input_file_format,
        n_points=n_points,
        probe_radius=probe_radius,
        **kwargs,
    )
