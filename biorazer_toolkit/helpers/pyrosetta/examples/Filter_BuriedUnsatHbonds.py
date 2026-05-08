import pyrosetta

setting_list = [
    "-ignore_unrecognized_res",
    "-ignore_zero_occupancy",
    "-holes:dalphaball ~/Applications/bin/DAlphaBall.gcc",
    "-corrections::beta_nov16 true",
]
pyrosetta.init(" ".join(setting_list))
import xml.etree.ElementTree as ET
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects


def create_buried_unsat_hbonds_filter():
    # 创建 1 个标签为 BuriedUnsatHbonds 的 xml 对象
    root = ET.Element("BuriedUnsatHbonds")
    # 设置 root 的 report_all_heavy_atom_unsats 属性为 true
    root.set("report_all_heavy_atom_unsats", "true")
    # root.set("scorefxn", "scorefxn")
    root.set("ignore_surface_res", "false")
    root.set("use_ddG_style", "true")
    root.set("dalphaball_sasa", "1")
    root.set("probe_radius", "1.1")
    root.set("burial_cutoff_apo", "0.2")
    root.set("confidence", "0")
    # 获取 root 的字符串表示
    ET.tostring(root)
    filter = XmlObjects.static_get_filter(ET.tostring(root))
    return filter


pose = pyrosetta.pose_from_pdb(
    "data/AF3_Repredictions/Complexes/bbf_14_b1/fold_bbf_14_b1_model_0.pdb"
)
filter = create_buried_unsat_hbonds_filter()
interface_delta_unsat_hbonds = filter.report_sm(
    pose
)  # 可以通过这个命令获取界面上的氢键数量
interface_delta_unsat_hbonds
