import pyrosetta
from xml.etree import ElementTree as ET
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects

pyrosetta.init()
pose = pyrosetta.pose_from_pdb(
    "../data/validative/pipeline-1.3.0/alpha-1/selection-1/00040_01.pdb"
)


def create_ddg_filter():
    root = ET.Element("Ddg")
    root.set("name", "ddg_filter")
    root.set("threshold", "-15")  # 阈值
    root.set("jump", "1")  # 分离第一对链间连接
    root.set("repeats", "5")  # 重复次数
    root.set("repack", "true")  # 允许重新包装
    root.set("translate_by", "100")  # 分离距离
    root.set("confidence", "0")  # 置信度
    xml_str = ET.tostring(root, encoding="unicode")
    filter = XmlObjects.static_get_filter(xml_str)
    return filter


ddg_filter = create_ddg_filter()
ddg = ddg_filter.report_sm(pose)
