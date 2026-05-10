from dataclasses import dataclass, field
from pathlib import Path

SUGGESTED_COMMANDS = {
    "ALLAA": "allow all 20 amino acids INCLUDING cysteine (same as ALLAAwc)",
    "ALLAAwc": "allow all 20 amino acids ( default )",
    "ALLAAxc": "allow all amino acids except cysteine",
    "POLAR": "allow only canonical polar amino acids (DEHKNQRST)",
    "APOLAR": "allow only canonical non polar amino acids (ACFGILMPVWY)",
    "PROPERTY <property>": "disallow any residue type that lacks the given property",
    "NOTAA <list of AAs>": "disallow only the specified amino acids (Use one letter codes, undelimited like ACFYRT. For NCAAs, use X[<full name>].)",
    "PIKAA <list of AAs>": "allow only the specified amino acids (Use one letter codes, undelimited like ACFYRT. For NCAAs, use X[<full name>].)",
    "NATAA": "allow only the native amino acid (NATive Amino Acid) - repack without design",
    "NATRO": "preserve the input rotamer ( do not pack at all) (NATive ROtamer)",
    "APOLAR NOTAA C": "allow only canonical non polar amino acids except cysteine (AFGILMPVWY)",
    "EX <i>": "allow extra rotamers level i (i=1,2,3... )",
}

COMMAND_SPLITORS = [
    "ALLAA",
    "ALLAAwc",
    "ALLAAxc",
    "POLAR",
    "APOLAR",
    "NATAA",
    "NATRO",
    "PROPERTY",
    "NOTAA",
    "PIKAA",
    "EX",
]

SUGGESTED_RESIDUE_IDENTIFIERS = {
    "* A": "All residues in chain A",
    "10 _": "Residue 10 in the unlabeled chain",
    "40 B": "Residue 40 in chain B",
    "40A Q": "Residue 40, insertion code A, on chain Q",
    "3C - 20 A": "From residues 3, insertion code C to residue 20 in chain A",
}


@dataclass
class HeaderLine:
    commands: list[str]

    def __str__(self) -> str:
        return " ".join(self.commands)


@dataclass
class BodyLine:
    residue_identifier: str
    commands: list[str]

    def __str__(self) -> str:
        return f"{self.residue_identifier} " + " ".join(self.commands)


@dataclass
class Resfile:
    header_lines: list[HeaderLine] = field(default_factory=list)
    body_lines: list[BodyLine] = field(default_factory=list)

    @classmethod
    def from_txt(cls, txt_file: str | Path):
        header_lines = []
        body_lines = []
        with open(txt_file, "r") as f:
            lines = f.readlines()
            is_header = True
            for line in lines:
                line = line.strip()
                if line == "" or line.startswith("#"):
                    continue
                if is_header:
                    if line.upper() == "START":
                        is_header = False
                    else:
                        commands = line.split()
                        header_lines.append(HeaderLine(commands))
                else:
                    parts = line.split()
                    not_command = True
                    residue_identifiers = []
                    commands = []
                    for part in parts:
                        if not_command:
                            if part in COMMAND_SPLITORS:
                                not_command = False
                                commands.append(part)
                            else:
                                residue_identifiers.append(part)
                        else:
                            commands.append(part)
                    residue_identifier = " ".join(residue_identifiers)
                    command = " ".join(commands)
                    body_lines.append(BodyLine(residue_identifier, [command]))
        return cls(header_lines, body_lines)

    def __str__(self) -> str:
        lines = []
        for header_line in self.header_lines:
            lines.append(str(header_line))
        lines.append("START")
        for body_line in self.body_lines:
            line = str(body_line)
            lines.append(line)
        return "\n".join(lines) + "\n"
