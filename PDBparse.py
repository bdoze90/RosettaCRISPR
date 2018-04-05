"""Has all the necessary functions for parsing PDB files"""


class PDB:

    def __init__(self, file):
        self.load_file(file)

    def load_file(self, file):
        f = open(file)
        for line in f:
            self.get_line_type(line)
            self.parse_line(line)

    def get_line_type(self, line):
        if line.find("ATOM") != -1:
            return "ATOM"
        elif line.find("HEADER") != -1:
            return "HEADER"
        elif line.find("TITLE") != -1:
            return "TITLE"
        elif line.find("REMARK") != -1:
            return "REMARK"
        elif line.find("COMPND") != -1:
            return "COMPND"
        elif line.find("SOURCE") != -1:
            return "SOURCE"
        elif line.find("KEYWDS") != -1:
            return "KEYWDS"
        elif line.find("EXPDATA") != -1:
            return "EXPDATA"
        elif line.find("AUTHOR") != -1:
            return "AUTHOR"
        elif line.find("REVDAT") != -1:
            return "REVDAT"
        elif line.find("JRNL") != -1:
            return "JRNL"
        elif line.find("DBREF") != -1:
            return "DBREF"
        elif line.find("SEQADV") != -1:
            return "SEQADV"
        elif line.find("SEQRES") != -1:
            return "SEQRES"
        elif line.find("") != -1:
            return ""

    def parse_line(self, l):


def ATOM_parse(l):
    name = l[:6]
    serial = int(l[6:11])
    a_name = l[12:16]
    altloc = l[16]
    resName = l[17:20]
    chainID = l[21]
    resSeq = l[22:26]
    iCode = l[26]
    x = l[30:38]
    y = l[38:46]
    z = l[46:54]
    occupancy = l[54:60]
    tempFactor =
    element =
    charge =

def