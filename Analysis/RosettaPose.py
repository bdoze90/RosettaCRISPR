"""File for manipulating Pose Data."""



class PoseVector:

    def __init__(self):
        self.type = str()
        self.algorithm = str()

        # Vectors of scores associated with each type of scoring algorithm and type of molecule
        self.dna_scrs = ["fa_atr", "fa_rep", "fa_sol", "fa"]


    def re_weight(self):
        poo =1