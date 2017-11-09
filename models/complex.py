import numpy as np
from scipy.spatial.distance import cdist

class Complex():
    """
    The class Complex stores two proteins(dockerasmus.pdb.Protein) treated in the program.
    """

    def __init__(self,rec,lig):
        self.rec = rec
        self.lig = lig

    def min_ca(self):
        """
        Return the minimum distance between receptor's and ligand's carbon alpha
        """
        # TODO : Use ann when matrix size is large
        return np.min(cdist(self.rec_ca_pos,self.lig_ca_pos)
