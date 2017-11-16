import numpy as np
from scipy.spatial.distance import cdist
import quaternion as quat
from functions import self_rot_quat, rot_around_rec_quat

class Complex():
    """
    The class Complex stores two proteins(dockerasmus.pdb.Protein) treated in the program.
    lig_quat = lig - rec
    """

    def __init__(self,rec,lig):
        self.rec = rec
        self.lig = lig
        self.center_rec = np.mean(rec.atom_positions(),axis=0)
        self.center_rec_quat = quat.quaternion(*self.center_rec)
        self.center_lig = np.mean(lig.atom_positions(),axis=0)
        self.axis = self.center_lig - self.center_rec
        self.lig_quat = quat.as_quat_array(np.c_[np.zeros(lig.atom_positions().shape[0]),lig.atom_positions()])


    def ca_dist(self,lig):
        """
        Return the minimum distance between receptor's and ligand's carbon alpha
        """
        # TODO : Use ann when matrix size is large
        return cdist(self.rec.get_ca(),lig[self.lig.get_ca_ind()])

    def rotations(self,theta,phi,alpha,beta,gamma):
        """
        Return ligand atoms positions after rotation around receptor and self-rotation.
        Rotation around receptor will be treated first, then self-rotation
        The first rotation is caracterized by two parameters, theta and phi (Eular angles)
        The self-rotation is caracterized by thee parameters, alpha, beta, and gamma.

        Keyword arguments:
        theta -- ligand's rotation angle around axis defined by z-axis X rec -> lig
        phi -- ligand's rotation angle around z-axis
        alpha -- self-rotation's angle around z-axis X (z-axis X lig -> rec)
        beta -- self-rotation's angle around z-axis X lig -> rec
        gamma -- self-rotation's angle around vector ligand -> receptor
        """
        # TODO : Simplify
        transq = rot_around_rec_quat(self.axis,theta,gamma)
        tmp_axis_quat = transq*quat.quaternion(*self.axis)*np.conjugate(transq)
        tmp_axis = quat.as_float_array(tmp_axis_quat)[1:]
        rotateq = self_rot_quat(-tmp_axis,alpha,beta,gamma)
        trans_quat = transq*(self.lig_quat-self.center_rec_quat)*np.conjugate(transq)
        rotate_quat = rotateq*(trans_quat-tmp_axis_quat)*np.conjugate(rotateq)
        return quat.as_float_array(rotate_quat+tmp_axis_quat+self.center_rec_quat)[:,1:]

