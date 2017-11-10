import numpy as np
import quaternion as quat

def self_rot_quat(axis,alpha,beta,gamma):
    """
    Return quaternion of self-rotation around a given center

    Keyword arguments:
    axis -- the vector from center to origin
    alpha, beta, gamma -- same as Complex.rotation
    """
    axis_2 = np.cross(np.array([0,0,1]),axis)
    axis_3 = np.cross(axis,axis_2)
    rot_1 = (gamma*0.5) * axis/np.linalg.norm(axis)
    rot_2 = (beta*0.5) * axis_2/np.linalg.norm(axis_2)
    rot_3 = (alpha*0.5) * axis_3/np.linalg.norm(axis_3)
    return np.exp(quat.quaternion(*rot_1))*np.exp(quat.quaternion(*rot_2))*np.exp(quat.quaternion(*rot_3))

def rot_around_rec_quat(axis,theta,phi):
    """
    Return quaternion of rotation of a given point around origin.

    Keyword arguments:
    axis -- the vector from origin to the point
    thata, phi -- same definition in spherical coordinations
    """
    phiq = np.exp(quat.quaternion(*(np.array([0.,0.,0.,phi*0.5]))))
    vec = quat.quaternion(*axis)
    v_prime = phiq*vec*np.conjugate(phiq)
    tmp_center_lig = (quat.as_float_array(v_prime))[1:]
    t_axis = np.cross(np.array([0,0,1]),tmp_center_lig)
    t_rot_axis = (theta*0.5) * t_axis/np.linalg.norm(t_axis)
    thetaq = np.exp(quat.quaternion(*t_rot_axis))
    return thetaq*phiq

