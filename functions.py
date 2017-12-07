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


def angles_generator(k,deg=0,rot_lim=np.pi/24,self_lim=np.pi/24,dist="uniform"):
    """
    Generate randomly k*5 rotation angles(2 for rotation around center and 3 for self rotation)
    thata and phi are selected from [-rot_lim,rot_lim]*(1-deg)
    alpha, beta, and gamma are selected from [-self_lim,self_lim]*(1-deg)
    Each element is a dim 1 numpy array [theta,phi,alpha,beta,gamma]

    Keyword arguments:
    k -- number of 5 rotation angles to generate
    deg -- value between 0 and 1. Default 0
    rot_lin -- maximum value of rotation angles around center in radian. Default pi/6
    self_lim -- maximum value of self-rotation angles. Default pi/6
    dist -- sampling distribution. Value "uniform" or "normal". Default "uniform"
            If dist == "normal", the probability of a rotation angle in [-lim, lim] is 0.997
    """
    for i in range(k):
        m = rot_lim * (1-deg)
        n = self_lim * (1-deg)
        if dist == "uniform":
            yield np.concatenate((np.random.sample(2)*m*2-m,np.random.sample(3)*2*n-n),axis=0)
        elif dist == "normal":
            yield np.concatenate((np.random.normal(scale=m/3,size=2),np.random.normal(scale=n/3,size=3)),axis=0)

def angles_random(deg=0,rot_lim=np.pi/24,self_lim=np.pi/24,dist='uniform'):
    """
    Generate a random set of rotation angles. Same as :
    angles_generator(1,deg,rot_lim,self_lim,dist)
    """
    m = rot_lim * (1-deg)
    n = self_lim * (1-deg)
    if dist == "uniform":
        return np.concatenate((np.random.sample(2)*m*2-m,np.random.sample(3)*2*n-n),axis=0)
    elif dist == "normal":
        return np.concatenate((np.random.normal(scale=m/3,size=2),np.random.normal(scale=n/3,size=3)),axis=0)

def vec_to_dist(a,b,d):
    """
    Return vector b->c such that ||ac|| = d and a->b // a->c
    """
    m = np.linalg.norm(b-a)
    return ((d-m)*b+(m-d)*a)/m

def rmsd(p1, p2):
    return p1.rmsd(p2)

if __name__ == '__main__':
    for l in angles_generator(10):
        print(l)
