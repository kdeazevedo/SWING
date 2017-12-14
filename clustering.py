import numpy as np
from scipy.cluster.hierarchy import linkage, cut_tree
from scipy.spatial.distance import pdist
import dockerasmus.pdb as pdb
import models.protein

def cluster(lst,threshold,distance,method='single'):
    """
    Return proteins clustering from a set of protein.
    The function uses hierarchical clustering(linkage) and cut_tree from SciPy to create protein clu
    sters.

    Keyword arguments:
    lst -- An iterable object, in which each element is an object of dockerasmus.pdb.Protein
    threshold -- height of cut_tree from where the hierarchical tree is cut
    distance -- 'rmsd' or 'center' for distance between geometrical center
    method -- method of linkage. 'single' and 'complete' are used for aligned ligands and 
              post-scoring clustering 
    """
    if distance == 'center':
        data = [np.mean(p.atom_positions(),axis=0) for p in lst]
        data_dist = pdist(data)
    elif distance == 'rmsd':
        data = [p.atom_positions() for p in lst]
        l1 = []
        l2 = []
        l = len(data)
        for i in range(l):
            l1 += [data[i]]*(l-i-1)
            l2 += data[i+1:]
        data_dist = np.sqrt(np.mean(np.apply_over_axes(np.sum,l**2,2),axis=1)).ravel()

    # Start hierarchy
    data_link = linkage(data_dist,method=method)
    return [t[0] for t in cut_tree(data_link,height=threshold)]
        
def cluster_from_files(lst,threshold,distance,method='single'):
    """
    Same as cluster except that the function takes an iterable object of pdb files' path as input 
    """
    return cluster(map(pdb.Protein.from_pdb_file,lst),threshold,distance,method=method)

def cut_to_cluster(lst):
    """
    Takes result from cluster.
    Retruns a dictionary in which key is the label of cluster and its value is a list of element in 
    the same cluster
    """
    d = {}
    for idx, v in enumerate(lst):
        try:
            d[v].append(idx)
        except:
            d[v] = [idx]
    return d

def interologs_cluster(dico):
    """
    Function to cluster aligned ligands using distance of geometric center
    Takes iterologs dictionary (same as sampling configuration file) and creates a new one from it 
    by keeping the first aligned ligand in each cluster
    The height to cut is set to 3A
    """
    dico_t = {}
    keys = list(dico.keys())
    files = [v['lig_aligned'] for v in dico.values()]
    clr = cluster_from_files(files,3,'center')
    for v in cut_to_cluster(clr).values():
        k = keys[v[0]]
        dico_t[k] = dico[k]
    return dico_t


