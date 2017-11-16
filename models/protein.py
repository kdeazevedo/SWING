import dockerasmus.pdb as pdb
import numpy as np
import models.residue

LETTERS = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G',
    'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T',
    'TRP':'W','TYR':'Y','VAL':'V'}


def iter_chain_atoms(self):
    """
    Yield every atom and its chain's id in ``self``. Similar to self.iteratoms()
    """
    for chain in self.itervalues():
        for residue in chain.itervalues():
            for atom in sorted(residue.itervalues(), key=lambda a: a.id):
                yield atom, chain.id

pdb.Protein.iter_chain_atoms = iter_chain_atoms

def write_atoms(self,path,pos=None):  
    """
    A naive protein's atoms' records printer.

    Keyword arguments:
    path -- the output file's path
    """
    atom_template = '{:7}{:4d}  {:4}{:3} {:1} {:3d}      {:5.3f}  {:5.3f}  {:5.3f}'
    out = open(path,'w')
    for idx, (atom, c) in enumerate(self.iter_chain_atoms()):
        if pos is not None:
            l = pos[idx]
        else:
            l = [atom.x,atom.y,atom.z]
        print(atom_template.format("ATOM",atom.id,atom.name,atom.residue.name,c,atom.residue.id,l[0],l[1],l[2]), file=out)
    print("END", file=out)
    out.close()

pdb.Protein.write_atoms = write_atoms

def get_ca_ind(self):
    """
    Return list(numpy.array) of carbon alpha index in self.atom_positions
    """
    try:
        return self._ca_ind
    except:
        self._ca_ind = np.where(np.array([a.name.lower() for a in self.iteratoms()])=='ca')[0]
        return self._ca_ind

pdb.Protein.get_ca_ind = get_ca_ind

def get_ca(self):
    try:
        return self._ca
    except:
        self._ca = self.atom_positions()[self.get_ca_ind()]
        return self._ca

pdb.Protein.get_ca = get_ca

def write_fasta(self,path):
    """
    Write protein's residus into a given file(format FASTA) 
    """
    f = open(path,"w")
    for chain in self.itervalues():
        print(">"+chain.id,file=f)
        print("".join([LETTERS[r.name] for r in chain.itervalues()]),file=f)
    f.close()

pdb.Protein.write_fasta = write_fasta
