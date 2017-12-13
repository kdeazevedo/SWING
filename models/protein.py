#Modified dockerasmus.pdb.protein.py
#Original Author: Martin Larralde (althonos)
#Modified by Hua-Ting Yao ([anthony]htyao)
import dockerasmus.pdb as pdb
import numpy as np
from dockerasmus.pdb import Chain, Atom, Residue
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
    atom_template = '{:6}{:5d}  {:3} {:3} {:1}{:4d}    {:8.3f}{:8.3f}{:8.3f}'
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
    """
    Return list of carbon alpha positions
    """
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


@classmethod
def from_pdb(cls, handle):
    """
    Create a new Protein object from a PDB file handle.
    Arguments:
    handle (file handle): a file-like object opened in
    binary read mode (must be line-by-line iterable).
    """
    protein = cls()
    for line in handle:
        if line.startswith(b"ATOM  "):
            atom = cls._parse_pdb_atom_line(line)
        elif line.startswith(b"HETATM") and line.split()[3] == b'MSE':
            atom = cls._parse_pdb_atom_line(line.replace(b"HETATM",b"ATOM  ").replace(b"MSE",b"MET"))
        else:
            continue
        if atom['chainID'] not in protein:
            protein[atom['chainID']] = Chain(atom['chainID'])
        if atom['resSeq'] not in protein[atom['chainID']]:
            protein[atom['chainID']][atom['resSeq']] = Residue(atom['resSeq'], atom['resName'])
        protein[atom['chainID']][atom['resSeq']][atom['name']] = Atom(
        atom['x'], atom['y'], atom['z'], atom['serial'], atom['name'],
             protein[atom['chainID']][atom['resSeq']],
        )
    return protein
pdb.Protein.from_pdb = from_pdb
