import dockerasmus.pdb as pdb

def write_atoms(self,path):  
    atom_template = '{:7}{:4d}  {:4}{:3} {:1} {:3d}      {:5.3f}  {:5.3f}  {:5.3f}'
    out = open(path,'w')
    for atom in self.iteratoms():
        print(atom_template.format("ATOM",atom.id,atom.name,"ARG","A",123,atom.x,atom.y,atom.z), file=out)
    print("END", file=out)
    out.close()

pdb.Protein.write_atoms = write_atoms

