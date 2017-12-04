import numpy as np
import re

LINE = r'(?P<date>\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3}) - +(?P<name>\w+) - +(?P<level>\w+) - (?P<detail>.+)\n'
FSTLINE = r'Start program\. Receptor: (.*) Ligand: (.*) Output: (.*)'


class Log():
    
    def __init__(self):
        self.lig = ""
        self.rec = ""
        self.failed = {}
        self.rot_angles = {}
        self.templates = []

    @classmethod
    def from_file(cls,path):
        log = cls()
        angles = {}
        failed = {}
        key = ""
        for line in open(path,'r'):
            dic = re.match(LINE,line).groupdict()
            if dic['name'] == 'main':
                # Get receptor and ligand
                if dic['detail'].startswith("Start program"):
                    l = re.match(FSTLINE,dic['detail']).groups()
                    log.rec = l[0]
                    log.lig = l[1]
                    log.out = l[2]
            elif dic['name'] == 'sampling':
                # Get template
                if dic['detail'].startswith('Start'):
                    key = re.findall(r'template (\w+)',dic['detail'])[0]
                    angles[key] = []
                    failed[key] = []
                    log.templates.append(key)
                # Get rotation angles
                else:
                    tmp = list(re.findall('-?\d+\.\d+',dic['detail']))
                    angles[key].append(tmp)
            elif dic['name'] == 'minimizer':
                ind = int(re.findall('No\. (\d+)',dic['detail'])[0])
                failed[key].append(ind)
        log.rot_angles = {k:np.array(v) for (k,v) in angles.items()}
        log.failed = failed
        return log
        
    def get_angles(self,key=None,failed=False):
        if not key == None:
            if failed:
                return self.rot_angles[key]
            else:
                return np.delete(self.rot_angles[key],self.failed[key],axis=0)
        else:
            if failed:
                return self.rot_angles
            else:
                return {k: np.delete(self.rot_angles[k],self.failed[k],axis=0) for k in self.templates}
    def get_sampling_ind(self,key=None):
        try:
            self._sampling_ind
        except:
            self._sampling_ind = {k:np.delete([i for i in range(v.shape[0])],self.failed[k]) for k,v in self.rot_angles.items()}
        if not key == None:
            return self._sampling_ind[key]
        else:
            return self._sampling_ind
            
if __name__ == '__main__':
    log = Log.from_file('out_1/log_20171203_140909.txt')
    print(log.rec)
    print(log.templates)
    print(log.get_angles())
    print(log.get_sampling_ind())
