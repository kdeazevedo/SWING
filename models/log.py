import numpy as np
import re

ROT = r'[()]'
LINE = r'(?P<date>\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3}) - +(?P<name>\w+) - +(?P<level>\w+) - (?P<detail>.+)\n'

class Log():
    
    def __init__(self):
        self.lig = ""
        self.rec = ""
        self.failed = {}
        self.rot_angles = {}

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
                    l = re.findall(r': (\w+)',dic['detail'])
                    log.rec = l[0]
                    log.lig = l[1]
            elif dic['name'] == 'sampling':
                # Get template
                if dic['detail'].startswith('Start'):
                    key = re.findall(r'template (\w+)',dic['detail'])[0]
                    angles[key] = []
                    failed[key] = []
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
        

            
if __name__ == '__main__':
    log = Log.from_file('out/log_20171202_160626.txt')
    print(log.failed)
    print(log.rot_angles)
