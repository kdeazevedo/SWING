import numpy as np
import re

ROT = r'[()]'
LINE = r'(?P<date>\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3}) - +(?P<name>\w+) - +(?P<level>\w+) - (?P<detail>.+)\n'

class Log():
    
    def __init__(self):
        self.failed = []

    @classmethod
    def from_file(cls,path):
        log = cls()
        angles = []
        for line in open(path,'r'):
            dic = re.match(LINE,line).groupdict()
            log_name = dic['name']
            if log_name == 'main':
                tmp = list(re.findall('-?\d+\.\d+',dic['detail']))
                if not len(tmp) == 0:
                    angles.append(tmp)


            elif log_name == 'minimizer':
                ind = int(re.findall('No\. (\d+)',dic['detail'])[0])
                log.failed.append(ind)
        log.rot_angles = np.array(angles)
        return log
        

            
if __name__ == '__main__':
    log = Log.from_file('out/log_20171201_104411.txt')
    print(log.failed)
    print(log.rot_angles)
