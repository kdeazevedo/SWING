import numpy as np
import pickle

def get(seed):
    """
    Return numpy randomstate from a given seed file if file exists
    Otherwise, create a new random state and store into seed.
    """
    try:
        with open(seed,'rb') as f:
            return pickle.load(f)
    except:
        st = np.random.get_state()
        with open(seed,'wb') as f:
            pickle.dump(st,f)
        return st

if __name__ == '__main__':
    print(get('test.seed'))
