'''
A function that performs block-boostrap resampling 
This script is based on Julien's block_bootstrap.m
'''

import numpy as np

def block_bootstrap_jeg(X, Lb, Nb):
    np.random.seed(3)
    nt = len(X)
    ns = int(np.ceil(nt/Lb))
    Xb = np.zeros((Nb, nt))

# naive sampling as in block_boostrap.m
    for b in range(Nb):
        sample = np.random.choice(ns, ns, replace = True)
        for j in range(ns):
            Xb[b,j*Lb:(j+1)*Lb] = X[sample[j]*Lb:(sample[j]+1)*Lb]
            
    return Xb
    
def block_bootstrap_yuxin(X, Lb, Nb):
    np.random.seed(3)
    nt = len(X)
    ns = int(np.ceil(nt/Lb))
    Xb = np.zeros((Nb, nt))
    Xblock = X.reshape(ns, Lb)

    # sample with replacement
    for b in range(Nb):
        # Xb[b,:] = np.random.permutation(Xblock).flatten()
        Xb[b,:] = Xblock[np.random.randint(ns,size=ns)].flatten()
        # the high boundary is added by 1 because the high end is open interval
        
    return Xb

def block_bootstrap_ET(X, Lb, Nb): 
    np.random.seed(3)
    nt = len(X)
    ns = int(np.ceil(nt/Lb))
    Xb = np.zeros((Nb, nt))

#  Implement Block Bootstrap as in: 
# http://nbviewer.ipython.org/github/welch/stats-notebooks/blob/master/SubsamplingBootstrap.ipynb
    for b in range(Nb):       
        for block_i, start in enumerate(np.random.randint(nt - Lb + 1, size=ns)):
            Xb[b,block_i*Lb:(block_i+1)*Lb] = X[start:start+Lb]
        
            
    return Xb
