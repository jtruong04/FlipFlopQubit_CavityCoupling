import numpy as np

s0 = np.array([[1,0],[0,1]])
s1 = np.array([[0,1],[1,0]])
s2 = np.array([[0,-1j],[1j,1]])
s3 = np.array([[1,0],[0,-1]])
sx = s1
sy = s2
sz = s3

pauli = (s0,s1,s2,s3)

def a(n):
    '''
        returns the annihilation operator for maximum n photons
    '''
    a = np.zeros((n+1,n+1))
    b = np.arange(1,n+1)
    np.fill_diagonal(a[:,1:],np.sqrt(b))
    return a

def adagger(n):
    return a(n).conj().T

def kron(ms):
    if len(ms) == 1:
        return ms[0]
    elif len(ms) == 2:
        if(ms[0].shape[0]==0):
            return ms[1]
        elif(ms[1].shape[0]==0):
            return ms[0]
        else:
            return np.kron(ms[0],ms[1])
    else:
        if(ms[0].shape[0]==0):
            return kron(ms[1:])
        else:
            return np.kron(ms[0],kron(ms[1:]))