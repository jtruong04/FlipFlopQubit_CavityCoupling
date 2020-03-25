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