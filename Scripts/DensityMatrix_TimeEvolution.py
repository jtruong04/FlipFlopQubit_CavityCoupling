import numpy as np
from Scripts.constants import *
from Scripts.pauli import *

def computeDensityMatrix(system, initialState, T, Nt, withNoise = False, dephaseOnly = True):
    H = system.getSystemHamiltonian()
    # Get eigenvalues and eigenvectors
    W, R = np.linalg.eig(H)
    Nstates = system.Nstates
    # Initial density matrix
    initial_rho = np.kron(initialState.transpose(), initialState)
    initial_rho_vec = initial_rho.flatten()
    # R matrices in Liouville-Fock space
    RR = np.kron(R, R)
    R_inv = np.linalg.inv(R)
    RR_INV = np.kron(R_inv, R_inv)
    # Get noise in eigenbasis.
    if withNoise:
        noise, noiseParameters = system.getNoiseHamiltonian()
        noise = [np.matmul( np.matmul( R_inv, hn ), R ) for hn in noise]
    # Oscillating term in L-F space. Is time dependent. We'll shape this as (N_t, ROW, COLS)
    timeRange = np.arange(0, T, T/Nt)
    W_mat = np.kron(W, np.ones(Nstates))-np.kron(np.ones(Nstates), W)
    eiwtList = np.exp(np.kron(timeRange[:, np.newaxis], W_mat)*1j)
    eiwt = np.zeros((Nt, Nstates**2, Nstates**2))
    for t in range(Nt):
        eiwt[t].flat[slice(0, None, 1+Nstates**2)] = eiwtList[t]

    # Build noise decay matrix
    noiseMatrix = np.zeros((Nt, Nstates**2, Nstates**2))

    if withNoise:
        # noiseMatrix = np.eye(RR.shape[0])
        if dephaseOnly:
            Kt = 0
            for i, hn in enumerate(noise):
                gamma_squared = np.power((np.kron(np.diag(np.diag(hn)), np.eye(Nstates)) - \
                    np.kron(np.eye(Nstates),np.diag(np.diag(hn))))/hbar,2)
                jt = np.asarray([noiseParameters[i]['type'].J(t, 0, 0)
                                 for t in timeRange])
                Kt += jt[:, np.newaxis, np.newaxis]*gamma_squared
            for i, K in enumerate(Kt):
                # print(K)
                noiseMatrix[i] = np.diag(np.exp(np.diag(-1.0*K)))
    #         for it, t in enumerate(timeRange):
    #             Klist = np.zeros((1, Nstates**2))
    #             for i, hn in enumerate(noise):
    #                 Klist = Klist + (np.kron(np.power(np.diag(hn), 2), np.ones(Nstates)) -
    #                                  np.kron(np.ones(Nstates), np.power(np.diag(hn), 2))) * noiseParameters[i]['type'].J(t, 0, 0)
    #             noiseMatrix[it].flat[slice(0, None, 1+Nstates**2)]=np.exp(Klist)
        else:
            pass
    else:
        noiseMatrix = np.eye(RR.shape[0])

    # print(noiseMatrix)

    rho = (np.matmul(np.matmul(np.matmul(np.matmul(RR, eiwt), noiseMatrix), RR_INV),
                     initial_rho_vec)).reshape(Nt, Nstates, Nstates)
    return timeRange, rho

def partialTrace(rho, Na, Nb):
    Nt = rho.shape[0]
    rho_tensor = rho.reshape(Nt, Na, Nb, Na, Nb)
    rho_a = np.trace(rho_tensor, axis1=2, axis2=4)
    rho_b = np.trace(rho_tensor, axis1=1, axis2=3)
    return rho_a, rho_b
