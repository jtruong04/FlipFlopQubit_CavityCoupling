import numpy as np
from Scripts.constants import *
from Scripts.pauli import *

def computeDensityMatrix(system, initialState, T, Nt, withNoise=False, dephaseOnly=True, approxInteraction=[1, 1, 1, 1, 1, 1, 1, 1, 1], approxNoise=[1, 1, 1, 1, 1, 1, 1, 1, 1]):
    H = system.getSystemHamiltonian(approxInteraction = approxInteraction)
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
        noise, noiseParameters = system.getNoiseHamiltonian(approxNoise=approxNoise)
        noise = [np.matmul( np.matmul( R_inv, hn ), R ) for hn in noise]
    timeRange = np.arange(0, T, T/Nt)

    W_mat = np.kron(W, np.ones(Nstates))-np.kron(np.ones(Nstates), W)
    rho = np.zeros((Nt, Nstates,Nstates))
    for it, t in enumerate(timeRange):
        eiwt = np.diag(np.exp(W_mat * 1j * t))
        noiseMatrix = np.zeros((Nstates**2, Nstates**2), dtype=complex)
        if withNoise:
            if dephaseOnly:
                Kt = np.zeros((Nstates**2, Nstates**2),dtype=complex)
                for i, hn in enumerate(noise):
                    gamma_squared = np.power((np.kron(np.diag(np.diag(hn)), np.eye(Nstates)) -
                                            np.kron(np.eye(Nstates), np.diag(np.diag(hn))))/hbar, 2)
                    jt = noiseParameters[i]['type'].J(t, 0, 0)
                    Kt = Kt + jt*gamma_squared
                noiseMatrix = np.diag(np.exp(np.diag(-1.0*Kt)))
            else:
                pass
        else:
            noiseMatrix = np.eye(RR.shape[0])
        # print(RR.shape, eiwt.shape, noiseMatrix.shape, RR_INV.shape, initial_rho_vec.shape)
        rho[it] = np.matmul(np.matmul(np.matmul(np.matmul(
                        RR,
                        eiwt),
                        noiseMatrix),
                        RR_INV),
                        initial_rho_vec).reshape(Nstates,Nstates)
    return timeRange, rho

def partialTrace(rho, Na, Nb):
    Nt = rho.shape[0]
    rho_tensor = rho.reshape(Nt, Na, Nb, Na, Nb)
    rho_a = np.trace(rho_tensor, axis1=2, axis2=4)
    rho_b = np.trace(rho_tensor, axis1=1, axis2=3)
    return rho_a, rho_b
