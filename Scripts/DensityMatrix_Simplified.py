import numpy as np
from Scripts.constants import *
from Scripts.pauli import *


def computeDensityMatrix(system, initialState, T, Nt):
    H = system.getSystemHamiltonian()
    # print(H)
    energies, states = np.linalg.eig(H)
    # print(energies)
    # print(states)
    states_inverse = np.linalg.inv(states)
    num_states = H.shape[0]
    # DONE Create time series
    timeRange = np.arange(0, T, T/Nt)
    # DONE Calculate Oscillation Matrix: Omega = Exp(-iwt)
    wjk = ((np.kron(energies, np.ones(num_states))
            - np.kron(np.ones(num_states), energies))/hbar).reshape(num_states, num_states)
    Omega = np.exp(-1j * wjk[np.newaxis, :] *
                   timeRange[:, np.newaxis, np.newaxis])
    # print(Omega)
    # DONE Calculate Decay Matrix: G = Exp(- sum J_i(t) Gamma_ijk^2 )
    noise, noiseParameters = system.getNoiseHamiltonian()
    # print(" Prerotated Noise: ", noise)
    noise = [np.diag(np.matmul(np.matmul(states_inverse, hn), states))/hbar
             for hn in noise]
    # print(" Rotated noise: ", noise)
    gamma_squared = [((np.kron(hn, np.ones(num_states)) - np.kron(np.ones(num_states), hn))
                      ** 2).reshape(num_states, num_states) for hn in noise]
    print(gamma_squared)
    Ki = [-np.vectorize(param['type'].J)(timeRange)[:, np.newaxis, np.newaxis] *
          gsquared[np.newaxis, :] for gsquared, param in zip(gamma_squared, noiseParameters)]
    G = np.exp(sum(Ki))
    # DONE Calculate Initial density Matrix in eigenbasis (rhoprime0 = R^-1 rho R)
    initial_rho = np.kron(initialState.transpose(), initialState)
    rhoprime0 = np.matmul(np.matmul(states_inverse, initial_rho), states)
    # DONE Calculate Density Matrix in eigenbasis (rhoprime = Omega * G * rhoprime0)
    rhoprime_woNoise = Omega * rhoprime0
    rhoprime_wNoise = rhoprime_woNoise * G
    # DONE Calculate Density Matrix in original basis (rho = R rhoprime R^-1)
    rho_woNoise = np.matmul(
        np.matmul(states, rhoprime_woNoise), states_inverse)
    rho_wNoise = np.matmul(np.matmul(states, rhoprime_wNoise), states_inverse)
    # Return timeRange, rho, rho_woNoise
    return timeRange, rho_wNoise, rho_woNoise


def partialTrace(rho, Na, Nb):
    Nt = rho.shape[0]
    rho_tensor = rho.reshape(Nt, Na, Nb, Na, Nb)
    rho_a = np.trace(rho_tensor, axis1=2, axis2=4)
    rho_b = np.trace(rho_tensor, axis1=1, axis2=3)
    return rho_a, rho_b
