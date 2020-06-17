import numpy as np
from Scripts.constants import *
from Scripts.pauli import *

# def computeDensityMatrix(system, initialState, T, Nt, dephaseOnly=True, approxInteraction=[1, 1, 1, 1, 1, 1, 1, 1, 1], approxNoise=[1, 1, 1, 1, 1, 1, 1, 1, 1]):
#     H = system.getSystemHamiltonian(approxInteraction = approxInteraction)
#     # Get eigenvalues and eigenvectors
#     W, R = np.linalg.eig(H)
#     print(" R = ", R)
#     Nstates = system.Nstates
#     # Initial density matrix
#     initial_rho = np.kron(initialState.transpose(), initialState)
#     initial_rho_vec = initial_rho.flatten()
#     # R matrices in Liouville-Fock space
#     RR = np.kron(R, R)
#     R_inv = np.linalg.inv(R)
#     print('Rinv = ', R_inv)
#     RR_INV = np.kron(R_inv, R_inv)
#     # Get noise in eigenbasis.
#     noise, noiseParameters = system.getNoiseHamiltonian(approxNoise=approxNoise)
#     print("hL: ", noise[1])
#     noise = [np.matmul( np.matmul( R_inv, hn ), R ) for hn in noise]

#     timeRange = np.arange(0, T, T/Nt)

#     gamma_squared = np.zeros((system.Nn,Nstates**2,Nstates**2),dtype = complex)
#     for i, hn in enumerate(noise):
#         # print(hn)
#         gamma_squared[i] = np.power((np.kron(np.diag(np.diag(hn)), np.eye(Nstates)) -
#                                             np.kron(np.eye(Nstates), np.diag(np.diag(hn))))/hbar, 2)

#     print("hL: ", noise[1])

#     print([np.diag(gamma2*hbar**2).reshape(Nstates,Nstates) for gamma2 in gamma_squared])
#     W_mat = np.kron(W, np.ones(Nstates))-np.kron(np.ones(Nstates), W)
#     print(W_mat)
#     rho = np.zeros((Nt, Nstates,Nstates))
#     rho_zeroNoise = np.zeros((Nt, Nstates,Nstates))
#     rrinvrho = np.matmul(RR_INV,initial_rho_vec)
#     for it, t in enumerate(timeRange):
#         # if ( it % 25 == 0 ):
#         #     print(t)
#         eiwt = np.diag(np.exp(W_mat * 1j * t / hbar))
#         noiseMatrix = np.zeros((Nstates**2, Nstates**2), dtype=complex)
#         if dephaseOnly:
#             Kt = np.zeros((Nstates**2, Nstates**2),dtype=complex)
#             for i, hn in enumerate(noise):
#                 # gamma_squared = np.power((np.kron(np.diag(np.diag(hn)), np.eye(Nstates)) -
#                 #                         np.kron(np.eye(Nstates), np.diag(np.diag(hn))))/hbar, 2)
#                 jt = noiseParameters[i]['type'].J(t, 0, 0)
#                 Kt = Kt + jt*gamma_squared[i]
#             noiseMatrix = np.diag(np.exp(np.diag(-1.0*Kt)))
#         else:
#             pass
#         # print(RR.shape, eiwt.shape, noiseMatrix.shape, RR_INV.shape, initial_rho_vec.shape)
#         rho[it] = np.matmul(np.matmul(np.matmul(
#                         RR,
#                         eiwt),
#                         noiseMatrix),
#                         rrinvrho).reshape(Nstates,Nstates)
#         rho_zeroNoise[it] = np.matmul(np.matmul(
#                 RR,
#                 eiwt),
#                 rrinvrho).reshape(Nstates,Nstates)
#     return timeRange, rho, rho_zeroNoise

def computeDensityMatrix(system, initialState, T, Nt):
    H = system.getSystemHamiltonian()
    W, R = np.linalg.eig(H)
    Nstates = system.Nstates
    # Initial density matrix
    initial_rho = np.kron(initialState.transpose(), initialState)
    # initial_rho_vec = initial_rho.flatten()
    # R matrices in Liouville-Fock space
    # RR = np.kron(R, R)
    R_inv = np.linalg.inv(R)
    # print('Rinv = ', R_inv)
    # RR_INV = np.kron(R_inv, R_inv)
    # Get noise in eigenbasis.
    noise, noiseParameters = system.getNoiseHamiltonian()
    # print("hL: ", noise[1])
    noise = [np.diag(np.matmul( np.matmul( R_inv, hn ), R )) for hn in noise]

    timeRange = np.arange(0, T, T/Nt)

    gamma_squared = np.zeros((system.Nn,Nstates**2,Nstates**2),dtype = complex)
    for i, hn in enumerate(noise):
        # print(hn)
        gamma_squared[i] = (np.kron(noise[i], np.ones(Nstates)) - np.kron(np.ones(Nstates), noise[i])**2).reshape(Nstates,Nstates)
        # gamma_squared[i] = np.power((np.kron(np.diag(np.diag(hn)), np.eye(Nstates)) -
        #                                     np.kron(np.eye(Nstates), np.diag(np.diag(hn))))/hbar, 2)

    # print("hL: ", noise[1])
    coefs = np.zeros((Nstates,Nstates,Nstates,Nstates))
    for a in range(Nstates):
        for b in range(Nstates):
            for j in range(Nstates):
                for k in range(Nstates):
                    coefs[a,b,j,k] = R[a,j] * R[b,k] * np.matmul(np.matmul(R_inv, initial_rho),np.transpose(R_inv))[j,k]
    # print([np.diag(gamma2*hbar**2).reshape(Nstates,Nstates) for gamma2 in gamma_squared])
    W_mat = np.kron(W, np.ones(Nstates))-np.kron(np.ones(Nstates), W)/hbar
    # print(W_mat)
    rho = np.zeros((Nt, Nstates,Nstates))
    rho_zeroNoise = np.zeros((Nt, Nstates,Nstates))
    # rrinvrho = np.matmul(RR_INV,initial_rho_vec)
    for it, t in enumerate(timeRange):
        # if ( it % 25 == 0 ):
        #     print(t)
        # eiwt = np.diag(np.exp(W_mat * 1j * t / hbar))
        # noiseMatrix = np.zeros((Nstates**2, Nstates**2), dtype=complex)
        # if dephaseOnly:
        #     Kt = np.zeros((Nstates**2, Nstates**2),dtype=complex)
        #     for i, hn in enumerate(noise):
        #         # gamma_squared = np.power((np.kron(np.diag(np.diag(hn)), np.eye(Nstates)) -
        #         #                         np.kron(np.eye(Nstates), np.diag(np.diag(hn))))/hbar, 2)
        #         jt = noiseParameters[i]['type'].J(t, 0, 0)
        #         Kt = Kt + jt*gamma_squared[i]
        #     noiseMatrix = np.diag(np.exp(np.diag(-1.0*Kt)))
        # else:
        #     pass
        # print(RR.shape, eiwt.shape, noiseMatrix.shape, RR_INV.shape, initial_rho_vec.shape)
        Kt = np.zeros((Nstates,Nstates))
        for i, hn in enumerate(noise):
            # gamma_squared = np.power((np.kron(np.diag(np.diag(hn)), np.eye(Nstates)) -
            #                         np.kron(np.eye(Nstates), np.diag(np.diag(hn))))/hbar, 2)
            jt = noiseParameters[i]['type'].J(t, 0, 0)
            Kt = Kt + jt*gamma_squared[i]/(hbar**2)        
        rho[it] = (coefs*np.cos(W_mat*t)*np.exp(-Kt)).reshape(Nstates,Nstates)
        rho_zeroNoise[it] = np.matmul(np.matmul(
                RR,
                eiwt),
                rrinvrho).reshape(Nstates,Nstates)
    return timeRange, rho, rho_zeroNoise

def partialTrace(rho, Na, Nb):
    Nt = rho.shape[0]
    rho_tensor = rho.reshape(Nt, Na, Nb, Na, Nb)
    rho_a = np.trace(rho_tensor, axis1=2, axis2=4)
    rho_b = np.trace(rho_tensor, axis1=1, axis2=3)
    return rho_a, rho_b
