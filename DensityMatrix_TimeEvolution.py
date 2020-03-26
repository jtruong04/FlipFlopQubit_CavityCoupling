import numpy as np
from constants import *
from pauli import *

class FlipFlopSystem:
    def __init__(self, Nm, Np, parameters_cavity, Nd, parameters_qubits, Nn = 0, parameters_noise = []):
        '''
        Units being used for analysis
            Energy: GHz
            Time  : ns
        Arguments:
            Nd : Number of donors
            parameters_qubits : list of size Nd dicts containing tunable parameters
                {'Vt','wB','eps'}
            Nm : Number of cavity modes
            Np : Max number of photons per mode
            parameters_cavity : list of size Nm dicts containing cavity parameters
                {'wc', 'gc'}
            Nn : Number of noise sources
            parameters_noise : list of size Nn
                {'type', 'wn', 'effected_donors'}
                type: 'pink', 'white', etc.
                wn: strength of noise
                effected_donors: boolean list of len Nd
        '''
        self.Nq = 4
        self.Nd = Nd
        self.Np = Np
        self.Nm = Nm
        self.Nn = Nn
        self.Nstates = (self.Np + 1)**Nm * (self.Nq**self.Nd)
        assert(len(parameters_qubits) == self.Nd)
        self.parameters_qubits = parameters_qubits
        assert(len(parameters_cavity) == self.Nm)
        self.parameters_cavity = parameters_cavity
        assert(len(parameters_noise) == self.Nn)
        self.parameters_noise = parameters_noise
        # print(parameters_cavity)
        # print(parameters_qubits)
        self.Z1q = []
        self.E1q = []

    def FlipFlopEnergies(self):
        '''
            Returns the eigenenergies for the flip-flip qubit given above
            arguments:
                parameters{
                    'Vt': tunnel coupling,
                    'wB': zeeman splitting,
                    'eps':applied electric field energy
                }
            returns:
                list of 4x4 matrices with energies along the diagonal
        '''
        self.E1q = []
        for i in range(self.Nd):
            Vt = self.parameters_qubits[i]['Vt']
            wB = self.parameters_qubits[i]['wB']
            eps = self.parameters_qubits[i]['eps']
            w0 = np.sqrt(eps**2 + Vt**2)
            eta = np.arctan2(Vt, eps)
            E_ff = np.zeros((4, 4))
            E_ff[0, 0] = -0.5*(w0+wB) - 0.125*hyperfine*(1-np.cos(eta))-0.25*Delta*wB*(1+np.cos(eta))-0.0625*(hyperfine**2/wB)*(
                (1-np.cos(eta))**2 + np.sin(eta)**2*(wB/(4*w0) + wB/(w0+wB) - Delta*wB**2/(hyperfine*w0) + Delta**2*wB**3/(w0*hyperfine**2)))
            E_ff[1, 1] = -0.5*(w0-wB) - 0.125*hyperfine*(1-np.cos(eta))+0.25*Delta*wB*(1+np.cos(eta))-0.0625*(hyperfine**2/wB)*(-(
                1-np.cos(eta))**2 + np.sin(eta)**2*(wB/(4*w0) + wB/(w0+wB) + Delta*wB**2/(hyperfine*w0) + Delta**2*wB**3/(w0*hyperfine**2)))
            E_ff[2, 2] = 0.5*(w0-wB) - 0.125*hyperfine*(1+np.cos(eta))-0.25*Delta*wB*(1-np.cos(eta))-0.0625*(hyperfine**2/wB)*(
                (1+np.cos(eta))**2 + np.sin(eta)**2*(-wB/(4*w0) - wB/(w0+wB) + Delta*wB**2/(hyperfine*w0) - Delta**2*wB**3/(w0*hyperfine**2)))
            E_ff[3, 3] = 0.5*(w0+wB) - 0.125*hyperfine*(1+np.cos(eta))+0.25*Delta*wB*(1-np.cos(eta))-0.0625*(hyperfine**2/wB) * \
                (-(1+np.cos(eta))**2 + np.sin(eta)**2*(-wB/(4*w0) - wB/(w0+wB) -
                                                    Delta*wB**2/(hyperfine*w0) - Delta**2*wB**3/(w0*hyperfine**2)))
            self.E1q.append(E_ff)
        return self.E1q

    def ElectronPosition(self):
        '''
            Returns the coefficients z_jk for the flip-flip qubit given above
            arguments:
                parameters{
                    'Vt': tunnel coupling,
                    'wB': zeeman splitting,
                    'eps':applied electric field energy
                }
            returns:
                list of 4x4 matrix
        '''
        self.Z1q = []
        for i in range(self.Nd):
            Vt = self.parameters_qubits[i]['Vt']
            wB = self.parameters_qubits[i]['wB']
            eps = self.parameters_qubits[i]['eps']
            w0 = np.sqrt(eps**2 + Vt**2)
            eta = np.arctan2(Vt, eps)
            Zcoef = np.zeros((4, 4))
            Zcoef[0, 1] = -(hyperfine*w0*Delta*wB*np.cos(eta) *
                            np.sin(eta)**2)/(4*wB*(w0**2-wB**2))
            Zcoef[0, 3] = (hyperfine**2*w0**3*np.cos(eta) *
                        np.sin(eta)**2)/(4*wB*(w0**2-wB**2)**2)
            Zcoef[1, 0] = np.sin(eta) + (hyperfine*np.cos(eta)*np.sin(eta))/(4*w0)
            Zcoef[1, 1] = -(hyperfine*w0*np.cos(eta)*np.sin(eta))/(2*(w0**2-wB**2))
            Zcoef[1, 3] = -(Delta*wB*np.cos(eta)*np.sin(eta))/(2*w0)
            Zcoef[2, 2] = -(hyperfine*w0**2*np.cos(eta) *
                            np.sin(eta))/(2*wB*(w0**2-wB**2))
            Zcoef[3, 0] = np.cos(eta) - (hyperfine*np.sin(eta)**2)/(4*w0)
            Zcoef[3, 1] = (hyperfine*w0*np.sin(eta)**2)/(2*(w0**2-wB**2))
            Zcoef[3, 3] = (Delta*wB*np.sin(eta)**2)/(2*w0)
            z = np.zeros((self.Nq, self.Nq))
            for j in np.arange(self.Nq):
                for k in np.arange(self.Nq):
                    z = z + np.kron(pauli[j], pauli[k])*Zcoef[j, k]
            self.Z1q.append(z)
        return self.Z1q

    def getQubitHamiltonian(self):
        H_qubit = np.zeros((np.power(self.Np+1,self.Nm) * np.power(self.Nq, self.Nd),
                            np.power(self.Np+1, self.Nm) * np.power(self.Nq, self.Nd)))
        ff_energies = self.FlipFlopEnergies()
        for donor in range(self.Nd):
            H_qubit = H_qubit + kron([np.eye(np.power(self.Np+1,self.Nm)), np.eye(
                np.power(self.Nq, donor)), ff_energies[donor], np.eye(np.power(self.Nq, (self.Nd-1-donor)))])
        return H_qubit

    def getCavityHamiltonian(self):
        H_cavity = np.zeros((np.power(self.Np+1, self.Nm) * np.power(self.Nq, self.Nd),
                            np.power(self.Np+1, self.Nm) * np.power(self.Nq, self.Nd)))
        for mode in range(self.Nm):
            # print(self.parameters_cavity[mode])
            wc = self.parameters_cavity[mode]['wc']
            H_cavity = H_cavity + kron([np.eye(mode),np.matmul(adagger(self.Np),a(self.Np)),np.eye(self.Nm-1-mode),np.eye(self.Nq**self.Nd)])
        return H_cavity

    def getInteractionHamiltonian(self):
        z_op = self.ElectronPosition()
        H_interaction = np.zeros((np.power(self.Np+1, self.Nm) * np.power(self.Nq, self.Nd),
                             np.power(self.Np+1, self.Nm) * np.power(self.Nq, self.Nd)))
        for mode in range(self.Nm):
            gc = self.parameters_cavity[mode]['gc']
            for donor in range(self.Nd):
                H_interaction = H_interaction + 0.5*gc * \
                    kron([np.eye(mode), adagger(self.Np)+
                          a(self.Np), np.eye(self.Nm-1-mode), np.eye(
                        np.power(self.Nq, donor)), z_op[donor]+np.eye(4), np.eye(np.power(self.Nq, (self.Nd-1-donor)))])
        return H_interaction

    def getNoiseHamiltonian(self):
        z_op = self.ElectronPosition()
        Hnoise = []
        for noise in range(self.Nn):
            wn = self.parameters_noise[noise]['wn']
            e_donors = self.parameters_noise[noise]['effected_donors']
            hn = np.zeros((np.power(self.Np+1, self.Nm) * np.power(self.Nq, self.Nd),
                           np.power(self.Np+1, self.Nm) * np.power(self.Nq, self.Nd)))
            for donor in range(self.Nd):
                if(e_donors[donor]):
                    hn = 0.5*wn*kron([np.eye(np.power(self.Np+1, self.Nm)), np.eye(
                        np.power(self.Nq, donor)), z_op[donor], np.eye(np.power(self.Nq, (self.Nd-1-donor)))])
            Hnoise.append(hn)
        return Hnoise

    def getSystemHamiltonian(self):
        return self.getCavityHamiltonian() + self.getQubitHamiltonian() + self.getInteractionHamiltonian()

def computeDensityMatrixZeroNoise(system, initialState, T, Nt):
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
    # Oscillating term in L-F space. This guy is time dependent. We'll shape this as (N_t, ROW, COLS)
    t = np.arange(0, T, T/Nt)
    W_mat = np.kron(W, np.ones(Nstates))-np.kron(np.ones(Nstates), W)
    eiwtList = np.exp(np.kron(t[:, np.newaxis], W_mat)*1j)
    eiwt = np.zeros((Nt, Nstates**2, Nstates**2))
    for i in range(Nt):
        eiwt[i].flat[slice(0, None, 1+Nstates**2)] = eiwtList[i]
    rho = (np.matmul(np.matmul(np.matmul(RR, eiwt), RR_INV),
                     initial_rho_vec)).reshape(Nt, Nstates, Nstates)
    return t, rho

def partialTraceSystem(rho, Np, Nd, Nq=4, Nm=1):
    Nt = rho.shape[0]
    rho_tensor = rho.reshape(Nt, (Np+1)**Nm, Nq**Nd, (Np+1)**Nm, Nq**Nd)
    rho_ff = np.trace(rho_tensor, axis1=1, axis2=3)
    rho_cav = np.trace(rho_tensor, axis1=2, axis2=4)
    return rho_cav, rho_ff
