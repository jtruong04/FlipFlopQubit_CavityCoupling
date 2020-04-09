from Scripts.constants import *
from Scripts.pauli import *

from Scripts.Noise.Pink import *

class FlipFlopSystem:
    def __init__(self, params = {
                                    'Nm': 1,
                                    'Np': 1,
                                    'Nd': 1,
                                    'Nn': 0,
                                    'parameters_cavity': [{'wc':0, 'gc':0}],
                                    'parameters_qubits': [{'eps':0, 'wB':0, 'Vt':0}],
                                    'parameters_noise': [{'type':Pink(),'wn':0,'effected_donors':[0]}]
                                }
                ):
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
        self.Nd = params['Nd']
        self.Np = params['Np']
        self.Nm = params['Nm']
        self.Nn = params['Nn']
        self.Nstates = (self.Np + 1)**self.Nm * (self.Nq**self.Nd)
        self.parameters_qubits = params['parameters_qubits']
        self.parameters_cavity = params['parameters_cavity']
        self.parameters_noise = params['parameters_noise']
        self.Z1q = []
        self.Zcoefs = []
        self.E1q = []
        assert(len(self.parameters_qubits) == self.Nd)
        assert(len(self.parameters_cavity) == self.Nm)
        assert(len(self.parameters_noise) == self.Nn)

    def Energies(self):
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
            Zcoef[1, 0] = np.sin(
                eta) + (hyperfine*np.cos(eta)*np.sin(eta))/(4*w0)
            Zcoef[1, 1] = -(hyperfine*w0*np.cos(eta) *
                            np.sin(eta))/(2*(w0**2-wB**2))
            Zcoef[1, 3] = -(Delta*wB*np.cos(eta)*np.sin(eta))/(2*w0)
            Zcoef[2, 2] = -(hyperfine*w0**2*np.cos(eta) *
                            np.sin(eta))/(2*wB*(w0**2-wB**2))
            Zcoef[3, 0] = np.cos(eta) - (hyperfine*np.sin(eta)**2)/(4*w0)
            Zcoef[3, 1] = (hyperfine*w0*np.sin(eta)**2)/(2*(w0**2-wB**2))
            Zcoef[3, 3] = (Delta*wB*np.sin(eta)**2)/(2*w0)
            # print(f'Donor {i} Zcoef:')
            # print(Zcoef)
            z = np.zeros((self.Nq, self.Nq))
            for j in np.arange(self.Nq):
                for k in np.arange(self.Nq):
                    z = z + np.kron(pauli[j], pauli[k])*Zcoef[j, k]
            self.Z1q.append(z)
            self.Zcoefs.append(Zcoef)
        return self.Z1q, self.Zcoefs

    def getQubitHamiltonian(self):
        H_qubit = np.zeros((np.power(self.Np+1, self.Nm) * np.power(self.Nq, self.Nd),
                            np.power(self.Np+1, self.Nm) * np.power(self.Nq, self.Nd)))
        ff_energies = self.Energies()
        for donor in range(self.Nd):
            H_qubit = H_qubit + kron([np.eye(np.power(self.Np+1, self.Nm)), np.eye(
                np.power(self.Nq, donor)), ff_energies[donor], np.eye(np.power(self.Nq, (self.Nd-1-donor)))])
        return H_qubit

    def getCavityHamiltonian(self):
        H_cavity = np.zeros((np.power(self.Np+1, self.Nm) * np.power(self.Nq, self.Nd),
                             np.power(self.Np+1, self.Nm) * np.power(self.Nq, self.Nd)))
        for mode in range(self.Nm):
            wc = self.parameters_cavity[mode]['wc']
            # print(wc)
            H_cavity = H_cavity + kron([np.eye(mode), np.matmul(adagger(self.Np), a(
                self.Np)), np.eye(self.Nm-1-mode), np.eye(self.Nq**self.Nd)])*wc
        # print(H_cavity)
        return H_cavity

    def getInteractionHamiltonian(self, approxInteraction=[1, 1, 1, 1, 1, 1, 1, 1, 1]):
        z_op, z_coefs = self.ElectronPosition()
        for donor in range(self.Nd):
            z_op[donor] = (
                z_coefs[donor][0, 1]*np.kron(pauli[0], pauli[1])*approxInteraction[0] +
                z_coefs[donor][0, 3]*np.kron(pauli[0], pauli[3])*approxInteraction[1] +
                z_coefs[donor][1, 0]*np.kron(pauli[1], pauli[0])*approxInteraction[2] +
                z_coefs[donor][1, 1]*np.kron(pauli[1], pauli[1])*approxInteraction[3] +
                z_coefs[donor][1, 3]*np.kron(pauli[1], pauli[3])*approxInteraction[4] +
                z_coefs[donor][2, 2]*np.kron(pauli[2], pauli[2])*approxInteraction[5] +
                z_coefs[donor][3, 0]*np.kron(pauli[3], pauli[0])*approxInteraction[6] +
                z_coefs[donor][3, 1]*np.kron(pauli[3], pauli[1])*approxInteraction[7] +
                z_coefs[donor][3, 3]*np.kron(pauli[3], pauli[3])*approxInteraction[8]
            )
        H_interaction = np.zeros((np.power(self.Np+1, self.Nm) * np.power(self.Nq, self.Nd),
                                  np.power(self.Np+1, self.Nm) * np.power(self.Nq, self.Nd)))
        for mode in range(self.Nm):
            gc = self.parameters_cavity[mode]['gc']
            for donor in range(self.Nd):
                H_interaction = H_interaction + 0.5*gc * \
                    kron([np.eye(mode), adagger(self.Np) +
                          a(self.Np), np.eye(self.Nm-1-mode), np.eye(
                        np.power(self.Nq, donor)), z_op[donor]+np.eye(4), np.eye(np.power(self.Nq, (self.Nd-1-donor)))])
        return H_interaction

    def getNoiseHamiltonian(self, approxNoise=[1, 1, 1, 1, 1, 1, 1, 1, 1]):
        z_op, z_coefs = self.ElectronPosition()
        for donor in range(self.Nd):
            z_op[donor] = (
                z_coefs[donor][0, 1]*np.kron(pauli[0], pauli[1])*approxNoise[0] +
                z_coefs[donor][0, 3]*np.kron(pauli[0], pauli[3])*approxNoise[1] +
                z_coefs[donor][1, 0]*np.kron(pauli[1], pauli[0])*approxNoise[2] +
                z_coefs[donor][1, 1]*np.kron(pauli[1], pauli[1])*approxNoise[3] +
                z_coefs[donor][1, 3]*np.kron(pauli[1], pauli[3])*approxNoise[4] +
                z_coefs[donor][2, 2]*np.kron(pauli[2], pauli[2])*approxNoise[5] +
                z_coefs[donor][3, 0]*np.kron(pauli[3], pauli[0])*approxNoise[6] +
                z_coefs[donor][3, 1]*np.kron(pauli[3], pauli[1])*approxNoise[7] +
                z_coefs[donor][3, 3]*np.kron(pauli[3], pauli[3])*approxNoise[8]
            )
        Hnoise = []
        for noise in range(self.Nn):
            wn = self.parameters_noise[noise]['wn']
            e_donors = self.parameters_noise[noise]['effected_donors']
            hn = np.zeros((np.power(self.Np+1, self.Nm) * np.power(self.Nq, self.Nd),
                           np.power(self.Np+1, self.Nm) * np.power(self.Nq, self.Nd)))
            # print(e_donors)
            for donor in range(self.Nd):
                if(e_donors[donor]==1):
                    hn = 0.5*wn*kron([np.eye(np.power(self.Np+1, self.Nm)), np.eye(
                        np.power(self.Nq, donor)), z_op[donor], np.eye(np.power(self.Nq, (self.Nd-1-donor)))])
            Hnoise.append(hn)
        return Hnoise, self.parameters_noise

    def getSystemHamiltonian(self, approxInteraction=[1, 1, 1, 1, 1, 1, 1, 1, 1]):
        H = self.getCavityHamiltonian() + self.getQubitHamiltonian() + \
            self.getInteractionHamiltonian(approxInteraction)
        # print(H)
        return H #self.getCavityHamiltonian() + self.getQubitHamiltonian() + self.getInteractionHamiltonian(approx)
