from Scripts.constants import *
from Scripts.pauli import *

from Scripts.Noise.Pink import *


class FlipFlopSWSystem:
    def __init__(self, params={
                            'n0': 4,
                            'parameters_cavity': {'wcav': 0, 'g': 0},
                            'parameters_qubits': [{'eps': 0, 'wB': 0, 'Vt': 0}, {'eps': 0, 'wB': 0, 'Vt': 0}],
                            'parameters_noise': [{'type': Pink(), 'wn': 0, 'effected_donors': [1,0]}, {'type': Pink(), 'wn': 0, 'effected_donors': [0,1]}]
                              }
    ):
        '''
        Units being used for analysis
            Energy: GHz
            Time  : ns
        Arguments:
            n0 : block >= 2
            parameters_qubits : list of size Nd dicts containing tunable parameters
                {'Vt','wB','eps'}
            parameters_cavity : list of size Nm dicts containing cavity parameters
                {'wcav', 'g'}
            parameters_noise : list of size Nn
                {'type', 'wn', 'effected_donors'}
                type: 'pink', 'white', etc.
                wn: strength of noise
                effected_donors: boolean list of len Nd
        '''
        self.parameters_qubits = params['parameters_qubits']
        self.parameters_cavity = params['parameters_cavity']
        self.parameters_noise = params['parameters_noise']
        self.Nstates=4
        self.Nn = 2
        self.n0 = params['n0']
        self.Zcoefs = []
        self.E1q = []

        for i in range(2):
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
            self.Zcoefs.append(Zcoef)
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

        for i in range(2):
            self.parameters_qubits[i]['gc'] = 0.5*self.parameters_cavity['g']*self.Zcoefs[i][1,0]
            self.parameters_qubits[i]['gf'] = 0.5*self.parameters_cavity['g']*self.Zcoefs[i][3,1]
            self.parameters_qubits[i]['dc'] = (self.E1q[i][2, 2]-self.E1q[i][0, 0])-self.parameters_cavity['wcav']
            self.parameters_qubits[i]['df'] = (self.E1q[i][1, 1]-self.E1q[i][0, 0])-self.parameters_cavity['wcav']

    def getNoiseHamiltonian(self, approxNoise=None):
        HnL = np.zeros((4, 4))
        HnR = np.zeros((4, 4))
        gcL = self.parameters_qubits[0]['gc']
        gcR = self.parameters_qubits[1]['gc']
        gfL = self.parameters_qubits[0]['gf']
        gfR = self.parameters_qubits[1]['gf']
        dcL = self.parameters_qubits[0]['dc']
        dcR = self.parameters_qubits[1]['dc']
        dfL = self.parameters_qubits[0]['df']
        dfR = self.parameters_qubits[1]['df']
        wnL = self.parameters_noise[0]['wn']
        wnR = self.parameters_noise[1]['wn']
        z11L = self.Zcoefs[0][1,1]
        z22L = self.Zcoefs[0][2,2]
        z11R = self.Zcoefs[1][1, 1]
        z22R = self.Zcoefs[1][2, 2]

        HnL[0, 0] = 0
        HnL[0, 1] = 0
        HnL[0, 2] = -gcL*wnL*((z11L-z22L)*np.sqrt(self.n0-1)*0+(z11L+z22L)*np.sqrt(self.n0))/(2.0*dcL)
        HnL[0, 3] = 0 #gcL*gfR*wnL*(z11L-z22L)/(2.0*dcL*(dfR-dcL))
        HnL[1, 0] = HnL[0,1]
        HnL[1, 1] = 0
        HnL[1, 2] = 0 #gcL*gfR*wnL*(z11L+z22L)/(2.0*dcL*(dfR-dcL))
        HnL[1, 3] = gcL*wnL*((-z11L+z22L)*np.sqrt(self.n0-2)*0-(z11L+z22L)*np.sqrt(self.n0-1))/(2.0*dcL)
        HnL[2, 0] = HnL[0,2]
        HnL[2, 1] = HnL[1,2]
        HnL[2, 2] = 0 #-gcL*gfL*wnL*(z11L+z22L)/(dcL*(dcL-dfL))
        HnL[2, 3] = 0
        HnL[3, 0] = HnL[0,3]
        HnL[3, 1] = HnL[1,3]
        HnL[3, 2] = HnL[2,3]
        HnL[3, 3] = 0 #-gcL*gfL*wnL*(z11L+z22L)/(dcL*(dcL-dfL))

        HnR[0, 0] = 0
        HnR[0, 1] = -gcR*wnR*((z11R-z22R)*np.sqrt(self.n0-1)*0+(z11R+z22R)*np.sqrt(self.n0))/(2.0*dcR)
        HnR[0, 2] = 0
        HnR[0, 3] = 0 # gcR*gfL*wnR*(z11R-z22R)/(2.0*dcR*(dfL-dcR))
        HnR[1, 0] = HnR[0, 1]
        HnR[1, 1] = 0 #-gcR*gfR*wnR*(z11R+z22R)/(dcR*(dcR-dfR))
        HnR[1, 2] = 0 #gcR*gfL*wnR*(z11R+z22R)/(2.0*dcL*(dfL-dcR))
        HnR[1, 3] = 0
        HnR[2, 0] = HnR[0, 2]
        HnR[2, 1] = HnR[1, 2]
        HnR[2, 2] = 0
        HnR[2, 3] = gcR*wnR*((-z11R+z22R)*np.sqrt(self.n0-2)*0-(z11R+z22R)*np.sqrt(self.n0-1))/(2.0*dcR)
        HnR[3, 0] = HnR[0, 3]
        HnR[3, 1] = HnR[1, 3]
        HnR[3, 2] = HnR[2, 3]
        HnR[3, 3] = 0 #-gcR*gfR*wnR*(z11R+z22R)/(dcR*(dcR-dfR))

        return [HnL,HnR], self.parameters_noise

    def getSystemHamiltonian(self, approxInteraction=None):
        H = np.zeros((4,4))
        gcL = self.parameters_qubits[0]['gc']
        gcR = self.parameters_qubits[1]['gc']
        gfL = self.parameters_qubits[0]['gf']
        gfR = self.parameters_qubits[1]['gf']
        dcL = self.parameters_qubits[0]['dc']
        dcR = self.parameters_qubits[1]['dc']
        dfL = self.parameters_qubits[0]['df']
        dfR = self.parameters_qubits[1]['df']

        H[0, 0] = -gcL*gcL/dcL - gcR*gcR/dcR - dfL/2.0 - dfR/2.0
        H[0, 1] = gfR*np.sqrt(self.n0)
        H[0, 2] = gfL*np.sqrt(self.n0)
        H[0, 3] = 0

        H[1, 0] = H[0,1]
        H[1, 1] = -dfL/2.0 + dfR/2.0
        H[1, 2] = 0
        H[1, 3] = gfL*np.sqrt(self.n0-1)

        H[2, 0] = H[0,2]
        H[2, 1] = 0
        H[2, 2] = dfL/2.0 - dfR/2.0
        H[2, 3] = gfR*np.sqrt(self.n0-1)

        H[3, 0] = 0
        H[3, 1] = H[1,3]
        H[3, 2] = H[2,3]
        H[3, 3] = gcL*gcL/dcL + gcR*gcR/dcR + dfL/2.0 + dfR/2.0

        # print(H)
        return H  # self.getCavityHamiltonian() + self.getQubitHamiltonian() + self.getInteractionHamiltonian(approx)
