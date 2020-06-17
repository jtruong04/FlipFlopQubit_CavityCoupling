from Scripts.constants import *
from Scripts.pauli import *

from Scripts.Noise.Pink import *

class CavityCoupled_16:
    def __init__(self, params={
                                'parameters_cavity': {'n0':4, 'w_cav': 10},
                                'parameters_qubits': {'w_c': 11.4, 'w_f': 10.8, 'g_c': 0.1, 'g_f': 0.01},
                                'parameters_noise': {'wp': 0, 'w03': 0, 'w30': 0, 'w33': 0}
                              }
                ):
        self.n0 = params['parameters_cavity']['n0']
        self.w_cav = params['parameters_cavity']['w_cav']
        self.w_c = params['parameters_qubits']['w_c']
        self.w_f = params['parameters_qubits']['w_f']
        self.g_c = params['parameters_qubits']['g_c']
        self.g_f = params['parameters_qubits']['g_f']
        self.wp = params['parameters_noise']['wp']
        self.w03 = params['parameters_noise']['w03']
        self.w30 = params['parameters_noise']['w30']
        self.w33 = params['parameters_noise']['w33']

        self.Nstates = 1 if self.n0 == 0 else (5 if self.n0 == 1 else (11 if self.n0 == 2 else (15 if self.n0 == 3 else 16)))
        self.states = [0] if self.n0 == 0 else ([0,1,2,4,8] if self.n0 == 1 else (
            [0, 1, 2, 4, 8, 3, 5, 6, 9, 10, 12] if self.n0 == 2 else ([0, 1, 2, 4, 8, 3, 5, 6, 9, 10, 12, 7, 11, 13, 14] if self.n0 == 3 else [0, 1, 2, 4, 8, 3, 5, 6, 9, 10, 12, 7, 11, 13, 14,15])))
        self.Nn = 2


    def getNoiseHamiltonian(self, approxNoise = None):

        w03 = self.w03
        w30 = self.w30
        w33 = self.w33
        wp = self.wp

        HL_Diag = np.diag([
             w03+w30+w33,  w30+w30+w33,  w03+w30+w33,  w03+w30+w33, 
            -w03+w30-w33, -w03+w30-w33, -w03+w30-w33, -w03+w30-w33,
             w03-w30-w33,  w03-w30-w33,  w03-w30-w33,  w03-w30-w33,
            -w03-w30+w33, -w03-w30+w33, -w03-w30+w33, -w03-w30+w33 
        ])
        HR_Diag = np.diag([
            w03+w30+w33, -w03+w30-w33, w03-w30-w33, -w03-w30+w33,
            w03+w30+w33, -w03+w30-w33, w03-w30-w33, -w03-w30+w33,
            w03+w30+w33, -w03+w30-w33, w03-w30-w33, -w03-w30+w33,
            w03+w30+w33, -w03+w30-w33, w03-w30-w33, -w03-w30+w33
        ])
        HL_Off = np.zeros((16, 16))
        HR_Off = np.zeros((16, 16))

        HL_Off[4,8] = wp
        HL_Off[5,9] = wp
        HL_Off[6,10] = wp
        HL_Off[7,11] = wp
        HR_Off[1,2] = wp
        HR_Off[5, 6] = wp
        HR_Off[9, 10] = wp
        HR_Off[13, 14] = wp

        # print([(HL_Diag+HL_Off+np.transpose(HL_Off))[self.states,:][:,self.states],
            #    (HR_Diag+HR_Off+np.transpose(HR_Off))[self.states,:][:,self.states]])
        return [(HL_Diag+HL_Off+np.transpose(HL_Off))[self.states,:][:,self.states], (HR_Diag+HR_Off+np.transpose(HR_Off))[self.states,:][:,self.states]], [{'type': Pink()}, {'type': Pink()}]

    def getSystemHamiltonian(self, approxInteraction = None):
        H0 = np.diag([
                        0 +        0+        0, -  self.w_cav +   self.w_f +        0,   -self.w_cav +        0 +   self.w_c, -2*self.w_cav +   self.w_f +   self.w_c,
              -self.w_cav + self.w_f +       0, -2*self.w_cav + 2*self.w_f +        0, -2*self.w_cav + self.w_f +   self.w_c, -3*self.w_cav + 2*self.w_f +   self.w_c,
              -self.w_cav +        0 +self.w_c, -2*self.w_cav +   self.w_f + self.w_c, -2*self.w_cav +        0 + 2*self.w_c, -3*self.w_cav +   self.w_f + 2*self.w_c,
            -2*self.w_cav + self.w_f +self.w_c, -3*self.w_cav + 2*self.w_f + self.w_c, -3*self.w_cav + self.w_f + 2*self.w_c, -4*self.w_cav + 2*self.w_f + 2*self.w_c
        ])
        Hc = np.zeros((16, 16))
        Hf = np.zeros((16, 16))

        # Charge coupling
        Hc[ 0,  2] = np.sqrt(self.n0-0) if self.n0 >= 0 else 0
        Hc[ 1,  3] = np.sqrt(self.n0-1) if self.n0 >= 1 else 0
        Hc[ 0,  8] = np.sqrt(self.n0-0) if self.n0 >= 0 else 0
        Hc[ 1,  9] = np.sqrt(self.n0-1) if self.n0 >= 1 else 0
        Hc[ 2, 10] = np.sqrt(self.n0-1) if self.n0 >= 1 else 0
        Hc[ 3, 11] = np.sqrt(self.n0-2) if self.n0 >= 2 else 0
        Hc[ 4,  6] = np.sqrt(self.n0-1) if self.n0 >= 1 else 0
        Hc[ 5,  7] = np.sqrt(self.n0-2) if self.n0 >= 2 else 0
        Hc[ 4, 12] = np.sqrt(self.n0-1) if self.n0 >= 1 else 0
        Hc[ 5, 13] = np.sqrt(self.n0-2) if self.n0 >= 2 else 0
        Hc[ 6, 14] = np.sqrt(self.n0-2) if self.n0 >= 2 else 0
        Hc[ 7, 15] = np.sqrt(self.n0-3) if self.n0 >= 3 else 0
        Hc[ 8, 10] = np.sqrt(self.n0-1) if self.n0 >= 1 else 0
        Hc[ 9, 11] = np.sqrt(self.n0-2) if self.n0 >= 2 else 0
        Hc[12, 14] = np.sqrt(self.n0-2) if self.n0 >= 2 else 0
        Hc[13, 15] = np.sqrt(self.n0-3) if self.n0 >= 3 else 0

        # Flip Flop Coupling
        Hf[ 0,  1] = np.sqrt(self.n0-0) if self.n0 >= 0 else 0
        Hf[ 2,  3] = np.sqrt(self.n0-1) if self.n0 >= 1 else 0
        Hf[ 0,  4] = np.sqrt(self.n0-0) if self.n0 >= 0 else 0
        Hf[ 1,  5] = np.sqrt(self.n0-1) if self.n0 >= 1 else 0
        Hf[ 2,  6] = np.sqrt(self.n0-2) if self.n0 >= 1 else 0
        Hf[ 3,  7] = np.sqrt(self.n0-1) if self.n0 >= 2 else 0
        Hf[ 4,  5] = np.sqrt(self.n0-1) if self.n0 >= 1 else 0
        Hf[ 6,  7] = np.sqrt(self.n0-2) if self.n0 >= 2 else 0
        Hf[ 8,  9] = np.sqrt(self.n0-1) if self.n0 >= 1 else 0
        Hf[10, 11] = np.sqrt(self.n0-2) if self.n0 >= 2 else 0
        Hf[ 8, 12] = np.sqrt(self.n0-1) if self.n0 >= 1 else 0
        Hf[ 9, 13] = np.sqrt(self.n0-2) if self.n0 >= 2 else 0
        Hf[10, 14] = np.sqrt(self.n0-2) if self.n0 >= 2 else 0
        Hf[11, 15] = np.sqrt(self.n0-3) if self.n0 >= 3 else 0
        Hf[12, 13] = np.sqrt(self.n0-2) if self.n0 >= 2 else 0
        Hf[14, 15] = np.sqrt(self.n0-3) if self.n0 >= 3 else 0
        # print((H0 + self.g_c * (Hc + np.transpose(Hc)) + self.g_f *
        #        (Hf + np.transpose(Hf))))
        # print((H0 + self.g_c * (Hc + np.transpose(Hc)) + self.g_f *
        #        (Hf + np.transpose(Hf)))[self.states, :][:, self.states])
        return (H0 + self.g_c * (Hc + np.transpose(Hc)) + self.g_f * (Hf + np.transpose(Hf)))[self.states,:][:,self.states]
