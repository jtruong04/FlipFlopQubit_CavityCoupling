from Scripts.constants import *
from Scripts.pauli import *

from Scripts.Noise.Pink import *


class CavityCoupled_4:
    def __init__(self, params={
        'parameters_cavity': {'n0': 4, 'w_cav': 10},
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

        self.df = self.w_f - self.w_cav
        self.dc = self.w_c - self.w_cav

        self.Nstates = 1 if self.n0 == 0 else (3 if self.n0 == 1 else 4)
        self.Nn = 2

    def getNoiseHamiltonian(self, approxNoise=None):
        hL = 0
        hR = 0
        n0 = self.n0
        gf = self.g_f
        gc = self.g_c
        df = self.df
        dc = self.dc
        wp = self.wp
        w03 = self.w03
        w30 = self.w30
        w33 = self.w33
        dzz = (self.g_c**2)*(w30+w33)/(self.dc**2)
        dxy = gc*gf*wp / (dc*(dc-df))
        if self.n0 == 0:
            hL = np.array([[w03+w30+w33]])
            hR = np.array([[w03+w30+w33]])
        elif self.n0 == 1:
            hL = np.array([
                [w03+w30+w33-2*dzz*n0,           0, -gc*wp/dc],
                [0, w03+w30+w33,               -dxy*n0],
                [-gc*wp/dc,     -dxy*n0, -w03+w30-w33-2*dxy*n0]
            ])
            hR = np.array([
                [w03+w30+w33-2*dzz*n0, -gc*wp/dc,           0],
                [-gc*wp/dc, -w03+w30-w33-2*dxy*n0,     -dxy*n0],
                [0,               -dxy*n0, w03+w30+w33]
            ])
        elif self.n0 == 2:
            hL = np.array([
                [w03+w30+w33-2*dzz*n0,                        0,          -
                    gc*np.sqrt(n0)*wp/dc,                         0],
                [0, w03+w30+w33-2*dzz,                            -dxy,  -gc*wp/dc],
                [-gc*np.sqrt(n0)*wp/dc,                  -dxy, -
                 w03+w30-w33-2*dzz-2*dxy,                         0],
                [0, -gc*wp/dc,
                 0, -w03+w30-w33-2*dxy]
            ])
            hR = np.array([
                [w03+w30+w33-2*dzz*n0,          -gc *
                    np.sqrt(n0)*wp/dc,                        0,                         0],
                [-gc*np.sqrt(n0)*wp/dc, -w03+w30-w33-2*dzz-2 *
                 dxy,                  -dxy,                         0],
                [0,                            -dxy,
                    w03+w30+w33-2*dzz,  -gc*wp/dc],
                [0,                               0, -
                 gc*wp/dc, -w03+w30-w33-2*dxy]
            ])
        else:
            hL = np.array([
                [w03+w30+w33-2*dzz*n0,                        0,          -
                    gc*np.sqrt(n0)*wp/dc,                                      0],
                [0, w03+w30+w33-2*dzz*(n0-1),                        -
                 dxy,               -gc*np.sqrt(n0-1)*wp/dc],
                [-gc*np.sqrt(n0)*wp/dc,                  -dxy, -w03+w30-w33 -
                 2*dzz*(n0-1)-2*dxy,                                      0],
                [0, -gc*np.sqrt(n0-1)*wp/dc,
                 0, -w03+w30-w33-2*dxy*(n0-1)-2*dzz*(n0-1)]
            ])
            hR = np.array([
                [w03+w30+w33-2*dzz*n0,          -gc*np.sqrt(
                    n0)*wp/dc,                        0,                                      0],
                [-gc*np.sqrt(n0)*wp/dc, -w03+w30-w33-2*dzz*(
                    n0-1)-2*dxy,                  -dxy,                                      0],
                [0,                        -dxy, w03+w30+w33-2*dzz *
                    (n0-1),               -gc*np.sqrt(n0-1)*wp/dc],
                [0,                               0, -gc *
                 np.sqrt(n0-1)*wp/dc, -w03+w30-w33-2*dxy*(n0-1)-2*dzz*(n0-1)]
            ])
        # print(hL)
        return [hL, hR], [{'type': Pink()}, {'type': Pink()}]

    def getSystemHamiltonian(self, approxInteraction=None):
        dd = 2*self.g_c*self.g_c/self.dc
        if self.n0 == 0:
            return np.array([[0]])
        elif self.n0 == 1:
            return np.array([[-self.df-dd, self.g_f, self.g_f],
                             [self.g_f,        0,        0],
                             [self.g_f,        0,        0]
                             ])
        elif self.n0 == 2:
            return np.array([[-self.df-dd, self.g_f*np.sqrt(self.n0), self.g_f*np.sqrt(self.n0),          0],
                             [self.g_f * np.sqrt(self.n0),
                              0,                         0,   self.g_f],
                             [self.g_f * np.sqrt(self.n0),
                              0,                         0,   self.g_f],
                             [0,                  self.g_f,
                                 self.g_f, self.df+dd]
                             ])
        else:
            return np.array([[-self.df-dd,   self.g_f*np.sqrt(self.n0),   self.g_f*np.sqrt(self.n0),                           0],
                             [self.g_f * np.sqrt(self.n0),                           0,
                              0, self.g_f*np.sqrt(self.n0-1)],
                             [self.g_f * np.sqrt(self.n0),                           0,
                              0, self.g_f*np.sqrt(self.n0-1)],
                             [0, self.g_f*np.sqrt(self.n0-1), self.g_f *
                              np.sqrt(self.n0-1),                  self.df+dd]
                             ])
        return
