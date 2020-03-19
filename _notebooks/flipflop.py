import numpy as np
from _notebooks.constants import *

def FlipFlopEnergies(parameters):
    '''
        Returns the eigenenergies for the flip-flip qubit given above
        arguments:
            parameters{
                'Vt': tunnel coupling,
                'wB': zeeman splitting,
                'eps':applied electric field energy
            }
        returns:
            4x4 matrix with energies along the diagonal
    '''
    Vt = parameters['Vt']
    wB = parameters['wB']
    eps = parameters['eps']
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
    return E_ff


def ElectronPositionCoefficients(parameters):
    '''
        Returns the coefficients z_jk for the flip-flip qubit given above
        arguments:
            parameters{
                'Vt': tunnel coupling,
                'wB': zeeman splitting,
                'eps':applied electric field energy
            }
        returns:
            4x4 matrix
    '''
    Vt = parameters['Vt']
    wB = parameters['wB']
    eps = parameters['eps']
    w0 = np.sqrt(eps**2 + Vt**2)
    eta = np.arctan2(Vt, eps)
    Zcoef = np.zeros((4, 4))
    Zcoef[0, 1] = -(hyperfine*w0*Delta*wB*np.cos(eta)*np.sin(eta)**2)/(4*wB*(w0**2-wB**2))
    Zcoef[0, 3] = (hyperfine**2*w0**3*np.cos(eta) *
                   np.sin(eta)**2)/(4*wB*(w0**2-wB**2)**2)
    Zcoef[1, 0] = np.sin(eta) + (hyperfine*np.cos(eta)*np.sin(eta))/(4*w0)
    Zcoef[1, 1] = -(hyperfine*w0*np.cos(eta)*np.sin(eta))/(2*(w0**2-wB**2))
    Zcoef[1, 3] = -(Delta*wB*np.cos(eta)*np.sin(eta))/(2*w0)
    Zcoef[2, 2] = -(hyperfine*w0**2*np.cos(eta)*np.sin(eta))/(2*wB*(w0**2-wB**2))
    Zcoef[3, 0] = np.cos(eta) - (hyperfine*np.sin(eta)**2)/(4*w0)
    Zcoef[3, 1] = (hyperfine*w0*np.sin(eta)**2)/(2*(w0**2-wB**2))
    Zcoef[3, 3] = (Delta*wB*np.sin(eta)**2)/(2*w0)
    return Zcoef
