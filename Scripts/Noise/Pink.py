# import scipy
import numpy as np
from scipy.special import sici
from Scripts.Noise.Noise import *

def si(x):
    return sici(x)[0]
def ci(x):
    return sici(x)[1]

# Pink noise specific functions
def Ciw(t, w, wl=2.0*np.pi*1e-9, wh=2.0*np.pi*1e3):
    return 2.0*(ci(np.abs(w*t+wh*t))+ci(np.abs(w*t-wh*t))-ci(np.abs(w*t+wl*t))-ci(np.abs(w*t-wl*t)))

def Siw(t, w, wl=2.0*np.pi*1e-9, wh=2.0*np.pi*1e3):
    return 2.0*(si(w*t+wh*t) + si(w*t-wh*t)-si(w*t+wl*t)-si(w*t-wl*t))


class Pink(Noise):
    def __init__(self, wl=2.0*np.pi*1e-9, wh=2.0*np.pi*1e3):
        self.wl = wl
        self.wh = wh

    def Sw(self,w):
        S = np.sqrt(np.pi/2)/np.log(self.wh/self.wl)
        if( np.abs(w) > self.wl and np.abs(w) < self.wh):
            return S/np.abs(w)
        else:
            return 0

    def St(self,t):
        return (ci(self.wh*t)-ci(self.wl*t)) / np.log(self.wh/self.wl)

    def f0(self,t):
        return t*self.St(t) - (np.sin(self.wh*t)/self.wh - np.sin(self.wl*t)/self.wl)/np.log(self.wh/self.wl)

    def ft(self,t):
        return 0.5*t*t*self.St(t) - t*(np.sin(self.wh*t)/self.wh - np.sin(self.wl*t)/self.wl)/(2.0*np.log(self.wh/self.wl)) + ((1-np.cos(self.wh*t))/self.wh**2 - (1-np.cos(self.wl*t))/self.wl**2)/np.log(self.wh/self.wl)

    def fsin(self,t, w):
        return -np.cos(w*t)*self.St(t)/w + (Ciw(t,w)-2*np.log(np.abs(((self.wh**2-w**2)/(self.wh**2))*((self.wl**2)/(self.wl**2-w**2)))))/(4*w*np.log(self.wh/self.wl))

    def fcos(self,t, w):
        return np.sin(w*t)*self.St(t)/w-Siw(t,w)/(4*np.log(self.wh/self.wl))

    def ftsin(self,t, w):
        return -self.St(t)*np.cos(w*t)*t/w+(self.St(t)*np.sin(w*t)/w-Siw(t, w)/(4*np.log(self.wh/self.wl)))/w+(((np.sin(w*t+self.wh*t))/(w+self.wh))+((np.sin(w*t-self.wh*t))/(w-self.wh))-((np.sin(w*t+self.wl*t))/(w+self.wl))-((np.sin(w*t-self.wl*t))/(w-self.wl)))/(2*w*np.log(self.wh/self.wl))

    def ftcos(self,t, w):
        return self.St(t)*np.sin(w*t)*t/w+(self.St(t)*np.cos(w*t)/w-(Ciw(t, w)-2*np.log(np.abs(((self.wh**2-w**2)/(self.wh**2))*((self.wl**2)/(self.wl**2-w**2)))))/(4*w*np.log(self.wh/self.wl)))/w - (((1-np.cos(w*t+self.wh*t))/(w+self.wh))+((1-np.cos(w*t-self.wh*t))/(w-self.wh))-((1-np.cos(w*t+self.wl*t))/(w+self.wl))-((1-np.cos(w*t-self.wl*t))/(w-self.wl)))/(2*w*np.log(self.wh/self.wl))
