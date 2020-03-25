# import scipy
import numpy as np
from scipy.special import sici

def si(x):
    return sici(x)[0]
def ci(x):
    return sici(x)[1]

# Pink noise specific functions

def Sw(w, param={'wl': 2.0*np.pi*1e-9, 'wh': 2.0*np.pi*1e3}):
    wl = param['wl']
    wh = param['wh']
    S = np.sqrt(np.pi/2)/np.log(wh/wl)
    if( np.abs(w) > wl and np.abs(w) < wh):
        return S/np.abs(w)
    else:
        return 0

def St(t, param={'wl': 2.0*np.pi*1e-9, 'wh': 2.0*np.pi*1e3}):
    wl = param['wl']
    wh = param['wh']
    return (ci(wh*t)-ci(wl*t)) / np.log(wh/wl)

def Ciw(t,w, param={'wl': 2.0*np.pi*1e-9, 'wh': 2.0*np.pi*1e3}):
    wl = param['wl']
    wh = param['wh']
    return 2.0*(ci(np.abs(w*t+wh*t))+ci(np.abs(w*t-wh*t))-ci(np.abs(w*t+wl*t))-ci(np.abs(w*t-wl*t)))

def Siw(t,w, param={'wl': 2.0*np.pi*1e-9, 'wh': 2.0*np.pi*1e3}):
    wl = param['wl']
    wh = param['wh']
    return 2.0*(si(w*t+wh*t) + si(w*t-wh*t)-si(w*t+wl*t)-si(w*t-wl*t))

def f0(t, param={'wl': 2.0*np.pi*1e-9, 'wh': 2.0*np.pi*1e3}):
    wl = param['wl']
    wh = param['wh']
    return t*St(t,param) - (np.sin(wh*t)/wh - np.sin(wl*t)/wl)/np.log(wh/wl)

def ft(t, param={'wl': 2.0*np.pi*1e-9, 'wh': 2.0*np.pi*1e3}):
    wl = param['wl']
    wh = param['wh']
    0.5*t*t*St(t, param) - t*(np.sin(wh*t)/wh - np.sin(wl*t)/wl)/(2.0*np.log(wh.wl)) + ((1-np.cos(wh*t))/wh**2 - (1-np.cos(wl*t))/wl**2)/np.log(wh/wl)

def fsin(t, w, param={'wl': 2.0*np.pi*1e-9, 'wh': 2.0*np.pi*1e3}):
    wl = param['wl']
    wh = param['wh']
    return -np.cos(w*t)*St(t,param)/w + (Ciw(t,w,param)-2*np.log(np.abs(((wh**2-w**2)/(wh**2))*((wl**2)/(wl**2-w**2)))))/(4*w*np.log(wh/wl))


def fcos(t, w, param={'wl': 2.0*np.pi*1e-9, 'wh': 2.0*np.pi*1e3}):
    wl = param['wl']
    wh = param['wh']
    return np.sin(w*t)*St(t,param)/w-Siw(t,w,param)/(2*np.log(wh/wl))


def ftsin(t, w, param={'wl': 2.0*np.pi*1e-9, 'wh': 2.0*np.pi*1e3}):
    wl = param['wl']
    wh = param['wh']
    return -St(t, param)*np.cos(w*t)*t/w+(St(t, param)*np.sin(w*t)/w-Siw(t, w, param)/(4*np.log(wh/wl)))/w+(((np.sin(w*t+wh*t))/(w+wh))+((np.sin(w*t-wh*t))/(w-wh))-((np.sin(w*t+wl*t))/(w+wl))-((np.sin(w*t-wl*t))/(w-wl)))/(2*w*np.log(wh.wl))

def ftcos(t, w, param={'wl': 2.0*np.pi*1e-9, 'wh': 2.0*np.pi*1e3}):
    wl = param['wl']
    wh = param['wh']
    return St(t, param)*np.sin(w*t)*t/w+(St(t, param)*np.cos(w*t)/w-(Ciw(t, w, param)-2*np.log(np.abs(((wh**2-w**2)/(wh**2))*((wl**2)/(wl**2-w**2)))))/(4*w*np.log(wh/wl)))/w - (((1-np.cos(w*t+wh*t))/(w+wh))+((1-np.cos(w*t-wh*t))/(w-wh))-((1-np.cos(w*t+wl*t))/(w+wl))-((1-np.cos(w*t-wl*t))/(w-wl)))/(2*w*np.log(wh/wl))

# General functions for all noise below

def JS(t, w1, w2, param):
    if(w1 == 0 and w2 == 0):
        return -(ft(t,param)-t*f0(t,param))
    elif (w1 == -w2 and w2 != 0):
        return -(ftcos(t,w1,param)-t*fcos(t,w1,param))
    elif (w1!=0 and w2==0):
        return -(np.exp(1j*w1*t/2)/w1)*(np.cos(w1*t/2)*fsin(t,w1,param)-np.sin(w1*t/2)*(fcos(t,w1,param)+f0(t,param)))
    elif (w1==0 and w2!=0):
        return -(np.exp(1j*w2*t/2)/w2)*(np.cos(w2*t/2)*fsin(t, w2, param)-np.sin(w2*t/2)*(fcos(t, w2, param)+f0(t, param)))
    else:
        return -(np.exp(1j*(w1+w2)*t/2)/(w1+w2))*(np.cos((w1+w2)*t/2)*(fsin(t, w1, param)+fsin(t,w2,param))-np.sin((w1+w2)*t/2)*(fcos(t, w1, param)+fcos(t,w2, param)))


def JA(t, w1, w2, param):
    if(w1 == 0 and w2 == 0):
        return 0
    elif (w1 == -w2 and w2 != 0):
        return -1j*(ftsin(t,w1,param)-t*fsin(t,w1,param))
    elif (w1!=0 and w2==0):
        return -1j*(np.exp(1j*w1*t/2)/w1)*(-np.sin(w1*t/2)*fsin(t, w1, param)-np.cos(w1*t/2)*(fcos(t, w1, param)-f0(t, param)))
    elif (w1==0 and w2!=0):
        return -1j*(np.exp(1j*w2*t/2)/w2)*(np.sin(w2*t/2)*fsin(t, w2, param)+np.cos(w2*t/2)*(fcos(t, w2, param)-f0(t, param)))
    else:
        return -1j*(np.exp(1j*(w1+w2)*t/2)/(w1+w2))*(np.cos((w1+w2)*t/2)*(fcos(t, w1, param)-fcos(t, w2, param))+np.sin((w1+w2)*t/2)*(fsin(t, w1, param)-fsin(t, w2, param)))

def J(t, w1, w2, param):
    return JS(t,w1,w2,param) + JA(t,w1,w2,param)
