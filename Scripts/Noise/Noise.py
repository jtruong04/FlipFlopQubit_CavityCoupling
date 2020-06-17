import numpy as np

class Noise:
    def Sw(self, w):
        return 0

    def St(self, t):
        return 0

    def f0(self, t):
        return 0

    def ft(self, t):
        return 0

    def fsin(self, t, w):
        return 0

    def fcos(self, t, w):
        return 0

    def ftsin(self, t, w):
        return 0

    def ftcos(self, t, w):
        return 0

    def JS(self, t, w1, w2):
        if(w1 == 0 and w2 == 0):
            return -(self.ft(t)-t*self.f0(t))
        elif (w1 == -w2 and w2 != 0):
            return -(self.ftcos(t, w1)-t*self.fcos(t, w1))
        elif (w1 != 0 and w2 == 0):
            return -(np.exp(1j*w1*t/2)/w1)*(np.cos(w1*t/2)*self.fsin(t, w1)-np.sin(w1*t/2)*(self.fcos(t, w1)+self.f0(t)))
        elif (w1 == 0 and w2 != 0):
            return -(np.exp(1j*w2*t/2)/w2)*(np.cos(w2*t/2)*self.fsin(t, w2)-np.sin(w2*t/2)*(self.fcos(t, w2)+self.f0(t)))
        else:
            return -(np.exp(1j*(w1+w2)*t/2)/(w1+w2))*(np.cos((w1+w2)*t/2)*(self.fsin(t, w1)+self.fsin(t, w2))-np.sin((w1+w2)*t/2)*(self.fcos(t, w1)+self.fcos(t, w2)))

    def JA(self, t, w1, w2):
        if(w1 == 0 and w2 == 0):
            return 0
        elif (w1 == -w2 and w2 != 0):
            return -1j*(self.ftsin(t, w1)-t*self.fsin(t, w1))
        elif (w1 != 0 and w2 == 0):
            return -1j*(np.exp(1j*w1*t/2)/w1)*(-np.sin(w1*t/2)*self.fsin(t, w1)-np.cos(w1*t/2)*(self.fcos(t, w1)-self.f0(t)))
        elif (w1 == 0 and w2 != 0):
            return -1j*(np.exp(1j*w2*t/2)/w2)*(np.sin(w2*t/2)*self.fsin(t, w2)+np.cos(w2*t/2)*(self.fcos(t, w2)-self.f0(t)))
        else:
            return -1j*(np.exp(1j*(w1+w2)*t/2)/(w1+w2))*(np.cos((w1+w2)*t/2)*(self.fcos(t, w1)-self.fcos(t, w2))+np.sin((w1+w2)*t/2)*(self.fsin(t, w1)-self.fsin(t, w2)))

    def J(self, t, w1, w2):
        if(t == 0):
            return 0
        return self.JS(t, w1, w2) + self.JA(t, w1, w2)

    def J(self, t):
        if(t == 0):
            return 0
        return self.JS(t, 0,0) + self.JA(t, 0,0)
