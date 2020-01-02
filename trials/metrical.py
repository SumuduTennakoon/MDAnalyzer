#CALCULATING GEOMETRIC RELATIONS FOR MINERALS USING THE METRICAL MATRIX"
#Sumudu Tennakoon (Created: 2017-10-15)
#Reference: G. V Gibbs, The Metrical Matrix in Teaching Mineralogy (1990).

from math import *
import numpy as np
 
class UnitCell:
    def __init__(self, a=1.0, b=1.0, c=1.0, alpha=90.0, beta=90.0, gamma=90.0):
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
    
    def v(self):
        return np.matrix([[a], [b], [c]])
    
    def v_T(self):
        return np.matrix([[a], [b], [c]]).T
    
    def G(self): 
        a = self.a
        b = self.b
        c = self.c
        alpha = radians(self.alpha)
        beta = radians(self.beta)
        gamma = radians(self.gamma)
        return np.matrix([[a * a, a * b * cos(gamma), a * c * cos(beta)], \
                         [a * b * cos(gamma), b * b, b * c * cos(alpha)], \
                         [a * c * cos(beta), b * c * cos(alpha), c * c]])
    def volume(self):
        a = self.a
        b = self.b
        c = self.c
        alpha = radians(self.alpha)
        beta = radians(self.beta)
        gamma = radians(self.gamma)
        return sqrt((a ** 2) * (b ** 2) * (c ** 2) * (1-(cos(alpha))**2-(cos(beta))**2-(cos(gamma))** 2 + 2*cos(alpha)*cos(beta)*cos(gamma)))

#Returns a position column vector.
def vector(x, y, z):
    return np.matrix([[x], [y], [z]])

#Angle between two vectors
def angle_theta(v1, v2, UC):
    G=UC.G() 
    cosTHETA = ((v1.T *G* v2).item(0)) / (sqrt((v1.T * G * v1).item(0)) * sqrt((v2.T * G * v2).item(0)))
    return degrees(acos(cosTHETA))

#Angle between two zones r1= [u1 v1 w1] and r2=[u2 v2 w2]
def angle_delta(r1, r2, UC):
    G=UC.G() 
    cosTHETA = ((r1.T *G* r2).item(0)) / (sqrt((r1.T * G * r1).item(0)) * sqrt((r2.T * G * r2).item(0)))
    return degrees(acos(cosTHETA))

#Angle between a face pole and and a zone
def angle_eps(s,r,UC): #s=face pole (h k l), r= zone [u v w] 
    G=UC.G()
    G_I=G.I
    cosEPS = ((s.T * r).item(0)) / (sqrt((s.T * G_I * s).item(0)) * sqrt((r.T * G * r).item(0)))
    return degrees(acos(cosEPS))

#Angle between two face poles
def angle_DELTA(s1,s2,UC):
    G_I=UC.G().I 
    cosDELTA = ((s1.T *G_I* s2).item(0)) / (sqrt((s1.T * G_I * s1).item(0)) * sqrt((s2.T * G_I * s2).item(0)))
    return degrees(acos(cosDELTA))

#Bond length (distance between two points/atoms)
def bond_length(r1, r2, UC):
    r = r2-r1
    G=UC.G() 
    return sqrt((r.T * G * r).item(0))

#d_hlk spacing
def d_hkl(s,UC):
    G_I=UC.G().I    
    Q_hkl=(s.T * G_I * s).item(0)
    return sqrt(1.0/Q_hkl)

