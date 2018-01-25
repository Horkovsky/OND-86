import numpy
from math import pi, sqrt
import matplotlib.pyplot as plt

X = 300
Y = 100


def c_m(M, H, D, w): #Расчёт максимальной концентрации, Опасной скорости ветра и
    v = (pi*(D**2)/4)*w
    print('V: %s' % v)
    dT = 10
    A = 140  #Для европейской части.
    f = 1000 * ((w**2)*D)/((H**2)*dT)
    print('f: %s' % f)
    Um = 0.65 * ((v*5)/H)**(1/3)
    Um2 = 1.3 * (w * D)/H
    print('Um: %s' % Um)
    print('Um2: %s' % Um2)
    if f < 100:
        m = 1/(0.67 + 0.1*sqrt(f) + 0.34*(f**(1/3)))
    else:
        m = 1.47/(f**(1/3))
    print('m: %s' % m)
    if Um >= 2:
        n = 1
    elif Um >= 0.5:
        n = 0.532*(Um**2) - 2.13*Um + 3.13
    else:
        n = 4.4*Um
    print('n: %s' % n)

    Cm = (A*M*m*n)/(H**2*(v*(dT)**(1/3)))
    print('Cm: %s' % Cm)
    return Cm, Um2, f

def Pollution_ond (x0, y0, H, Cm, um, d = 6 ):
    Num = numpy.zeros((Y, X))
    F = 1
    Xm = ((5-F)/4)*d*H
    print('Xm: %s' % Xm)
    for j in range(x0+1, X):
        dx = j - x0
        if dx/Xm <= 1:
            S = 3*(dx/Xm)**4-8*(dx/Xm)**3+6*(dx/Xm)**2
        elif dx/Xm <=8:
            S = 1.113/(0.13*(dx/Xm)**2+1)
        else:
            S = 1/(0.1*(dx/Xm)**2+2.47*(dx/Xm)-17.8)
        Num[y0, j] = Cm * S

        for i in range(Y):
            dy = i - y0
            if um <= 5:
                ty = (um*(dy**2))/(dx**2)
            else:
                ty = (5*(dy**2))/(dx**2)
            S = ((1+5*ty+12.8*ty**2+17*ty**3+45.1*ty**4)**2)**(-1)
            if i != y0:
                Num[i, j] = Num[y0, j] * S
    return Num
tau = c_m(0.511, 5, 0.2, 10)
um2 = tau[1]
if um2 < 0.5:
    d = 5.7
elif um2 <=2:
    d = 11.4 * um2
else:
    d = 16 * sqrt(um2)

print('d: %s' % d)

o1 = Pollution_ond(0, 50, 1, tau[0], 150, d)
o2 = Pollution_ond(0, 60, 1, tau[0], 150, d)
o3 = o1+o2
plt.figure()
plt.imshow(o3)
plt.colorbar()
plt.show()