import numpy
from math import pi, sqrt
import matplotlib.pyplot as plt
import matplotlib as mpl

X = 300                                 # Границы расчётной области
Y = 100


def c_m(M, H, D, w):  # Расчёт максимальной концентрации. Принимает на вход мощность выброса (М), высоту источника (H),
                        # диаметр устья(D), скорость выхода ГВС (w)

    v = ((pi*D**2)/4)*w                 # расчёт объёмного расхода ГВС
    print('V: %s' % v)
    F = 1                               # Коэффициент оседания 1 для газа. 2 для мелкодисперсных аэрозолей
    dT = 10                            # Разница температур ГВС и окружающей среды
    A = 140                             # Коэффициент стратификации. Принят для Москвы и области.
    f = 1000 * ((w**2)*D)/((H**2)*dT)   # вспомогательные
    vm = 0.65 * ((v*dT)/H)**(1/3)       #    параметры, характеризующие источник
    vm2 = 1.3 * (w * D)/H               #                для рассчёта опасной скорости ветра
    fe = 800 * (vm2**3)
    print('f: %s' % f)                  #                       и коэффициентов m и n в конечной формуле
    print('fe: %s' % fe)
    print('vm: %s' % vm)
    print('vm2: %s' % vm2)
    if f < 100:
        m = 1/(0.67 + 0.1*sqrt(f) + 0.34*(f**(1/3)))
    else:
        m = 1.47/(f**(1/3))
    print('m: %s' % m)
    if vm >= 2:
        n = 1
    elif vm >= 0.5:
        n = 0.532*(vm**2) - 2.13*vm + 3.13
    else:
        n = 4.4*vm
    print('n: %s' % n)
    if f >= 100 and vm2 >= 0.5:
        K = D/8*v
        Cm = ((A*M*n*F)/(H**(4/3)))*K
    elif f < 100 and vm2 < 0.5:
        m = m * 2.86
        Cm = (A * M * m * F) / (H**(7/3))
    elif f >= 100 and vm2 < 0.5:
        m = 0.9
        Cm = (A * M * m * F) / (H ** (7 / 3))
    else:
        Cm = (A*M*m*n*F)/(H**2 *((v*dT)**(1/3)))

    if f >= 100:
        if vm2 < 0.5:
            um = 0.5
            d = 5.7
        elif vm2 <= 2:
            d = 11.4 * vm2
            um = vm
        else:
            d = 16 * sqrt(vm2)
            um = 2.2 * vm2
    else:
        if vm2 < 0.5:
            d = 2.48 * (1 + 0.28*(fe**(1/3)))
            um = 0.5
        elif vm2 <= 2:
            d = 4.95 * vm * (1 + 0.28 * (f**(1/3)))
            um = vm
        else:
            d = 7 * sqrt(vm) * (1 + 0.28 * (f ** (1 / 3)))
            um = vm * (1 + 0.12*sqrt(f))
    print('Cm: %s' % Cm)
    print('d: %s' % d)
    print('um: %s' % um)
    return Cm, um, d


def Pollution_ond (x0, y0, H, Cm, um, d = 6 ):  # Собственно рассчёт рассеивания. Вычисление концентрации по оси факела
                                                # с последующим расчётом в точках, перпендиклярным оси
    Num = numpy.zeros((Y, X))
    F = 1
    Xm = ((5-F)/4)*d*H
    print('Xm: %s' % Xm)
    for j in range(x0+1, X):
        dx = j - x0
        if dx/Xm <= 1:
            if H >= 2 and H < 11:
                S = 0.125 * (10 - H) + 0.125 * (H - 2)
            else:
                S = 3*(dx/Xm)**4 - 8*(dx/Xm)**3 + 6*(dx/Xm)**2
        elif dx/Xm <= 8:
            S = 1.113/(0.13*(dx/Xm)**2 + 1)
        elif dx/Xm > 8 and dx/Xm <= 100:
            if F <= 1.5:
                S = (dx/Xm)/(3.556*(dx/Xm)**2 - 35.2*(dx/Xm) + 120)
            else:
                S = 1 / (0.1 * (dx / Xm) ** 2 + 2.456 * (dx / Xm) - 17.8)
        else:
            if F <= 1.5:
                S = 144.3 *(dx/Xm)**(-7/3)
            else:
                S = 37.76 *(dx/Xm)**(-7/3)
        Num[y0, j] = (Cm * S)/0.2

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

o1 = Pollution_ond(0, 50, 5, tau[0], tau[1], tau[2])

# o2 = Pollution_ond(0, 60, 1, tau[0], tau[1], d)
# o3 = o1+o2
plt.figure()
cmap = mpl.colors.ListedColormap(['blue', 'green', 'xkcd:pumpkin', 'red', 'black'])
bounds = [0.1, 0.5, 1, tau[0]/0.2/2, tau[0]/0.2-0.25, tau[0]/0.2]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
plt.imshow(o1, cmap=cmap, norm=norm)
plt.show()