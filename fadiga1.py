from dataclasses import dataclass
from signal import raise_signal
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

b = 2
d = 1
D = 1.125
r = 0.5
a = 5
l = 6
print("b = ", b, "in")
print("d = ", d, "in")
print("D = ", D, "in")
print("r = ", r, "in")
print("a = ", a, "in")
print("l = ", l, "in")

#Aço SAE 1040
Sut = 80000
E = 3e7
print("Sut = ", Sut, "psi")
print("E = ", E)
print("\n")

F = 500 #lb
print("F = ", F, "lb")
R = 500 #lb
print("R = ", R, "lb")
M = R*l - F*(l-a)
print("M = ",M," lb*in")
I = (b*d**3)/12
print("I = ",I," in^4")
c = d/2 
sigma_a_nom = M*c/I
print("sigma_a_nom = ",sigma_a_nom," psi")

D_d_k  = np.array([3, 2, 1.3, 1.2, 1.1, 1.05])
A_k = np.array([0.90720, 0.93232, 0.95880, 0.99590, 1.01650, 1.02260])
b_k = [-0.33333, -0.30304, -0.27269, -0.23829, -0.21548, -0.19156]
Ak_interp = interpolate.interp1d(D_d_k ,A_k, 'linear')
bk_interp = interpolate.interp1d(D_d_k ,b_k, 'linear')
x = D/d
ya = Ak_interp(x)
yb = bk_interp(x)
print("A = ",ya)
print("b = ",yb)
Kt = ya*(r/d)**yb
print("Kt = ",Kt)

Sut_q = [50, 55, 60, 70, 80, 90, 100, 110, 120, 130, 140, 160, 180, 200, 220, 240]
raiz_a = [0.130, 0.118, 0.108, 0.093, 0.080, 0.070, 0.062, 0.055, 0.049, 0.044, 0.039, 0.031, 0.024, 0.018, 0.013, 0.009]
raiza_interp = interpolate.interp1d(Sut_q, raiz_a, 'linear')
x2 = Sut/1000
y2 = raiza_interp(x2)
print("raiz_a = ",y2)
q = 1/(1+(y2/np.sqrt(r)))
print("q = ",q)

Kf = 1 + q*(Kt - 1)
print("Kf = ",Kf)
sigma_a = Kf*sigma_a_nom
print("sigma_a = ",sigma_a, "psi")
sigma_x = sigma_a
sigma_y = 0
tau_xy = 0
tau_ab = np.sqrt(((sigma_x+sigma_y)/2)**2 + tau_xy**2)
print("tau_ab = ",tau_ab, "psi")
sigma1_a = (sigma_x+sigma_y)/2 + tau_ab
sigma2_a = 0
sigma3_a = (sigma_x+sigma_y)/2 - tau_ab
sigma_vm = np.sqrt(sigma1_a**2 - sigma1_a*sigma2_a + sigma2_a**2)
print("sigma_vm = ",sigma_vm, "psi")

if Sut < 200000:
    Se_linha = 0.5*Sut
if Sut > 200000:
    Se_linha = 100000
print("Se_linha = ",Se_linha, "psi")
A95 = 0.05*d*b
print("A95 = ",A95, "in²")
d_equi = np.sqrt(A95/0.0766)
print("d_equi = ",d_equi, "in")
if d < 0.3:
    C_tamanho = 1
if d > 0.3 and d < 10:
    C_tamanho = 0.869*(d_equi)**(-0.097)
print("C_tamanho = ",C_tamanho)

C_carreg = 1  # 6.7a 1 = flexão  0.70 = força normal
print("C_carreg = ",C_carreg)
C_superf = 2.7*(Sut/1000)**(-0.265) # 6.7e Usinado a frio A = 2,7  b = –0,265
print("C_superf = ",C_superf)
T = 28  # 6.7f
if T < 450:
    C_temp = 1
if T > 450 and T < 550:
    C_temp = 1 - 0.0058*(T - 450)
print("C_temp = ",C_temp)
Conf_i = [50, 90, 95, 99, 99.9, 99.99, 99.999, 99.9999]
C_conf_i = [1, 0.897, 0.868, 0.814, 0.753, 0.702, 0.659, 0.620]
C_conf_interp = interpolate.interp1d(Conf_i, C_conf_i, 'linear')
Conf = 99.9
C_conf = C_conf_interp(Conf)
print("C_conf = ",C_conf)

Se = C_carreg*C_tamanho*C_superf*C_temp*C_conf*Se_linha
print("Se = ", Se, "psi")
Nf = Se/sigma_vm
print("Nf = ", Nf)

y = (F/(6*E*I))*(l**3 - 3*a*(l**2) - (l - a)**3)
print("Y = ", y, "in")
