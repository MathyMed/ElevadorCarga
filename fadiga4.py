from dataclasses import dataclass
from re import A
from signal import raise_signal
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from fadiga1 import I

from fadiga2 import Mmax

l = 6
a = 8
OD = 2
ID = 1.5
r = 0.25
print("a = ", a, "in")
print("l = ", l, "in")
print("OD = ", OD, "in")
print("ID = ", ID, "in")
print("r = ", r, "in")

#Aluminio 2024-T4
Sut = 68000
Sy = 47000
E = 0
print("Sut = ", Sut, "psi")
print("Sy = ", Sy, "psi")
print("E = ", E)
print("\n")

N = 6e7

if Sut < 48000:
    Se_linha = 0.4*Sut
if Sut > 48000:
    Se_linha = 19000
print("Se_linha = ", Se_linha, "psi")
A95 = 0.010462*(OD**2)
print("A95 = ",A95, "in²")
d_equi = np.sqrt(A95/0.0766)
print("d_equi = ",d_equi, "in")
if OD < 0.3:
    C_tamanho = 1
if OD > 0.3 and OD < 10:
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

Sm = 0.9*Sut
N2 = [1e6, 5e6, 1e7, 5e7, 1e8, 5e8, 1e9, 5e9]
z_i = [-3, -3.699, -4, -4.699, -5, -5.699, -6, -6.699]
z_interp = interpolate.interp1d(N2, z_i, 'linear')
N_a = 5e8
z = z_interp(N_a)
print("z = ",z)
b = (1/z)*np.log10(Sm/Se)
print("b = ",b)
a_sn = 10**(np.log10(Sm) - 3*b)
print("a = ",a)
Sn = a_sn*(N)**b
print("Sn = ",Sn)

Sut_q = [15, 20, 30, 40, 50, 60, 70, 80, 90]
raiz_a = [0.475, 0.380, 0.278, 0.219, 0.186, 0.162, 0.144, 0.131, 0.122]
raiza_interp = interpolate.interp1d(Sut_q, raiz_a, 'linear')
x2 = Sut/1000
y2 = raiza_interp(x2)
print("raiz_a = ",y2)
q = 1/(1+(y2/np.sqrt(r)))
print("q = ",q)
Kt = 1.7
Kts = 1.35
Kf = 1 + q*(Kt - 1)
Kfs = 1 + q*(Kts - 1)
print("Kf = ", Kf)
print("Kfs = ", Kfs)

Fmax = 340
Fmin = -200
Fm = (Fmax + Fmin)/2
Fa = (Fmax - Fmin)/2
print("Fm = ", Fm, "lb")
print("Fa = ", Fa, "lb")
Ma = Fa*l
Mm = Fm*l
Ta = Fa*a
Tm = Fm*a
print("Ma = ", Ma, "lb*in")
print("Mm = ", Mm, "lb*in")
print("Ta = ", Ta, "lb*in")
print("Tm = ", Tm, "lb*in")
I = (np.pi*(OD**4 - ID**4))/64
J = (np.pi*(OD**4 - ID**4))/32
print("I = ", I, "in^4")
print("J = ", J, "in^4")
c = OD/2
sigma_max = Kf*((Fmax*l*c)/I)
tau_max = Kfs*((Fmax*a*c)/J)
print("sigma_max = ", sigma_max, "psi")
print("tau_max = ", tau_max, "psi")
sigma_vm = np.sqrt(sigma_max**2 + 3*(tau_max)**2)
print("sigma_vm = ", sigma_vm, "psi")
if sigma_vm < Sy:
    Kfm = Kf
    Kfsm = Kfs
if sigma_vm > Sy:
    Kfm = (Sy - Kf*sigma_max)/(np.abs(sigma_max))
    Kfsm = (Sy - Kf*tau_max)/(np.abs(tau_max))
print("Kfm = ", Kfm)
print("Kfsm = ", Kfsm)

sigma_a = Kf*(Ma*c/I)
sigma_ay = 0
tau_torcao_a = Kfs*(Ta*c/J)
sigma_m = Kfm*(Mm*c/I)
sigma_my = 0
tau_torcao_m = Kfsm*(Tm*c/J)
print("sigma_a = ", sigma_a, "psi")
print("tau_torcao_a = ", tau_torcao_a, "psi")
print("sigma_m = ", sigma_m, "psi")
print("tau_torcao_m = ", tau_torcao_m, "psi")
sigma_a_vm = np.sqrt(sigma_a**2 + sigma_ay**2 -sigma_a*sigma_ay + 3*(tau_torcao_a**2))
sigma_m_vm = np.sqrt(sigma_m**2 + sigma_my**2 -sigma_m*sigma_my + 3*(tau_torcao_m**2))
print("sigma_a_vm = ", sigma_a_vm, "psi")
print("sigma_m_vm = ", sigma_m_vm, "psi")

sigma_m_vm_Q = (1 - (sigma_a_vm/Sy))*Sy
Nf1 = sigma_m_vm_Q/sigma_m_vm
print("Nf1 = ", Nf1) 
sigma_a_vm_P = (1 - (sigma_m_vm/Sut))*Sn
Nf2 = sigma_a_vm_P/sigma_a_vm
print("Nf2 = ", Nf2)
Nf3 = (Sn*Sut)/(sigma_a_vm*Sut + sigma_m_vm*Sn)
print("Nf3 = ", Nf3)
sigma_m_vm_S = (Sut*(Sn**2 - Sn*sigma_a_vm + Sut*sigma_m_vm))/(Sn**2 + Sut**2)
sigma_a_vm_S = -(Sn/Sut)*sigma_m_vm_S + Sn
ZS = np.sqrt((sigma_m_vm - sigma_m_vm_S)**2 + (sigma_a_vm - sigma_a_vm_S)**2)
OZ = np.sqrt((sigma_a_vm)**2 + (sigma_m_vm)**2)
Nf4 = (OZ + ZS)/OZ
print("Nf4 = ", Nf4)

A = (np.pi/4)*(OD**2 - ID**2)
tau_trans_a = Kfs*(2*Fa/A)
print("tau_trans_a = ", tau_trans_a, "psi")
tau_trans_m = Kfsm*(2*Fm/A)
print("tau_trans_m = ", tau_trans_m, "psi")
tau_total_a = tau_torcao_a + tau_trans_a
tau_total_m = tau_torcao_m + tau_trans_m
print("tau_total_a = ", tau_total_a, "psi")
print("tau_total_m = ", tau_total_m, "psi")
sigma_a_vm_B = np.sqrt(3*(tau_total_a**2))
sigma_m_vm_B = np.sqrt(3*(tau_total_m**2))
print("sigma_a_vm_B = ", sigma_a_vm_B, "psi")
print("sigma_m_vm_B = ", sigma_m_vm_B, "psi")

sigma_m_vm_Q = (1 - (sigma_a_vm_B/Sy))*Sy
Nf1 = sigma_m_vm_Q/sigma_m_vm_B
print("Nf1 = ", Nf1) 
sigma_a_vm_P = (1 - (sigma_m_vm_B/Sut))*Sn
Nf2 = sigma_a_vm_P/sigma_a_vm_B
print("Nf2 = ", Nf2)
Nf3 = (Sn*Sut)/(sigma_a_vm_B*Sut + sigma_m_vm_B*Sn)
print("Nf3 = ", Nf3)
sigma_m_vm_S = (Sut*(Sn**2 - Sn*sigma_a_vm_B + Sut*sigma_m_vm_B))/(Sn**2 + Sut**2)
sigma_a_vm_S = -(Sn/Sut)*sigma_m_vm_S + Sn
ZS = np.sqrt((sigma_m_vm_B - sigma_m_vm_S)**2 + (sigma_a_vm_B - sigma_a_vm_S)**2)
OZ = np.sqrt((sigma_a_vm_B)**2 + (sigma_m_vm_B)**2)
Nf4 = (OZ + ZS)/OZ
print("Nf4 = ", Nf4)