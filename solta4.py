#Caso 4 - Carregamento estático de torção
#Cisalhamento direto na solda
#Cisalhamnto devido ao giro do componente

import numpy as np

b = 50/1000
d = 100/1000
r = d/2
Sy = 331e6
fs = 2.5
F = 10*1000
L = 150/1000

Ju = ((8*(b**3) + 6*b*(d**2) + d**3)/12) - (b**4/(2*b + d))
print(Ju)
h = np.sqrt(((0.707*Sy/fs*F)**(-2))*((1/L**2) + (r/Ju)**2 + ((2*np.cos(53.13*np.pi/180)*r)/(L*Ju))))
print(h)