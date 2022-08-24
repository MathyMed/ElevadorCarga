#Caso 3 - Carregamento estático de flexão
#Cisalhamento direto na solta
# + Cisalhamento devido à flexão no componente

import numpy as np

V = 10*1000*9.81
L = 0.2
M = V*L
e = 0.1
c = e/2
Iu = ((e)**3)/6

Sy = 331e6
fs = 2.5

h = np.sqrt(((((0.707*Sy)/fs)**2)**(-1))*((V/L)**2+(M*c/Iu)**2))
print("h = ", h,"m")
print("h = ", h*1000,"mm")

