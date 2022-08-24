#Caso 1 - Carregamento estático paralelo ao cordão da solda
#Apenas cisalhemento direto na solda

F = 10*1000
P = F*9.81
fs = 2.5
Sy = 331*10**6
L = 0.2 

h = (1.41*P*fs)/(Sy*L)
print("h = ", h)

F = 19.2*1000
P = F*9.81
fs = 2.5
Sy = 331*10**6
L = 0.30

h = (1.41*P*fs)/(Sy*L)
print("h = ", h)