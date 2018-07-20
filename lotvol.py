import numpy as np
from scipy.integrate import odeint
from equation import lotVol
import matplotlib.pyplot as plt


# initial condition
vibrio_0 = 1
ahl_0 = 0
ecoli_0 = 1
nisin_0 = 0

# parameters 
Vibrio_GrothRate = 0.1
Ahl_Vibrio_Ratio = 1
Ahl_DecayRate = 0
Ahl_Secrete_Threshod = 1

Ecoli_GrothRate = 0.1
Nsin_Ecoli_Ratio = 1
Nisin_DecayRate = 0

Nisin_Vibrio_Binding_Strength = 0.05
Ahl_Ecoli_Binding_Strength = 0.5

# time interval
t = np.arange(0,100,0.1)

# load data
y0 = [vibrio_0, ahl_0, ecoli_0, nisin_0]
params = (
    Vibrio_GrothRate, 
    Nisin_Vibrio_Binding_Strength, 
    Ecoli_GrothRate, 
    Ahl_Ecoli_Binding_Strength,
    Ahl_DecayRate, 
    Nisin_DecayRate,
    Ahl_Vibrio_Ratio,
    Nsin_Ecoli_Ratio,
    Ahl_Secrete_Threshod)


sol = odeint(lotVol,y0,t,args=params)

# create a figure and add a subplot(111)
fig ,ax = plt.subplots()

ax.set_title('Revised model')
ax.plot(t, sol[:,0], 'r', label='vibrio')
ax.plot(t, sol[:,2], 'b', label='e-coli')
ax.set(xlabel='time ($hr^{-1}$)')
ax.set(ylabel='concentration (M)')
# open the legend box
ax.legend(loc='best')
ax.grid()

plt.show()
