import numpy as np
from scipy.integrate import odeint
from equation import lotVolRevised
import matplotlib.pyplot as plt


# initial condition
vibrio_0 = 0.5
ahl_0 = 0
ecoli_0 = 0.1
nisin_0 = 0

# parameters 
## Vibrio and AHL associated
Vibrio_GrothRate = 0.02
Vibrio_MaxCapacity = 1
Nisin_Vibrio_Binding_Strength = 0.05

Ahl_Secrete_Ratio = 0.02
Ahl_DecayRate = 0.01
Ahl_Secrete_Threshod = 0

## Ecoli and Nisin associated
Ecoli_GrothRate = 0.02
Ecoli_MaxCapacity = 1
Ahl_Ecoli_Binding_Strength = 0.05

Nsin_Sectete_Ratio = 0.03
Nisin_DecayRate = 0.01

## Mechenical associated
Dilution_Rate = 0


# time interval
t = np.arange(0,1000,0.1)

# load data
y0 = [vibrio_0, ahl_0, ecoli_0, nisin_0]
params = (
    Vibrio_GrothRate,
    Vibrio_MaxCapacity,
    Nisin_Vibrio_Binding_Strength,
    Ahl_Secrete_Ratio,
    Ahl_DecayRate,
    Ahl_Secrete_Threshod,
    Ecoli_GrothRate,
    Ecoli_MaxCapacity,
    Ahl_Ecoli_Binding_Strength,
    Nsin_Sectete_Ratio,
    Nisin_DecayRate,
    Dilution_Rate
)


sol = odeint(lotVolRevised,y0,t,args=params)
plt.plot(t, sol[:,0], 'r', label='vibrio')
plt.plot(t, sol[:,2], 'b', label='e-coli')
plt.legend(loc='best')
plt.grid()
plt.show()
