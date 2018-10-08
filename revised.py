import numpy as np
from scipy.integrate import solve_ivp
from equation import lotVolRevised
import matplotlib.pyplot as plt

# initial condition
vibrio_0 = 0.2
ahl_0 = 0
ecoli_0 = 0.01
nisin_0 = 0
suicide_0 = 0

params = {
   ## Vibrio and AHL associated:
    # Vibrio_GrothRate
    'rV': 0.02,
    # Vibrio_MaxCapacity
    'Vm': 1,
    # Nisin_Vibrio_Binding_Strength
    'a': 0.05,
    # Ahl_Secrete_Ratio
    'kA': 0.02,
    # Ahl_DecayRate
    'lmA': 0.01,
    
   ## Ecoli, Nisin, Suicide associated:
    # Ahl_Ecoli_Binding_Strength
    'b': 0.01,
    # Suicide_Ecoli_Binding_Strength
    'c': 0.04,
    # Ecoli_DeathRate
    'lmE': 0.01,
    # Nsin_Sectete_Ratio
    'kN': 0.02,
    # Nisin_DecayRate
    'lmN': 0.01,
    # Suicide_Secrete_Ratio
    'kS': 0.04,
    # Suicide_DeacyRate
    'lmS': 0.01,
    # Vibrio_Threshold
    'thesV': 0.4,
    
   ## Mechnical:
    # Dilution_Rate
    'D': 0
}


# load data
t_span = [0,10000] 
y_0 = [vibrio_0, ahl_0, ecoli_0, nisin_0, suicide_0]


sol = solve_ivp(
    lambda t_span, y_0: lotVolRevised(t_span, y_0,**params),t_span, y_0)
    # ,t_eval=np.arange(0,10000,0.01) )


# create a figure and add two subplot
fig, (ax1,ax2) = plt.subplots(1,2)

# title of subplot 1
ax1.set_title('Revolution of time')
ax1.plot(sol.t, sol.y[0], 'tomato', label='vibrio')
ax1.plot(sol.t, sol.y[2], 'dodgerblue', label='e-coli')
ax1.set(ylabel='concentration (M)')
ax1.set(xlabel='time ($min^{-1}$)')
# open the legend box and auto adjust location
ax1.legend(loc='best')
ax1.grid()

# title of subplot 2
ax2.set_title('Phase diagram')
ax2.plot(sol.y[0], sol.y[2],'mediumseagreen')
ax2.set(ylabel='Vibrio (M)')
ax2.set(xlabel='Ecoli (M)')
ax1.grid()
# prevention of firgure overlap
plt.tight_layout()
plt.savefig('result.png')