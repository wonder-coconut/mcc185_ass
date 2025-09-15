import numpy as np
import matplotlib.pyplot as plt
import qutip
from qutip import Bloch, QobjEvo, basis, sesolve, sigmax, sigmaz
from qutip.measurement import measure,measurement_statistics

#hbar = 1.05457182e-34
hbar = 1

def HRWA(w0,w,omega0):
    return 0.5*((w-w0)*sigmaz()  + omega0*sigmax())

w0 = 20
omega0 = 1
detuning = [0,1,2,5]

timesteps = np.linspace(0,10*np.pi/omega0,1000)

psi0 = basis(2,1)
e_state = basis(2,1)
results = []

for delta in detuning:

    H = HRWA(w0,w0+delta,omega0)

    res = sesolve(H, psi0, timesteps)
    results.append(res)

epops = [[],[],[],[]]

i = 0
while(i < len(detuning)):
    delta = detuning[i]
    res = results[i]
    op_states = res.states
    
    for state in op_states:
        epops[i].append(np.abs(state.overlap(e_state))**2)
    
    plt.plot(timesteps/1000,epops[i])

    i+=1


plt.show()