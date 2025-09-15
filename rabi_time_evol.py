import numpy as np
import matplotlib.pyplot as plt
import qutip
from qutip import Bloch, QobjEvo, basis, sesolve, sigmax, sigmaz
from qutip.measurement import measure,measurement_statistics

#hbar = 1.05457182e-34
hbar = 1

def HRWA(w0,w,omega0):
    return hbar * 0.5*((w-w0) * sigmaz()  + omega0 * sigmax())

def H_time_function(t,args):
    return np.cos(args['w'] * t)

def H_time_independent(w0,omega0):
    H0 = -hbar * 0.5 * w0 * sigmaz()
    H1 = hbar * omega0 * sigmax()
    return [H0,H1]

w0 = 20
omega0 = 1
detuning = [0]

timesteps = np.linspace(0,2*np.pi/omega0,1000)

psi0 = basis(2,1)
e_state = basis(2,1)
results = []

for delta in detuning:

    #H = HRWA(w0,w0+delta,omega0)
    H_temp = H_time_independent(w0,omega0)
    H0 = H_temp[0]
    H1 = H_temp[1]

    args = {'w':w0+delta}    

    H = QobjEvo([H0,[H1,H_time_function]],args = args,tlist = timesteps)


    res = sesolve(H, psi0, timesteps)

    results.append(res)

epops = [[],[],[],[]]

i = 0

b = [Bloch(),Bloch(),Bloch(),Bloch()]
colors = ['r','g','b','y']

while(i < len(detuning)):
    delta = detuning[i]
    res = results[i]
    op_states = res.states

    for state in op_states:
        epops[i].append(np.abs(state.overlap(e_state))**2)
    
    j = 0
    while(j < len(op_states)):
        b[i].add_states(op_states[j],'point',colors[i])
        j += 1

    i+=1


i = 0
while(i < len(detuning)):
    b[i].show()
    plt.show()
    i+=1