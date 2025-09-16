import numpy as np
import matplotlib.pyplot as plt
import qutip
from qutip import Bloch, QobjEvo, basis, sesolve, sigmax, sigmaz
from qutip.measurement import measure,measurement_statistics

#hbar = 1.05457182e-34
hbar = 1
choose_hamiltonian = 1 #0 for RWA, 1 for non-RWA
choose_bloch = 1 #0 for population plot, 1 for bloch sphere

def HRWA(w0,w,omega0):
    return hbar * 0.5*((w-w0) * sigmaz()  + omega0 * sigmax())

def H_time_function(t,args):
    return np.cos(args['w'] * t)

def H_time_independent(w0,omega0):
    H0 = -hbar * 0.5 * w0 * sigmaz()
    H1 = hbar * omega0 * sigmax()
    return [H0,H1]

w0 = 20 #GHz
omega0 = 1 #GHz
detuning = [0,2]

timesteps = np.linspace(0,2*np.pi/omega0,1000) #nanosecond scale

psi0 = basis(2,1)
e_state = basis(2,1)
results = []
for delta in detuning:

    if(choose_hamiltonian):
        H_temp = H_time_independent(w0,omega0)
        H0 = H_temp[0]
        H1 = H_temp[1]
        args = {'w':w0+delta}    
        H = QobjEvo([H0,[H1,H_time_function]],args = args,tlist = timesteps)
    else:
        H = HRWA(w0,w0+delta,omega0)


    res = sesolve(H, psi0, timesteps)

    results.append(res)



epops = [[],[]]

i = 0

b = [Bloch(),Bloch()]
colors = ['r','g']

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

if(choose_bloch):
    while(i < len(detuning)):
        b[i].show()
        plt.savefig(f'../op_figs/bloch_{choose_hamiltonian}_{detuning[i]}.png')
        i += 1
else:
    while(i < len(detuning)):
        plt.plot(timesteps,epops[i], label=r'$\Delta$ = ' + f'{detuning[i]}' + r'$\cdot \Omega_0$')
        i+=1
    plt.xlabel('Time (ns)')
    plt.ylabel(r'Population of $|e\rangle$')
    plt.legend(loc='lower right')
    plt.title(r'$\omega_0$ = ' + f'{w0} GHz' + r', $\Omega_0$ = ' + f'{omega0} GHz')
    plt.savefig(f'../op_figs/population_{choose_hamiltonian}.png')