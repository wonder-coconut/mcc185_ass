import numpy as np
import matplotlib.pyplot as plt
from qutip import basis, sigmaz, sigmax, sesolve

# Parameters
Omega = 2.0                   # Rabi frequency
detunings = [0.0, 1.0, 2.0, 5.0]  # list of detunings to try
tlist = np.linspace(0, 10, 200)

# Initial state |0>
psi0 = basis(2, 0)

# Container for results
populations = {}

for Delta in detunings:
    # Hamiltonian
    H = 0.5 * Delta * sigmaz() + 0.5 * Omega * sigmax()
    
    # Solve Schrodinger equation
    result = sesolve(H, psi0, tlist)
    
    # Population in |1>
    p_excited = [abs(state.overlap(basis(2,1)))**2 for state in result.states]
    populations[Delta] = p_excited

plt.figure(figsize=(8,6))

for Delta, pops in populations.items():
    plt.plot(tlist, pops, label=fr'Delta = {Delta}')

plt.xlabel("Time")
plt.ylabel("Excited state population P(|1>)")
plt.title("Rabi oscillations vs detuning")
plt.legend()
plt.show()
