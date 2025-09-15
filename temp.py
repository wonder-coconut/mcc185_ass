import numpy as np
import matplotlib.pyplot as plt
from qutip import basis, sigmaz, sigmax, sesolve, Bloch

# Parameters
Omega = 2.0
detunings = [0.0, 1.0, 2.0, 5.0]
tlist = np.linspace(0, 10, 200)

# Initial state |0>
psi0 = basis(2, 0)

# Bloch sphere
b = Bloch()

# Custom colors for different detunings
colors = ['r', 'g', 'b', 'm']

for Delta, c in zip(detunings, colors):
    # Hamiltonian
    H = 0.5 * Delta * sigmaz() + 0.5 * Omega * sigmax()
    
    # Solve
    result = sesolve(H, psi0, tlist)
    
    # Add trajectory (list of states) as a line
    b.add_states(result.states)
    
    # Assign color for this trajectory
    b.line_color = [c]

# Show
b.show()
plt.show()