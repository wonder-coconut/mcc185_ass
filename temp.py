import qutip as qt

# Create Bloch sphere
b = qt.Bloch()

# Define a qubit state
psi = (qt.basis(2,0) + qt.basis(2,1)).unit()

# Convert the state to Bloch vector coordinates
vec = qt.bloch_vector(psi)  # returns [x, y, z]

# Add as a point (note: needs to be 2D: shape 3xN)
b.add_points([vec])  # wrap in list to make it 3x1

# Optional: customize
b.point_color = ['r']
b.point_marker = ['o']
b.point_size = [50]

# Show
b.show()
