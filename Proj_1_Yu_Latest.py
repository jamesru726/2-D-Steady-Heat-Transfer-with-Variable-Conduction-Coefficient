import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import scipy as sp
import sys
import time

# Project 1 4933 CFD
# James Yu 9/23/24


################ Helper Functions ################ Helper Functions ################ Helper Functions ################ Helper Functions ################ Helper Functions 

# Building Outer Nodes
def BC_N(A, b, n_x, n_y, Dir_or_Neu, North_wall, South_wall, East_wall, West_wall, inv_delta_x, inv_delta_y, k1, k2):

      # Prepare index ranges
      north_indices = np.arange(n_x)
      south_indices = np.arange(n_x) + n_x * (n_y - 1)
      east_indices = (np.arange(2, n_y) * n_x - 1)
      west_indices = (np.arange(1, n_y - 1) * n_x)

      # Initialize north face using Neumann
      if Dir_or_Neu[0] == 'N':

            # Center node
            A[north_indices, north_indices] = -1     

            # Bottom node
            A[north_indices, north_indices + n_x] = 1

      # Initialize north face using Dirichelt
      elif Dir_or_Neu[0] == 'D':

            # North boundary condition
            b[north_indices] = North_wall
        
            # North boundary condition
            A[north_indices, north_indices] = 1


      # Initialize South face using Neumann
      if Dir_or_Neu[1] == 'N':

            # Center Node
            A[south_indices, south_indices] = -1

            # Top Node
            A[south_indices, south_indices - n_x] = 1

      # Initialize south face, excluding top left element, using Dirichelt
      elif Dir_or_Neu[1] == 'D':

            # South boundary condition
            b[south_indices] = South_wall

            # South boundary condition
            A[south_indices, south_indices] = 1

      # Initialize east face using Neumann
      if Dir_or_Neu[2] == 'N':
                  
            # Center Node
            A[east_indices, east_indices - 1] = 1

            # Left Node
            A[east_indices, east_indices] = -1

      # Initialize east face using Dirichelt
      elif Dir_or_Neu[2] == 'D':

            # East boundary condition
            b[east_indices] = East_wall

            # East boundary condition
            A[east_indices, east_indices] = 1

      # Initialize west face, excluding bottom left corner, using Neumann
      if Dir_or_Neu[3] == 'N':
                  
            # West boundary condition
            A[west_indices, west_indices] = -1
      
            # West boundary condition
            A[west_indices, west_indices + 1] = 1

      # Initialize west face, excluding bottom left corner, using Dirichelt
      elif Dir_or_Neu[3] == 'D':
            
            # West boundary condition
            b[west_indices] = West_wall
      
            # West boundary condition
            A[west_indices, west_indices] = 1

      # Building Internal Nodes
      for j in range(1, n_y - 1):

            # If statement used to check if we're in the top material
            if j < np.ceil(n_y/2) - 1:
                  k = k1
                  k_out = 0

            # If statement used to check if we're in the bottom material
            elif j > np.floor(n_y/2):
                  k = k2
                  k_out = 0

            # Use ceil and floor to deal with even and odd cases of grid length n_y
            elif j == np.ceil(n_y/2) - 1:
                  k = k1
                  k_out = (k2-k1)/4

            elif j == np.floor(n_y/2):
                  k = k2
                  k_out = (k2-k1)/4

            # Inner for loop used to iterate through the current row and assign coefficients
            for i in range(j * n_x + 1, (j + 1) * n_x - 1):

                  # Top node
                  A[i, i - n_x] = inv_delta_y * (k - k_out)

                  # Left and right node
                  A[i, i - 1] = A[i, i + 1] = inv_delta_x * k 
                  
                  # Center node
                  A[i, i] = -2*(inv_delta_x + inv_delta_y) * k

                  # Bottom node
                  A[i, i + n_x] = inv_delta_y * (k + k_out)

################ Helper Functions ################ Helper Functions ################ Helper Functions ################ Helper Functions ################ Helper Functions 





################ Adjustable parameters ################ Adjustable parameters ################ Adjustable parameters################ Adjustable parameters

# Grid size
n_x = 1000 # Grid Size in the x-direction
n_y = 1000 # Grid size in the y-direction

# Grid length
L_x = 1 # Length in the x-dir (meters)
L_y = 1 # Length in the y-dir (meters)

# Thermal conductivities (W/m*K)
k1 = 15  # Thermal conductivity of stainless steel
k2 = 0.1 # Thermal conductivity of some other material

# Choose the boundary condition for each face, Dirichlet or Neumann
# D for Dirichlet - Initializes outer nodes to specified value
# N for Neumann - Sets outer nodes equal to adjacent inner node
# [North, South, East, West]
Dir_or_Neu = ['D', 'D', 'N', 'D']

# Boundary conditions
# Temperature is converted from C to K for calculations
North_wall = 273 + 0 # north face
South_wall = 273 + 50 # south face
East_wall = 273 + 50 # east face
West_wall = 273 + 20 # west face

# Error handling: If the user enters anything beside D and N, exit code
if not all([x in ['N', 'D'] for x in Dir_or_Neu]):
      
      print('\nError: One of the strings in Dir_or_Neu are not D or N\nPlease enter D or N for Dir_or_Neu\n')
      sys.exit(1)

################ Adjustable parameters ################ Adjustable parameters ################ Adjustable parameters ################ Adjustable parameters





################ Initialization ################ Initialization ################ Initialization ################ Initialization ################ Initialization 
# Initialize values for calculations
# Solve for delta_x and delta_y
delta_x = L_x/(n_x - 1)
delta_y = L_y/(n_y - 1)

# Initialize inv_delta so they wouldn't be repeatedly calculated in the for loops
inv_delta_x = 1/delta_x**2
inv_delta_y = 1/delta_y**2

# Initialize sparse matrix for memmory and speed
# Used sparse type lil for efficient interative building of A coefficients and b values
b = sp.sparse.lil_matrix((n_x * n_y,1))
A = sp.sparse.lil_matrix((n_x * n_y, n_y * n_x))


################ Initialization ################ Initialization ################ Initialization ################ Initialization ################ Initialization 





################ Main ################ Main ################ Main ################ Main ################ Main ################ Main ################ Main 

num_tests = 50
Time2 = np.zeros(num_tests) 
avg_time1 = np.zeros((num_tests, 1))

for i in range(num_tests):
      start_time = time.time()
      # Building A matrix
      BC_N(A, b, n_x, n_y, Dir_or_Neu, North_wall, South_wall, East_wall, West_wall, inv_delta_x, inv_delta_y, k1, k2)
      end_time1 = time.time()
      avg_time1[i] = end_time1 - start_time

avg_time1 = np.average(avg_time1)
print(f'Average Time for building Matrix: {avg_time1}')
'''

# Solve for phi using SciPy spaarse spsolve
# Convert Matrix A to csr for efficient solving
start_time = time.time()
phi = sp.sparse.linalg.spsolve(A.tocsr(), b)

# Temperatures are converted back to C
phi = phi - 273

end_time2 = time.time()

Time2 = end_time2 - start_time

avg_time2 = np.average(Time2)
print(f'Average Time for building Matrix: {avg_time1}')
print(f'Average Time for solving Matrix: {avg_time2}')

################ Main ################ Main ################ Main ################ Main ################ Main ################ Main ################ Main 





################ Countor Plots ################ Countor Plots ################ Countor Plots ################ Countor Plots ################ Countor Plots 

# Setting up values for contour plot
# Starting at zero to represent the grid positions
x = np.linspace(0, L_x, n_x)
y = np.linspace(L_y, 0, n_y)
X, Y = np.meshgrid(x, y)
Z = np.round(phi.reshape(n_y, n_x), decimals = 2)

# Contour Plot
plt.contourf(x, y, Z, levels=np.max([n_x, n_y]), cmap='jet')
# Colobar label, tick, and temperature range adjustments
contour = plt.colorbar(label='Temperature (C)')
contour.formatter.set_scientific(False)  # Disable scientific notation

# Title and Axis labels
plt.title('Discretization of 2-D Heat Conduction with Varying Thermal Conductivity')
plt.xlabel('x, i')
plt.ylabel('y, j')
plt.show()

################ Countor Plots ################ Countor Plots ################ Countor Plots ################ Countor Plots ################ Countor Plots 


'''