This code uses the discretization of the 2-D heat conduction equation with variable thermal conductivity using central difference scheme.
This is under the assumption that this block has two different materials at the top and bottom and is split horizontally in the middle.

User can change the temperature, grid size (number of nodes), grid size (actual size in meters), and thermal conductivity of each material.

Flexibility in boundary conditions
  - can choose between Dirichlet and Neumann for each face
  - can choose temperatures of each face

This code utlizes sparse matrices for speed and memory.
For grid size 1000x1000
Average Time for building Matrix: 5.71251127243042
Average Time for solving Matrix: 18.544085969924925
Note: This varies from PC to PC

