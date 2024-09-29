This code uses the discretization of the 2-D heat conduction equation with variable heat transfer coeffcient using central difference scheme.
This is under the assumption that this block has two different materials at the top and bottom and is split perfectly in the middle.

User can change the temperature, grid size (number of nodes), grid size (actual size in meters), and thermal conductivity of each material.

Flexibility in boundary conditions
  - can choose between Dirichlett and Neumann
  - can choose temperatures of each face

Can handle different grid sizes that are not square.
