# COE352-Spring-SVD
This project contains both an singul value decomposition solver. It also contains a simple spring solver that uses the single value decomposition function to solve the spring system. Outputs condition number of K by calculating and printing the singular values (and eigenvalues) using the SVD algorithm, and calculates and prints a l2-condition number. It also returns the displacement values of the masses, elongation of the springs, internal forces, and the K inverse.

## Files
spring_solver.py - User inputted spring system solver
SVD.py -  Single value decomposition function. Outputs both the calculated SVD and the numpy blackbox solution

## Using the Spring Solver

Open your terminal and clone the repository using the command 
```
$ git clone https://github.com/aJimenez19037/COE352-Spring-SVD.git
$ cd COE352_Spring-SVD.git
$ python3 spring_solver.py (BC) (# of springs) [spring constants] [masses] 
```
The first argument BC is where you specify your boundary conditions. The options are: 

0 - Free free spring system  

1 - Fixed free spring system  

2 - Fixed fixed spring system  


Number of springs is where you specify the number of springs your system will have. Depending on your boundary condition, the number of masses will also be calculated. The logic is as follows:
```
if bc == 0:
    num_masses = num_springs + 1
elif bc == 1:
    num_masses = num_springs
elif bc == 2: 
    num_masses = num_springs - 1
```
The following two arguments you enter the spring constants and the masses. An example of a spring mass system with a boundary condition of 0 is shown below:
```
$ python3 spring_solver.py 2 4 1 1 1 1 1 1 1
Singular values = [3.41421356 2.         0.58578644]
Condition number = 5.828427124746185
l2 condition number = 4.0
++++++++++++
K_INV = [[0.75 0.5  0.25]
 [0.5  1.   0.5 ]
 [0.25 0.5  0.75]]
f = [9.81 9.81 9.81]
U = [14.715 19.62  14.715]
E = [ 14.715   4.905  -4.905 -14.715]
W = [ 14.715   4.905  -4.905 -14.715]
```
The spring solver outputs 5 variables. K_inv is the inverse of the K matrix used. f is a vector showing the forces acting on the system. U is the displacements of the masses. E is the elongation of the springs. And W is the internal forces of the system.
## Using the SVD.py file
Within the repository run:
```
$ python3 SVD.py
U = [[-0.5825206   0.18000222 -0.79263418]
 [-0.35945334 -0.93167999  0.05258985]
 [-0.72901512  0.31554967  0.60742519]]
S = [6.94389547 3.52224158 0.61329432]
V = [[-0.67015897 -0.63803475 -0.37920788]
 [-0.10695016 -0.42256798  0.89999887]
 [-0.73447166  0.64369866  0.2149498 ]]
A_INV = [[ 1.          0.         -0.66666667]
 [-0.8         0.2         0.66666667]
 [-0.2        -0.2         0.33333333]]
---------------- BLACK BOX -------------------
U = [[-0.5825206  -0.18000222  0.79263418]
 [-0.35945334  0.93167999 -0.05258985]
 [-0.72901512 -0.31554967 -0.60742519]]
S = [6.94389547 3.52224158 0.61329432]
V = [[-0.67015897 -0.63803475 -0.37920788]
 [ 0.10695016  0.42256798 -0.89999887]
 [ 0.73447166 -0.64369866 -0.2149498 ]]
A_INV = [[ 1.          0.         -0.66666667]
 [-0.8         0.2         0.66666667]
 [-0.2        -0.2         0.33333333]]
----------------------------------------------

```
If you wish to run your own custom matrix through the SVD.py file simply the A matrix within the main function. 

