# COE352-Spring-SVD

## Using the Spring Solver

Open your terminal and clone the repository using the command 
```
$ git clone https://github.com/aJimenez19037/COE352-Spring-SVD.git
$ cd COE352_Spring-SVD.git
$ python3 main.py (BC) (# of springs) [spring constants] [masses] 
```
The first argument BC is where you specify you boundary conditions. The options are: 
0 - Free free spring system
1 - Fixed free spring system
2 - Fixed fixed spring system

"#" of springs is where you specify the number of springs your system will have. Depending on your boundary condition, the number of masses will also be calculated. The logic is as follows:
```
if bc == 0:
    num_masses = num_springs + 1
elif bc == 1:
    num_masses = num_springs
elif bc == 2: 
    num_masses = num_springs - 1
```
