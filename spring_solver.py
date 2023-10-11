import SVD
import math
import numpy as np
import sys
import warnings
from scipy.linalg import svd
from numpy import linalg


def spring_solver(bc, num_springs, num_masses, spring_const, masses):

    # Elongation equations
    # u =  displacements
    U = [1]*num_masses
    if bc == 0: # free free system
        A = np.zeros((num_springs,num_masses))
        for i in range(num_springs):
            # Set the diagonal elements (e1 = u1 - u0, e2 = u2 - u1, e3 = u3 - u2)
            A[i, i] = 1
            
            # Set the off-diagonal elements (e1 = u1 - u0, e2 = u2 - u1, e3 = u3 - u2)
            A[i, i + 1] = -1
        print("A = " + str(A))
        check = np.matmul(A,U)
        for ii in range(num_springs):
            if check[ii]!=0:
                print("Matrix does not have a solution")
                return 0

        # Spring Const Matrix
        C = np.diag(spring_const)

        # SVD of K matrix
        K = np.matmul(np.matmul(np.transpose(A), C),A)
        print("K = " + str(K))
        eig_vec_u, sigma, eig_vec_v, condition_number, K_INV = SVD.singular_value_decomposition(K)

        if (type(K_INV) == int):
            return 0
        print("Singular values = " + str(sigma))
        print("Condition number = " + str(condition_number))

        # l2 norm
        sum = 0
        for ii in range (K.shape[0]):
            for jj in range (K.shape[1]):
                sum = sum + pow(K[ii][jj],2)
        l2 = math.sqrt(sum)
        print("l2 condition number = " + str(l2))

        #f vector 
        f = np.array(masses)
        f = f * SVD.g
        
        # k_inv * f = u 
        U = np.matmul(K_INV,f)
        print("++++++++++++")
        print(K_INV)
        print(f)
        print("U = " + str(U))

        # calculate Elongation
        E = np.matmul(A, U)
        print("E = " + str(E))

        # calculate W aka internal force 
        print("C = " + str(C))
        W = np.matmul(C, E)
        print("W = " + str(W))
        


        

        # A check



    elif bc == 1: #1 fixed ends
        U = [0]*num_masses
        E = [0]*num_masses

        #Create A matrix
        A = np.zeros((num_masses, num_masses))
        for i in range(num_masses):
            # Set the diagonal elements (e1 = u1 - u0, e2 = u2 - u1, e3 = u3 - u2)
            A[i, i] = 1
            
            # Set the off-diagonal elements (e1 = u1 - u0, e2 = u2 - u1, e3 = u3 - u2)
            if i > 0:
                A[i, i - 1] = -1

        print("A = " + str(A))

        # Spring Const Matrix
        C = np.diag(spring_const)

        # SVD of K matrix
        K = np.matmul(np.matmul(np.transpose(A), C),A)
        eig_vec_u, sigma, eig_vec_v, condition_number, K_INV=SVD.singular_value_decomposition(K)
        print("Singular values = " + str(sigma))
        print("Condition number = " + str(condition_number))

        # l2 norm
        sum = 0
        for ii in range (K.shape[0]):
            for jj in range (K.shape[1]):
                sum = sum + pow(K[ii][jj],2)
        l2 = math.sqrt(sum)
        print("l2 condition number = " + str(l2))

        #f vector 
        f = np.array(masses)
        f = f * SVD.g
        
        # k_inv * f = u 
        U = np.matmul(K_INV,f)
        print("++++++++++++")
        print(K_INV)
        print(f)
        print("U = " + str(U))

        # calculate Elongation
        E = np.matmul(A, U)
        print("E = " + str(E))

        # calculate W aka internal force 
        print("C = " + str(C))
        W = np.matmul(C, E)
        print("W = " + str(W))
        

        return 0
    elif bc == 2: # 2 fixed ends
        U = [0]*num_masses
        E = [0]*num_masses

        #Create A matrix
        A = np.zeros((num_springs, num_masses))
    
        if num_masses == 1:
            A = np.array([[1], [-1]])
        else:
            A[0, 0] = 1
            A[0, 1] = 0
            # Set other e values = u2 - u1 except shifted for their respective u values
            for i in range(1, num_springs - 1):
                A[i, i - 1] = -1
                A[i, i] = 1
                # Set e4 = -u3
                A[num_springs - 1, num_springs - 2] = -1    

        # Spring Const Matrix
        C = np.diag(spring_const)

        # SVD of K matrix
        K = np.matmul(np.matmul(np.transpose(A), C),A)
        eig_vec_u, sigma, eig_vec_v, condition_number, K_INV=SVD.singular_value_decomposition(K)
        print("Singular values = " + str(sigma))
        print("Condition number = " + str(condition_number))

        # l2 norm
        sum = 0
        for ii in range (K.shape[0]):
            for jj in range (K.shape[1]):
                sum = sum + pow(K[ii][jj],2)
        l2 = math.sqrt(sum)
        print("l2 condition number = " + str(l2))

        #f vector 
        f = np.array(masses)
        f = f * SVD.g
        
        # k_inv * f = u 
        U = np.matmul(K_INV,f)
        print("++++++++++++")
        print("K_INV = " + str(K_INV))
        print("f = " + str(f))
        print("U = " + str(U))

        # calculate Elongation
        E = np.matmul(A, U)
        print("E = " + str(E))

        # calculate W aka internal force 
        W = np.matmul(C, E)
        print("W = " + str(W))
        


        return 0
    

def main():

    arguments = sys.argv
    bc = int(arguments[1])
    num_springs = int(arguments[2])
    if num_springs == 0:
        print("[ERROR]: Cannot enter 0 springs")
        return 0
    if bc == 0:
        num_masses = num_springs + 1
    elif bc == 1:
        num_masses = num_springs
    elif bc == 2: 
        num_masses = num_springs - 1
        if num_springs == 1:
            print("[ERROR]: BC of 2 requires at least 2 springs")
            return 0
    else:
        print("[ERROR] Please input a boundary condition of 1 or 2.")
   
    
    spring_const = []
    masses = []
    for ii in range(3,num_springs+3):
        spring_const.append(int(arguments[ii]))
    for ii in range(3+num_springs, 3+num_springs+num_masses):
        masses.append(int(arguments[ii]))

    if (3+num_springs+num_masses != len(arguments)):
        print("[ERROR] Number of arguments does not match given the BC and number of springs")

    
    spring_solver(bc,num_springs,num_masses,spring_const,masses)

    return 0




if __name__=="__main__": 
    main() 
