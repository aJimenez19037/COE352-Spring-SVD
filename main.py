import math
import numpy as np
from scipy.linalg import svd
from numpy import linalg

import sys

g = 9.81

def SVD(A):

    A_T = np.transpose(A)

    eig_val_u, eig_vec_u = np.linalg.eig(np.matmul(A,A_T))
    _, n = eig_vec_u.shape
    #sort U
    sorted_indices = np.argsort(eig_val_u)[::-1]
    eig_val_u = eig_val_u[sorted_indices]
    eig_vec_u = eig_vec_u[:, sorted_indices]
    #sort V_T
    eig_val_v, eig_vec_v = np.linalg.eig(np.matmul(A_T,A))
    _, m = eig_vec_v.shape
    sorted_indices = np.argsort(eig_val_v)[::-1]
    eig_val_v = eig_val_v[sorted_indices]
    eig_vec_v = eig_vec_v[:, sorted_indices]
    #calc sigma
    sigma = np.sqrt(eig_val_v[:min(m,n)])

    same_sign = np.sign((A @ eig_vec_v)[0] * (eig_vec_u @ np.diag(eig_val_u))[0])
    eig_vec_v = eig_vec_v * same_sign.reshape(1, -1)
    
    #Condition number 
    condition_number = sigma[0]/sigma[len(sigma)-1]
    
    #Calculating Inverse 
       #Inverse check
    for ii in range(len(sigma)):
        if sigma[ii] == 0:
            print( "[ERROR] Inverse of matrix A does not exist.")
            return 0

    A_inv = np.matmul(np.matmul(eig_vec_v, np.diag(pow(sigma, -1))), np.transpose(eig_vec_u))
   
    return eig_vec_u, sigma, eig_vec_v, condition_number, A_inv


def spring_solver(bc, num_springs, num_masses, spring_const, masses):

    # Elongation equations
    # u =  displacements
    print(num_masses)
    if bc == 1: #1 fixed ends
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
        eig_vec_u, sigma, eig_vec_v, condition_number, K_INV=SVD(K)
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
        f = f * g
        
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
        print("\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\")
        U = [0]*num_masses
        E = [0]*num_masses

        #Create A matrix
        A = np.zeros((num_springs, num_masses))
    
        A[0, 0] = 1
        A[0, 1] = 0
        
        # Set other e values = u2 - u1 except shifted for their respective u values
        for i in range(1, num_springs - 1):
            A[i, i - 1] = -1
            A[i, i] = 1
        
        # Set e4 = -u3
        A[num_springs - 1, num_springs - 2] = -1    

        print("A = " + str(A))

        # Spring Const Matrix
        C = np.diag(spring_const)

        # SVD of K matrix
        K = np.matmul(np.matmul(np.transpose(A), C),A)
        eig_vec_u, sigma, eig_vec_v, condition_number, K_INV=SVD(K)
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
        f = f * g
        
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
    

def main():
    A = np.array([[3,2,2],[2,3,-2],[4,5,6]])
    eig_vec_u, sigma, eig_vec_v, condition_number, A_inv = SVD(A)
    print("U = " + str(eig_vec_u))
    print("S = " + str(sigma))
    print("V = " + str(eig_vec_v))
    print("A_INV = " + str(A_inv))

    print("---------------- BLACK BOX -------------------")
    U, S, V = svd(A)
    print("U = " + str(U))
    print("S = " + str(S))
    print("V = " + str(V))
    print("A_INV = " + str(np.linalg.inv(A)))

    print("----------------------------------------------")

    arguments = sys.argv
    bc = int(arguments[1])
    num_springs = int(arguments[2])
    if bc == 1:
        num_masses = num_springs
    elif bc == 2: 
        num_masses = num_springs - 1
    else:
        print("[ERROR] Please input a boundary condition of 1 or 2.")
   
    
    spring_const = []
    masses = []
    for ii in range(3,num_springs+3):
        spring_const.append(int(arguments[ii]))
    for ii in range(3+num_springs, 3+num_springs+num_masses):
        masses.append(int(arguments[ii]))
    
    spring_solver(bc,num_springs,num_masses,spring_const,masses)

    return 0




if __name__=="__main__": 
    main() 
