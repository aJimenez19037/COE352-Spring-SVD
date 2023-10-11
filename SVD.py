import SVD
import math
import numpy as np
import sys
from scipy.linalg import svd
from numpy import linalg

g = 9.81
tolerance = .00001

def singular_value_decomposition(A):

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
    #Inverse check
    if np.any(sigma <= tolerance):
        print("[ERROR] Inverse does not exist. Singular value = 0")
        return 0,0,0,0,0



    same_sign = np.sign((A @ eig_vec_v)[0] * (eig_vec_u @ np.diag(eig_val_u))[0])
    eig_vec_v = eig_vec_v * same_sign.reshape(1, -1)
    
    #Condition number 
    condition_number = sigma[0]/sigma[len(sigma)-1]
    
    #Calculating Inverse 
  
    A_inv = np.matmul(np.matmul(eig_vec_v, np.diag(pow(sigma, -1))), np.transpose(eig_vec_u))
   
    return eig_vec_u, sigma, eig_vec_v, condition_number, A_inv

def main():

    A = np.array([[3,2,2],[2,3,-2],[4,5,6]])
    eig_vec_u, sigma, eig_vec_v, condition_number, A_inv = singular_value_decomposition(A)
    print("U = " + str(eig_vec_u))
    print("S = " + str(sigma))
    print("V = " + str(np.transpose(eig_vec_v)))
    print("A_INV = " + str(A_inv))

    print("---------------- BLACK BOX -------------------")
    U, S, V = svd(A)
    print("U = " + str(U))
    print("S = " + str(S))
    print("V = " + str(V))
    print("A_INV = " + str(np.linalg.inv(A)))

    print("----------------------------------------------")



if __name__ == "__main__":
    main()
