import numpy as np
from scipy.linalg import svd
from numpy import linalg
# define a matrix


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
    # correct signs
    # print(eig_vec_u)
    # print(sigma)
    # print(np.transpose(eig_vec_v))

    same_sign = np.sign((A @ eig_vec_v)[0] * (eig_vec_u @ np.diag(eig_val_u))[0])
    eig_vec_v = eig_vec_v * same_sign.reshape(1, -1)
  
    
    #Condition number - working
    condition_number = sigma[0]/sigma[len(sigma)-1]
    # print("++++++++++++++++++++")
    # print(condition_number)
    # print("++++++++++++++++++++")
    #Calculating Inverse - working 
       #Inverse check
    for ii in range(len(sigma)):
        if sigma[ii] == 0:
            print( "[ERROR] Inverse of matrix A does not exist.")
            return 0



    A_inv = np.matmul(np.matmul(eig_vec_v, np.diag(pow(sigma, -1))), np.transpose(eig_vec_u))
    print(A_inv)
    print(np.linalg.inv(A))
    # print("***********************")
    return 0

def SVD_Python(A):
    U, S, V = svd(A)
    print(U)
    print(S)
    print(V)





def main():
    A = np.array([[3,2,2],[2,3,-2],[4,5,6]])
    SVD(A)
    print("----------------------------------------------")
    SVD_Python(A)

if __name__=="__main__": 
    main() 
