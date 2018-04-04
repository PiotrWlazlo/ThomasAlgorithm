##################################################################
#                                                                #
#   Python 3.6 program which solve tridiagonal sparse 128x128    #
#   matrix using Thomas algorithm. This algorithm is simplified  #
#   form of Gauss elimination. His time complexity is O(n).      #
#                                                                #
#   Piotr Wlazlo                                    04.12.2017   #
##################################################################


# -*- coding: cp1250 -*-
import numpy as np
import copy

#Create matrices
A = np.diag([2.0,1.0,1.0,2.0,1.0,2.0,1.0,2.0],-1)\
    +np.diag([4.0,8.0,4.0,3.0,4.0,5.0,3.0,8.0,5.0])\
    +np.diag([2.0,1.0,1.0,2.0,1.0,2.0,1.0,2.0],1)

#b column
b = np.array([1.,2.,3.,4.,5.,6.,7.,8.,9.])


#Thomas algorithm - LU decomposition
L=np.diag([1.0]*(9),0)
U=np.zeros(shape=[9,9])
U[0][0]=A[0][0] #because f1 = 1*f1'
#LU decomposition
for i in range(0,8):
    U[i][i+1]=A[i][i+1]
    L[i+1][i] = (A[i+1][i]/U[i][i])
    U[i+1][i+1]=A[i+1][i+1]-L[i+1][i]*U[i][i+1]

#Tworzę kopię wektora do obliczeń
#Creating copy of vector for calculations
b_temp=copy.copy(b)
#Creating vector x
x = [0.,0.,0.,0.,0.,0.,0.,0.,0.]
#Making forward-substitution
for i in range(1,9):
    b_temp[i] = b[i]-((L[i][i-1])*(b_temp[i-1]))
#Making back-substitution
x[8]=b_temp[8]/U[8][8]
for i in range(7,-1,-1):
    x[i]=(b_temp[i]-((A[i][i+1])*(x[i+1])))/U[i][i]
for i in range(9):
    print("x%d"%(i+1)+" = %f"%x[i])
    

