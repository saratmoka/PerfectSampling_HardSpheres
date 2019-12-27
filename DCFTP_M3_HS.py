# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 12:55:29 2018

@author: uqsbabum
"""

# Program to generate Strauss process on [0,1]^2 with parameter Gamma and abosolutely continuous with respect to 
# the Poisson point process with intensity Kappa = Kappa_0/(pi*IntRange*IntRange), where IntRange is 
# the interaction range of the model.


def reset():
    try:
        get_ipython().magic('reset -sf')  #analysis:ignore
    except NameError:
        pass
reset()


import numpy as np
#import time
import datetime
#import math
import matplotlib.pyplot as plt




###### Euclidean distance between two points on R^2 
def euclidean2dim(yy, xx):
    return np.sqrt((xx[0] - yy[0])**2 + (xx[1] - yy[1])**2)

###### For checking overlap under the assumption the underlying space is a torus
def overlap(C, NC, IntRange, torus_flag = True): # returns True if there is a overlap. Otherwise, False
    if IntRange > 1:
        print('Error: Interaction range is bigger than 1')
        exit(0)
        
    NC = NC[0:2]
    if torus_flag == False: 
        for x in C:
            if euclidean2dim(x, NC) < IntRange:
                return True
            
    elif  (IntRange < NC[0] < 1 - IntRange and IntRange < NC[1] < 1 - IntRange):
        for x in C:
            if euclidean2dim(x, NC) < IntRange:
                return True
            
    else:
        t = [0,0]
        
        for x in C:
            if euclidean2dim(x, NC) < IntRange:
                return True        
            
            
            t[0] = x[0] - 1
            t[1] = x[1] + 1
            if euclidean2dim(t, NC) < IntRange:
                return True
        
            t[0] = x[0]
            t[1] = x[1] + 1
            if euclidean2dim(t, NC) < IntRange:
                return True    
        
            t[0] = x[0] + 1
            t[1] = x[1] + 1
            if euclidean2dim(t, NC) < IntRange:
                return True     
        
            t[0] = x[0] - 1
            t[1] = x[1]
            if euclidean2dim(t, NC) < IntRange:
                return True      
        
            t[0] = x[0] + 1
            t[1] = x[1]
            if euclidean2dim(t, NC) < IntRange:
                return True      
        
            t[0] = x[0] - 1
            t[1] = x[1] - 1
            if euclidean2dim(t, NC) < IntRange:
                return True      
        
            t[0] = x[0]
            t[1] = x[1] - 1
            if euclidean2dim(t, NC) < IntRange:
                return True      
        
            t[0] = x[0] + 1
            t[1] = x[1] - 1
            if euclidean2dim(t, NC) < IntRange:
                return True  

    return False
 

###### Subroutine to count the number of overlaps the points in C have with the new point NC
def OverlapCount(C, NC, IntRange, torus_flag = True):

    if IntRange > 1:
        print('Error: Interaction range is bigger than 1')
        exit(0)
        
    NC = NC[0:2]
    count = 0        
    if (torus_flag == False): 
        for x in C:
            if euclidean2dim(x, NC) < IntRange:
                count = count + 1
                
    elif  IntRange < NC[0] < 1 - IntRange and IntRange < NC[1] < 1 - IntRange:
        for x in C:
            if euclidean2dim(x, NC) < IntRange:
                count = count + 1
    
    else:
        t = [0,0]
        for x in C:
            OL = False
            if euclidean2dim(x, NC) < IntRange:
                OL = True       
            
            
            t[0] = x[0] - 1
            t[1] = x[1] + 1
            if euclidean2dim(t, NC) < IntRange:
                OL = True
        
            t[0] = x[0]
            t[1] = x[1] + 1
            if euclidean2dim(t, NC) < IntRange:
                OL = True   
        
            t[0] = x[0] + 1
            t[1] = x[1] + 1
            if euclidean2dim(t, NC) < IntRange:
                OL = True    
        
            t[0] = x[0] - 1
            t[1] = x[1]
            if euclidean2dim(t, NC) < IntRange:
                OL = True      
        
            t[0] = x[0] + 1
            t[1] = x[1]
            if euclidean2dim(t, NC) < IntRange:
                OL = True      
        
            t[0] = x[0] - 1
            t[1] = x[1] - 1
            if euclidean2dim(t, NC) < IntRange:
                OL = True      
        
            t[0] = x[0]
            t[1] = x[1] - 1
            if euclidean2dim(t, NC) < IntRange:
                OL = True      
        
            t[0] = x[0] + 1
            t[1] = x[1] - 1
            if euclidean2dim(t, NC) < IntRange:
                OL = True  
            
            if OL == True:
                count = count + 1
                

    return count

###### Subroutine to find the indecies of the circles in A that are blocking the new circle NC. Used in dCFTP function. See Huber's swap based dCFTP for details
def BlockInd(C, event, IntRange, torus_flag = True):
    NC = event[0:2]
    if IntRange > 1:
        print('Error: Interaction range is bigger than 1')
        exit(0)
        
    Ind = []   
    lC = len(C)
    NC = NC[0:2]
    
    if (torus_flag == False):
        for i in range(lC):
            if euclidean2dim(C[i], NC) < IntRange:
                Ind.append(i)
                
    elif (IntRange < NC[0] < 1 - IntRange and IntRange < NC[1] < 1 - IntRange):
        for i in range(lC):
            if euclidean2dim(C[i], NC) < IntRange:
                Ind.append(i)
                
    else:
        t = [0,0]
        for i in range(lC):
            OL = False
            x = C[i][0:2]
            
            if euclidean2dim(x, NC) < IntRange:
                OL = True       
            
                        
            t[0] = x[0] - 1
            t[1] = x[1] + 1
            if euclidean2dim(t, NC) < IntRange:
                OL = True
        
            t[0] = x[0]
            t[1] = x[1] + 1
            if euclidean2dim(t, NC) < IntRange:
                OL = True   
        
            t[0] = x[0] + 1
            t[1] = x[1] + 1
            if euclidean2dim(t, NC) < IntRange:
                OL = True    
        
            t[0] = x[0] - 1
            t[1] = x[1]
            if euclidean2dim(t, NC) < IntRange:
                OL = True      
        
            t[0] = x[0] + 1
            t[1] = x[1]
            if euclidean2dim(t, NC) < IntRange:
                OL = True      
        
            t[0] = x[0] - 1
            t[1] = x[1] - 1
            if euclidean2dim(t, NC) < IntRange:
                OL = True      
        
            t[0] = x[0]
            t[1] = x[1] - 1
            if euclidean2dim(t, NC) < IntRange:
                OL = True      
        
            t[0] = x[0] + 1
            t[1] = x[1] - 1
            if euclidean2dim(t, NC) < IntRange:
                OL = True  
            
            if OL == True:
                Ind.append(i)

    return Ind
####### Checks whether a given configuration is overlapping or not
def Overlap_check(C, IntRange, torus_flag = True):
    len_C = len(C)
    if len_C ==0:
        return False
    
    for i in range(len_C-1):

        if overlap(C[0:i+1], C[i+1], IntRange, torus_flag) == True:
            return True
           
    return False

#### Removes event from the list A
def remove(A, event):
    for i in range(len(A)):
        if event[2] == A[i][2]:
            del A[i]
            return
    return
## Returns if L == U
def check_stop(L, U):
    len_L = len(L)
    len_U = len(U)
    if len_L < len_U:
        return False
    if len_L > len_U:
        print('Error: L is bigger than U')
    if set(L) == set(U):
        return True
    return False

def update(D, Kappa_0, count_pts):
    no_cirD = len(D)
    X_len = 1
    Y_len = 1
    Kappa = Kappa_0*X_len*Y_len
    
    if  np.random.random_sample() < np.divide(Kappa, Kappa + no_cirD): # True if it is an arrival
        count_pts = count_pts + 1 
        NewCen = (X_len*np.random.random_sample(), Y_len*np.random.random_sample(), count_pts)
        D.append(NewCen)
        event = (-1, NewCen)
    else:
        k = np.random.randint(0, no_cirD) ## Index of the departuring element
        vv = D.pop(k)
        event = (1,vv)  
    return (event, count_pts)                  

############## Ploting the circles in C
def plot_points(C, IntRange):
    for x in C:
        circle1=plt.Circle(x, IntRange/2, color='r')
        plt.gcf().gca().add_artist(circle1) 
        
    plt.show()
    
###### Huber's dominated CFTP method for generating a perfect sample of Strauss process on Cell numbered 'Cell_no' with parameter gamma and intensity Kappa_cell
def dCFTP_swap(Kappa_0, IntRange, torus_flag = True):

    X_len = 1
    Y_len = 1


    Kappa = Kappa_0*X_len*Y_len
    

    events = [] #Records the events happened
    N = np.random.poisson(Kappa)
    D = [(X_len*np.random.random_sample(), Y_len*np.random.random_sample(), i + 1) for i in range(N)] # Intial state 
    


    count_pts = N + 1

    (event, count_pts) = update(D, Kappa_0, count_pts)
    events.append(event)

    n = 1

    while True: ## This while Loop ends only when there is a coalscence
        L = [] # initialization of Lower bounding process
        U = D[:] # initialization of Upper bounding process

        for j in range(n):
            if events[n-j-1][0] == -1: ## It is a deprture in the forward chain (since it is an arrival backward in D)
                remove(L, events[n-j-1][1])
                remove(U, events[n-j-1][1])
                
            else:
                #### One can optimize the following steps 
                Block_Ind_L = BlockInd(L, events[n - j - 1][1], IntRange, torus_flag)
                Block_Ind_U = BlockInd(U, events[n - j - 1][1], IntRange, torus_flag)
                len_BIL = len(Block_Ind_L)
                len_BIU = len(Block_Ind_U)
#                if len_BIU < len_BIL:
#                    print('Error: Upper bound has more spheres than the lower bound\n')
#                    exit(1)
                
                if  len_BIU == 0:
                    L.append(events[n - j - 1][1])
                    U.append(events[n - j - 1][1])
                    
                elif len_BIU == 1:
                    pt = U[Block_Ind_U[0]]
                    remove(L, pt)
                    del U[Block_Ind_U[0]]
                    
                    L.append(events[n - j - 1][1])
                    U.append(events[n - j - 1][1])                        
                    
                elif len_BIL == 0:
                    U.append(events[n - j - 1][1])
                    
                    
                elif len_BIL == 1:

                    del L[Block_Ind_L[0]]                        
                    U.append(events[n - j - 1][1])
            ### Remove this if loop after the debugging
#            if (set(L) <= set(U)) == False:
#                print("L is not a subset of U")
#                exit()
        
        if len(L) == len(U):
            break

        for j in range(n): #Loop for doubling the length of the dominated process
            (event, count_pts) = update(D, Kappa_0, count_pts)
            events.append(event)

        
        n = 2*n            

        
    return (L, count_pts)




###### Model parameters
Kappa = 20
r = 1
eta = 0.25
IntRange = 2*np.divide(r, Kappa**eta)

Itot = 1 ### Number of samples generated

###### Estimation parameters
Est_time = 0
Exp_while_loops = 0
Exp_total_pts = 0
npoints = 0

torus_flag = True ## When the torus_flag is true, the underlying cube is treated as a torus

print('Intensity of PPP, Kappa = ', Kappa, ', \n Interaction Range = ', IntRange,'# iterations = ', Itot)


np.random.seed(0)
print('Program Starting time: ', datetime.datetime.now().time())
print('++++++++++ dominated CFTP M3 for Gibbs HS Process (with swaps) +++++++++++++++')



for b in range(Itot):
    print('---------- DCFTP M3 Iteration: ', b+1, '---------')
    

    (State, count_pts) = dCFTP_swap(Kappa, IntRange, torus_flag)
    npoints = (b/(b+1))*npoints + (1/(b+1))*len(State)
    Exp_total_pts = (b/(b+1))*Exp_total_pts + (1/(b+1))*count_pts
#    print(len(State))
#    print(State)
#    print(Overlap_check(State, IntRange, torus_flag))
    print(' Exp no of points :', npoints, Exp_total_pts)
#    print(' Exp number of pooints generated per a sample:', Exp_total_pts)
#    plot_points(State, IntRange)

#### Conclusions
print('++++++++++ dominated CFTP M3 for Gibbs HS Process (with swaps) +++++++++++++++')
print('\n \n Parameters: \n Kappa_0 = ', Kappa, '\n Interaction Range = ', IntRange,'\n # iterations = ', Itot)
print('Results:')
print('Estimated # points in HS Model:', npoints)
print('Exp number of pooints generated per a sample:', Exp_total_pts)
print('Program Ending time: ', datetime.datetime.now().time())



     