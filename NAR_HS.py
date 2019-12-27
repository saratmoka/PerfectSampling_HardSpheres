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
 


####### Checks whether a given configuration is overlapping or not
def Overlap_check(C, IntRange, torus_flag = True):
    len_C = len(C)
    if len_C ==0:
        return False
    
    for i in range(len_C-1):

        if overlap(C[0:i+1], C[i+1], IntRange, torus_flag) == True:
            return True
           
    return False




               

############## Ploting the circles in C
def plot_points(C, IntRange):
    for x in C:
        circle1=plt.Circle(x, IntRange/2, color='r')
        plt.gcf().gca().add_artist(circle1) 
        
    plt.show()
    
###### Huber's dominated CFTP method for generating a perfect sample of Strauss process on Cell numbered 'Cell_no' with parameter gamma and intensity Kappa_cell
def NAR(Kappa_0, IntRange, torus_flag = True):

    X_len = 1
    Y_len = 1


    Kappa = Kappa_0*X_len*Y_len
    
    count_pts = 0
    while True:

        N = np.random.poisson(Kappa)
        count_pts += N
        C = [(X_len*np.random.random_sample(), Y_len*np.random.random_sample(), i + 1) for i in range(N)] # Intial state 
        if not Overlap_check(C, IntRange, torus_flag):
            return (C, count_pts)
            




###### Model parameters
Kappa = 55
r = 0.5
eta = 0.75
IntRange = 2*np.divide(r, Kappa**eta)

Itot = 200 ### Number of samples generated

###### Estimation parameters
Est_time = 0
Exp_while_loops = 0
Exp_total_pts = 0
npoints = 0

torus_flag = True ## When the torus_flag is true, the underlying cube is treated as a torus

print('Intensity of PPP, Kappa = ', Kappa, ', \n Interaction Range = ', IntRange,'# iterations = ', Itot)


np.random.seed(0)
print('Program Starting time: ', datetime.datetime.now().time())
print('++++++++++ NAR for Gibbs HS Process  +++++++++++++++')



for b in range(Itot):
    print('---------- NAR Iteration: ', b+1, '---------')
    

    (State, count_pts) = NAR(Kappa, IntRange, torus_flag)
    npoints = (b/(b+1))*npoints + (1/(b+1))*len(State)
    Exp_total_pts = (b/(b+1))*Exp_total_pts + (1/(b+1))*count_pts
#    print(len(State))
#    print(State)
#    print(Overlap_check(State, IntRange, torus_flag))
    print(' Exp no of points :', npoints, Exp_total_pts)
#    print(' Exp number of pooints generated per a sample:', Exp_total_pts)
#    plot_points(State, IntRange)

#### Conclusions
print('++++++++++ NAR for Gibbs HS Process +++++++++++++++')
print('\n \n Parameters: \n Kappa_0 = ', Kappa, '\n Interaction Range = ', IntRange,'\n # iterations = ', Itot)
print('Results:')
print('Estimated # points in HS Model:', npoints)
print('Exp number of pooints generated per a sample:', Exp_total_pts)
print('Program Ending time: ', datetime.datetime.now().time())



     