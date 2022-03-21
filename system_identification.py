# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 11:13:38 2021

@author: 502001
"""

import numpy as np
import matplotlib.pyplot as plt


samptime=4000
sampfreq=1

#      a
#----------------
#  b  +  c * s
def _1zi_okure(a,b,c,t,_in,inzm1,outzm1):
    return (_in*a*t+inzm1*a*t-b*t*outzm1+2*outzm1*c)/(b*t+2*c)


t=np.linspace(0,samptime,(samptime*sampfreq)+1)

step=np.zeros((samptime*sampfreq)+1).astype(float)
step[2:]=135
step[:2]=0
deviation=np.zeros((samptime*sampfreq)+1).astype(float)

Gin=np.zeros((samptime*sampfreq)+1).astype(float)
output_1zi=np.zeros((samptime*sampfreq)+1).astype(float)
output=np.zeros((samptime*sampfreq)+1).astype(float)
output_solve=np.zeros((samptime*sampfreq)+1).astype(float)
#step[sampfreq:]=1.0
#step[:sampfreq-1]=0.0
M_sig_elem=8
M_signal=np.zeros(M_sig_elem).astype(int)
M_signal[:]=1

if(1):
    b1=4.504e-05
    b2=0.0005987100000000001
    a11=-1.90024278
    a12=0.15881931999999999
    a21=0.0116751
    a22=-0.0009749400000000001
if(0):    
    b1=0.0
    b2=1.0

    a11=0.0
    a12=1.0
    a21=0.0
    a22=0.0


A=[[a11,a12],
   [a21,a22]]
b=[1,0]
c=[[0],[1]]
delta_x1=np.zeros((samptime*sampfreq)+1).astype(float)
delta_x2=np.zeros((samptime*sampfreq)+1).astype(float)
x1=np.zeros((samptime*sampfreq)+1).astype(float)
x2=np.zeros((samptime*sampfreq)+1).astype(float)

power=np.array([1e+2,1e+1,1e+0,1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8])
mult=np.array([-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9])

def state_system(u,a11,a12,a21,a22,b1,b2,x1,x2):
    delta_x1=a11*x1+a12*x2+u*b1
    delta_x2=a21*x1+a22*x2+u*b2
    return delta_x1,delta_x2


fdbk_bit=0
Msig=0
for i in range(1,((samptime*sampfreq)+1),1):
    M_signal[1:]=M_signal[:-1]
    for j in range(1,M_sig_elem):
        M_signal[0]^=M_signal[j]
    Gin[i]=M_signal[0]
Gin=step
for i in range(1,((samptime*sampfreq)+1),1):    
    #output_1zi[i]=_1zi_okure(1,1,1,1/sampfreq,Gin[i],Gin[i-1],output_1zi[i-1])
   # output[i]=_1zi_okure(1,1,1,1/sampfreq,output_1zi[i],output_1zi[i-1],output[i-1])
   
   output[i]=_1zi_okure(1,0.0418,19950,1/sampfreq,Gin[i],Gin[i-1],output[i-1])

def system_solve(a11,a12,a21,a22,b1,b2):
    x1_sum=0
    x2_sum=0  
    for i in range(1,((samptime*sampfreq)+1),1):
        delta_x1[i]=state_system(Gin[i],a11,a12,a21,a22,b1,b2,x1[i-1],x2[i-1])[0]
        delta_x2[i]=state_system(Gin[i],a11,a12,a21,a22,b1,b2,x1[i-1],x2[i-1])[1]
        x1_sum+=delta_x1[i]/sampfreq #integral
        x2_sum+=delta_x2[i]/sampfreq #integral
        x1[i]=x1_sum
        x2[i]=x2_sum    
    return x1

if(0):
    for k in range(200):
        if(1):  #b1
            final_value = 0
            pow2_min=10e+20
            pow2=0
            pv=0
            for j in power:
                for i in mult:
                    b1=i*j+pv
                    pow2=np.sum((output-system_solve(a11,a12,a21,a22,b1,b2))**2)
                    
                    if(pow2_min > pow2) :   
                        pow2_min=pow2
                        final_value=b1
                b1=final_value
                pv=b1
            
        if(1):  #b2
            final_value = 0
            pow2_min=10e+20
            pow2=0
            pv=0
            for j in power:
                for i in mult:
                    b2=i*j+pv
                    pow2=np.sum((output-system_solve(a11,a12,a21,a22,b1,b2))**2)
                    
                    if(pow2_min > pow2) :   
                        pow2_min=pow2
                        final_value=b2
                b2=final_value
                pv=b2
            
            
        if(1):    #a11
            final_value = 0
            pow2_min=10e+20
            pow2=0
            pv=0
            for j in power:
                for i in mult:
                    a11=i*j+pv
                    pow2=np.sum((output-system_solve(a11,a12,a21,a22,b1,b2))**2)
                    
                    if(pow2_min > pow2) :   
                        pow2_min=pow2
                        final_value=a11
                a11=final_value
                pv=a11
        if(1):   #a12
            final_value = 0
            pow2_min=10e+20
            pow2=0
            pv=0
            for j in power:
                for i in mult:
                    a12=i*j+pv
                    pow2=np.sum((output-system_solve(a11,a12,a21,a22,b1,b2))**2)
                    
                    if(pow2_min > pow2) :   
                        pow2_min=pow2
                        final_value=a12
                a12=final_value
                pv=a12
        if(1):    #a21
            final_value = 0
            pow2_min=10e+20
            pow2=0
            pv=0
            for j in power:
                for i in mult:
                    a21=i*j+pv
                    pow2=np.sum((output-system_solve(a11,a12,a21,a22,b1,b2))**2)
                    
                    if(pow2_min > pow2) :   
                        pow2_min=pow2
                        final_value=a21
                a21=final_value
                pv=a21
        if(1):    #a22
            final_value = 0
            pow2_min=10e+20
            pow2=0
            pv=0
            for j in power:
                for i in mult:
                    a22=i*j+pv
                    pow2=np.sum((output-system_solve(a11,a12,a21,a22,b1,b2))**2)
                    
                    if(pow2_min > pow2) :   
                        pow2_min=pow2
                        final_value=a22
                a22=final_value
                pv=a22
        print('pow2='+str(pow2))
        print('k='+str(k))
    
#output=np.append([0],output)
#plt.plot(t,step)
#plt.plot(t,Gin)
x1=system_solve(a11,a12,a21,a22,b1,b2)
plt.plot(t,output)
plt.plot(t,x1)
plt.show()
print('b1='+str(b1))
print('b2='+str(b2))
print('a11='+str(a11))
print('a12='+str(a12))
print('a21='+str(a21))
print('a22='+str(a22))
#print(output[sampfreq*6000])