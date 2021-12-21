#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 18:30:06 2021

@author: peterspencer-smith
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps

c=3E8
e=8.85418782E-12
k=(2*np.pi)/(1E-6)
h=(2E-5)/100
l=1e-6
        
#Part One code  
MyInput = '0'
while MyInput != 'q':
    MyInput = input('In this exercise the diffraction intensity pattern resulting from a 1-D single slit and 2-D \
square aperture is to be modelled. Please enter one of the following choices, "a" for the 1-D case, "b" for the 2-D case or "q" to quit: ')
    print('User has entered the choice: ',MyInput)

    if MyInput == 'a':
        
        input_X= input('Please enter a value for X, the aperture width where X also sets the number of iterations \
  over which the aperture is integrated (an integer of around 100) and X scaled by 2e-5 m: ')
        X=int(input_X)
        
        input_z= input('Please enter a value for Z, the distance to the screen in meters (a floating point \
number of around 2e-2): ')
        z=float(input_z)
        
        a=k/(2*np.pi*z)

        x=np.arange(-0.005,0.00505,5E-5)
        
        F1first=np.cos((k/2*z)*(x)**2)
        
        
        def F1odd(x):
            F1odd=0
            for i in np.arange(1,X,2):
                F1odd=F1odd+np.cos((k/(2*z))*(x-(i*h))**2)
            return (4*F1odd)
            
        def F1even(x):
            F1even=0
            for i in np.arange(2,X+1,2):
                F1even=F1even+np.cos((k/(2*z))*(x-(i*h))**2)
            return (2*F1even)
            
        F1last=np.cos((k/2*z)*(x-(X*h))**2)
            
        Re=(a*h/3)*(F1first+F1odd(x)+F1even(x)+F1last)
            
            
        F2first=np.sin((k/2*z)*(x)**2)
            
        def F2odd(x):
            F2odd=0
            for i in np.arange(1,X,2):
                F2odd=F2odd+np.sin((k/(2*z))*(x-(i*h))**2)
            return (4*F2odd)
            
        def F2even(x):
            F2even=0
            for i in np.arange(2,X+1,2):
                F2even=F2even+np.sin((k/(2*z))*(x-(i*h))**2)
            return (2*F2even)
            
        F2last=np.sin((k/2*z)*(x-(X*h))**2)
            
        Im=(a*h/3)*(F2first+F2odd(x)+F2even(x)+F2last)
             
        I = (e*c*((Re)**2 +(Im)**2))
        plt.plot(x,I)
        plt.ylabel('Relative Intensity')
        plt.xlabel('Screen position (m)')
        plt.show()

#Part Two code    
    elif MyInput == 'b':
        
        input_x_aperture= input('Please enter a value for x_aperture (m) ( a floating point number of around 2e-5): ')
        x_aperture=float(input_x_aperture)
        delta_x = x_aperture/100
        
        input_y_aperture= input('Please enter a value for y_aperture (m) ( a floating point number of around 2e-5): ')
        y_aperture=float(input_y_aperture)
        delta_y = y_aperture/100

        input_z= input('Please enter a value for Z, the distance to the screen in meters (a floating point number of around 2e-2): ')
        z=float(input_z)

        a=k/(2*np.pi*z)

        x_scale = 1e-4
        y_scale = 1e-4
#        x and y translations are integers that represent the percentage of the screen eg: 50 = 50% translation. x and y scales are both in meters.        
        x_translation = -50
        y_translation = -50
        
        
        def simp(x,y,z):
            x_answer=0
            x_iter = np.arange(0,x_aperture+delta_x,delta_x)
            x_total = []
            for xp in x_iter:
                x_total.append(np.exp(((complex(0,1)*k)/(2*z))*(((x+(x_translation))*x_scale)-(xp))**2))
            x_answer = simps(x_total)
            y_answer=0
            y_iter = np.arange(0,y_aperture+delta_y,delta_y)
            y_total = []
            for yp in y_iter:
                y_total.append(np.exp(((complex(0,1)*k)/(2*z))*(((y+(y_translation))*y_scale)-(yp))**2))
            y_answer = simps(y_total)
            
            Z=a*x_answer*y_answer
            return abs(e*c*Z*Z.conjugate())
        
        E = np.zeros((100,100))
        
        for x in np.arange(0,100,1):
            for y in np.arange(0,100,1):
                E[x,y] = simp(x,y,z)
        
        plt.imshow(E)
        plt.ylabel('y Screen position 1e-4(m)')
        plt.xlabel('x Screen position 1e-4(m)')
        plt.show()
        


    elif MyInput != 'q':
        print('Sorry, this is not a recognised choice')
        
print('Goodbye dear friend, just remember, if things get a little hairy.. everything is Physics, but Physics isnt everything.')

