# -*- coding: utf-8 -*-
"""
Created on Sun Apr  1 13:45:59 2018

@author: wilson
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Apr  1 13:14:15 2018

@author: wilson
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Apr  1 10:24:05 2018

@author: wilson
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
from scipy.integrate import odeint
from droplet import *
import matplotlib.pyplot as plt

def dT_dt(T,t,droplet, Tenv):#heat transfer to droplet, droplet heating and vaporization    
    Xsurf = np.exp(-droplet.hfg/8315.0*droplet.MW*(1.0/T - 1.0/droplet.ABP))
    Ysurf = Xsurf*droplet.MW/(droplet.MW + 28.85)
    B = Ysurf/(1 - Ysurf)
    mdot = 4.0*np.pi*droplet.k*droplet.r/droplet.cp*np.log(1.0+B)
    
    #calculate the heat transfer coefficient
    Pr = 0.7
    #calculate Grashof #
    beta = 0.003695 #thermal expansion coefficient, 1/K ??of air or fuel droplet??
    L = 2*droplet.r #droplet diameter, m
    g = 9.81 #gravitational acceleration, m/s^2
    rho = 0.25 #density of air, approx., kg/m^3
    del_T = Tenv - droplet.T
    nu = 530e-7 #air viscosity?
    Gr = L**3*rho**2*g*del_T*beta/nu**2

    #Rayleigh # based on Grashof and Prandtl #'s
    Ra = Gr*Pr

    #Nusselt #, -- T.Yuge
    try:
        Nu = 2 + 0.43*Ra**0.25
    except:
        print('error occured at finding nusselt number')
    #convective heat transfer coefficient
    k = 91e-3 #thermal conductivity for air at ?1600K?
    h = Nu*k/L
    
    #establish the ordinary differential equation
    A = 4.0*np.pi*droplet.r**2.0
    V = 4.0/3.0*np.pi*droplet.r**3.0
    dTdt = (h*A*(Tenv - droplet.T) - mdot*droplet.hfg)/(droplet.rho*V*droplet.cp)
    return dTdt


def vaporize(T_kernel, d_array, dt):
    qdot = 0
    mdot = 0
    
    
    for droplet in d_array:
        #calculate the wetbulb temperature
        BT = droplet.cp*(T_kernel-droplet.ABP)/droplet.hfg
        YS = BT/(1+BT)
        XS = YS*28.85/(droplet.MW - YS*(droplet.MW - 28.85))
        T_wetbulb = 1/(1.0/droplet.ABP - 8314/droplet.MW/droplet.hfg*np.log(XS))
        
        if droplet.T > T_wetbulb:
            print('set temp')
            droplet.T = T_wetbulb
            pass
     
        if (droplet.T < T_wetbulb):
            ts = np.linspace(0,dt,2)
            Ts = odeint(dT_dt,droplet.T,ts,args=(droplet, T_kernel))
            Tavg = (droplet.T + Ts[1])/2
            Xsurf = np.exp(-droplet.hfg/8315.0*droplet.MW*(1.0/Tavg - 1.0/droplet.ABP))
            Ysurf = Xsurf*droplet.MW/(droplet.MW + 28.85)
            B = Ysurf/(1 - Ysurf)
            mdot += 4.0*np.pi*droplet.k*droplet.r/droplet.cp*np.log(1.0+B)
             #calculate the heat transfer coefficient
            Pr = 0.7
            #calculate Grashof #
            beta = 0.003695 #thermal expansion coefficient, 1/K ??of air or fuel droplet??
            L = 2*droplet.r #droplet diameter, m
            g = 9.81 #gravitational acceleration, m/s^2
            rho = 0.25 #density of air, approx., kg/m^3
            del_T = T_kernel - droplet.T
            nu = 530e-7 #air viscosity?
            Gr = L**3*rho**2*g*del_T*beta/nu**2
        
            #Rayleigh # based on Grashof and Prandtl #'s
            Ra = Gr*Pr
        
            #Nusselt #, -- T.Yuge
            try:
                Nu = 2 + 0.43*Ra**0.25
            except:
                print('error occured at finding nusselt number')
            #convective heat transfer coefficient
            k = 91e-3 #thermal conductivity for air at ?1600K?
            h = Nu*k/L
            A = 4.0*np.pi*droplet.r**2.0
            qdot += h*A*(T_kernel - Tavg)
            droplet.T = Ts[1]
            droplet.m = droplet.m - 4.0*np.pi*droplet.k*droplet.r/droplet.cp*np.log(1.0+B)*dt
            droplet.r = (droplet.m*3/4/np.pi/droplet.rho)**(1.0/3.0)
            pass
        

            
        if droplet.T == T_wetbulb:
            print('in vaporization')
            Xsurf = np.exp(-droplet.hfg/8315.0*droplet.MW*(1.0/droplet.T - 1.0/droplet.ABP))
            Ysurf = Xsurf*droplet.MW/(droplet.MW + 28.85)
            B = Ysurf/(1 - Ysurf)
            
            #calculate the heat transfer coefficient
            Pr = 0.7
            #calculate Grashof #
            beta = 0.003695 #thermal expansion coefficient, 1/K ??of air or fuel droplet??
            L = 2*droplet.r #droplet diameter, m
            g = 9.81 #gravitational acceleration, m/s^2
            rho = 0.25 #density of air, approx., kg/m^3
            del_T = T_kernel - droplet.T
            nu = 530e-7 #air viscosity?
            Gr = L**3*rho**2*g*del_T*beta/nu**2
        
            #Rayleigh # based on Grashof and Prandtl #'s
            Ra = Gr*Pr
        
            #Nusselt #, -- T.Yuge
            try:
                Nu = 2 + 0.43*Ra**0.25
            except:
                print('error occured at finding nusselt number')
            #convective heat transfer coefficient
            k = 91e-3 #thermal conductivity for air at ?1600K?
            h = Nu*k/L
            A = 4.0*np.pi*droplet.r**2.0
            qdot += h*A*(T_kernel - droplet.T)   
            droplet.m = droplet.m - qdot/droplet.hfg*dt
            if droplet.m < 0:
                pass
            else:
                droplet.r = (droplet.m*3.0/4.0/np.pi/droplet.rho)**(1.0/3.0)    
                mdot += qdot/droplet.hfg
            pass
               
        if droplet.m < 0:
            d_array.remove(droplet)
            pass
    return mdot, qdot