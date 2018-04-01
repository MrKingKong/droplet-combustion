# -*- coding: utf-8 -*-
"""
Created on Sat Mar 31 10:11:24 2018

@author: wilson
"""

import cantera as ct

gas = ct.Solution('gri30.cti')
gas2 = ct.Solution('gri30.cti')
r = ct.ConstPressureReactor(gas)
r.volume = 1
env = ct.Reservoir(gas2)
w = ct.Wall(r, env)
w.set_heat_flux(2)# w.set_heat_flux is the function to use

net = ct.ReactorNet({r})

print(r.thermo.h*r.mass)
net.advance(1)
print(r.thermo.h*r.mass)