#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 17:57:51 2021
@author: cadecastro.com
"""
import numpy as np
Q=float(input('Caudal [m³/s]='))
D=float(input('Diámetro ducto [mm]='))
L=float(input('Longitud ducto [m]='))
rho=float(input('Densidad fluido [kg/m³]='))
mu=float(input('Viscosidad dinámica [Pa*s]='))
ks=float(input('Rugosidad pared [mm]='))
#Parámetros:
D=D/1000 #Pasar D a metros
V=4*Q/(np.pi*D*D)
print("RESULTADOS:")
print('V =',np.format_float_positional(V,precision=3),'m/s')
Re=rho*V*D/mu##Número de Reynolds
print('Re =',np.format_float_scientific(Re,precision=3))
e_d=0.001*ks/D##Rugosidad relativa
print('Rugosidad relativa =',np.format_float_scientific(e_d,precision=3))
#Solución numérica:
if Re<2300:
    f=64/Re#
    print('Flujo laminar.')
    print('f =',np.format_float_positional(f,precision=5))
else:
   itmax=50#
   emax=1e-6#
   it=0#
   cdc=1#
   f0=0.015#
   while cdc==1:
       it=it+1#
       f=np.power(-2*np.log10(e_d/3.7+2.51/(Re*np.sqrt(f0))),-2)#
       res=abs(f-f0)#
       if res<=emax or it>=itmax:
           cdc=0#
       else:
           f0=f#
   print('Flujo turbulento.')
   print('f =',np.format_float_positional(f,precision=5))
   print('Información solución numérica:')
   print('Residual absoluto=',np.format_float_scientific(res,precision=3))
   print('Iteraciones =',it)
print('Pérdida por fricción hf=',np.format_float_positional(f*L*V*V/(2*9.81*D),precision=2),'m')