#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22/2021
@author: cdc
"""
import numpy as np
Re=float(input('Número de Reynolds='))
e_d=float(input('Rugosidad relativa='))
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