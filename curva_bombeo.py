#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 12:43:45 2021

@author: cadecastro.com
"""
import numpy as np
import matplotlib.pyplot as plt
#Caudal máximo:
Qmax=0.050#[m^3/s] Caudal máximo a simular
N=int(100)#Puntos de caudal a simular

#Geometría ducto:
D=0.1016#[m] Diámetro (hidráulico si no es circular)
L=50.0#[m] Longitud ducto
z=10.0#[m] Altura a elevar fluido

#Propiedades fluido:
rho=999.8#[kg/m^3] Densidad
mu=0.001#[Pa*s] Viscosidad dinámica
g=9.81#[m/s^2] Aceleración de la gravedad

#Material ducto:
ks=0.5#[mm] Rugosidad

#Pérdidas locales:
K=3.5#Suma de coef. de pérdidas locales en el ducto

#Curva de la bomba (elíptica)
H_max=25 #[m] Cabeza máxima (estática)
Q_max=0.050 #[m^3/s] Caudal máximo bomba (sin carga)

#SOLUCIÓN:
Q=np.linspace(Qmax/N,Qmax,N)#
hreq=np.zeros(N)#
hp=np.zeros(N)#
Preq=np.zeros(N)#
Pp=np.zeros(N)#
e_d=ks*0.001/D#Rugosidad relativa

for i in range(0,N):
    V=4*Q[i]/(np.pi*D*D)#[m/s] Velocidad media flujo
    Re=rho*V*D/mu#Número de Reynolds

    #Pérdidas por fricción:
    if Re<2300:
        f=64/Re
    else:
        itmax=50#
        emax=1e-6#
        it=0#
        cdc=1#
        f0=float(0.015)#
        while cdc==1:
            it=it+1#
            f=np.power(-2*np.log10(e_d/3.7+2.51/(Re*np.sqrt(f0))),-2)
            res=abs(f-f0)
            if res<=emax or it>=itmax:
                cdc=0
            else:
                f0=f
    hL=f*(L/D)*V*V/(2*g)+K*V*V/(2*g)#[m] Pérdidas
    hreq[i]=z+hL#[m] Cabeza de bomba requerida
    hp[i]=H_max*np.sqrt(1-Q[i]*Q[i]/(Q_max*Q_max))#[m] Cabeza entregada por la bomba
    Preq[i]=hreq[i]*rho*g*Q[i]/1000#[kW] Potencia requerida
    Pp[i]=hp[i]*rho*g*Q[i]/1000#[kW] Potencia de la bomba
plt.figure(1)
plt.plot(Q,hreq,'b-')
plt.plot(Q,hp,'r-')
plt.grid(True,'both','both')
plt.legend(['Cabeza requerida','Cabeza de la bomba'])
plt.xlabel('Q [m^3/s]')
plt.ylabel('h [m]')
plt.title('Cabeza de bombeo')
plt.figure(2)
plt.plot(Q,Preq,'b-')
plt.plot(Q,Pp,'r-')
plt.grid(True,'both','both')
plt.legend(['Potencia requerida','Potencia de la bomba'])
plt.xlabel('Q [m^3/s]')
plt.ylabel('P [kW]')
plt.title('Potencia de bombeo')