#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 12:42:14 2021
CANAL ABIERTO RECTANGULAR
Autor: Carlos Armando De Castro (cadecastro.com)
"""
import numpy as np
import matplotlib.pyplot as plt
#Datos canal rectangular:
b=5.0#[m] Base del canal
h=2.0#[m] Altura del canal
S0=0.001000# Pendiente del canal
n=0.013#Coef. rugosidad de Manning
Q_des=10.0#[m^3/s] Caudal deseado

#Análisis de flujo gradualmente variado:
y0=1.5#[m] Altura inicial
L=1000#[m] Longitud del canal
Nx=10000#Intervalos análisis FGV.

#Parámetros análisis:
N=1000#Número de pasos curvas normales
g=9.81#[m/s^2] Gravedad

#SOLUCIÓN:
#Parámetros solución numérica niveles:
yn1=h/2#Suposición inicial nivel normal
tol=1e-4#[m] Tolerancia absoluta nivel normal
it_max=1000#Iteraciones máximas
y=np.linspace(h/N,h,N)
V=np.zeros(N)
Q=np.zeros(N)
Fr=np.zeros(N)

#CURVAS EN FUNCIÓN DEL CAUDAL:
A=b*y
P=b+2*y
Rh=A/P
V=np.power(Rh,2/3)*np.sqrt(S0)/n
Q=A*V
Fr=V/np.sqrt(g*y)
del A
del P
del Rh

#PROFUNDIDAD NORMAL:
it=1
cdc=1
while cdc==1:
    P=b+2*yn1
    A=np.power(n*Q_des*np.power(P,2/3)/np.sqrt(S0),3/5)
    yn=A/b
    residual=abs(yn-yn1)
    if residual<=tol or it>=it_max:
       cdc=0
    it=it+1
    yn1=yn

#Valores de flujo normal:
An=b*yn
Vn=Q_des/An
Fr_n=Vn/np.sqrt(g*yn)
print('--------------------------')
print('VALORES FLUJO NORMAL:')
print('Nivel normal yn =',np.format_float_positional(yn,precision=3),'m')
print('Velocidad normal = ',np.format_float_positional(Vn,precision=3),'m/s')
print('Número de Froude normal = ',np.format_float_positional(Fr_n,precision=3))
print('Iteraciones totales normal = ',it)
print('Residual absoluto normal =',np.format_float_scientific(residual,precision=3),'m')  
#Dibujo del canal nivel normal:
x1=np.array([-b/2,-b/2,b/2,b/2])
y_wall=np.array([h,0,0,h])
x2=np.array([-b/2,b/2])
ysurf=np.array([yn,yn])

#PROFUNDIDAD CRÍTICA:
yc=np.power(Q_des/(b*np.sqrt(g)),2/3)#[m]
Vc=np.sqrt(g*yc)#[m/s]
print('--------------------------')
print('VALORES FLUJO CRÍTICO:')
print('Nivel crítico yc =',np.format_float_positional(yc,precision=3),'m')
print('Velocidad crítica Vc [m/s] =',np.format_float_positional(Vc,precision=3),'m/s')
y3=np.array([yc,yc])
#FLUJO GRADUALMENTE VARIADO:
dx=L/Nx#[m]
x=np.linspace(0,L,Nx+1)#[m]
yx=np.zeros(len(x))
zpiso=np.zeros(len(x))
Vx=np.zeros(len(x))
Frx=np.zeros(len(x))
if Fr_n>=1:
    yx[0]=y0
    for i in range(0,len(x)-1):
        A=b*yx[i]
        P=b+2*yx[i]
        Rh=A/P
        Vx[i]=Q_des/A
        Frx[i]=Vx[i]/np.sqrt(g*yx[i])
        S=n*n*Vx[i]*Vx[i]/np.power(Rh,4/3)
        yx[i+1]=yx[i]+dx*(S0-S)/(1-Frx[i]*Frx[i])
        zpiso[i+1]=zpiso[i]-S0*dx
else:
    yx[len(x)-1]=y0
    for i in range(1,len(x)-1):
        A=b*yx[len(x)-i]
        P=b+2*yx[len(x)-i]
        Rh=A/P
        Vx[len(x)-i]=Q_des/A
        Frx[len(x)-i]=Vx[len(x)-i]/np.sqrt(g*yx[len(x)-i])
        S=n*n*Vx[len(x)-i]*Vx[len(x)-i]/np.power(Rh,4/3)
        yx[len(x)-i-1]=yx[len(x)-i]-dx*(S0-S)/(1-Frx[len(x)-i]*Frx[len(x)-i])
        zpiso[i+1]=zpiso[i]-S0*dx
    yx[0]=yx[1]-dx*(S0-S)/(1-Frx[1]*Frx[1])
del A
A=b*yx
Vx=Q_des/A
Frx=Vx/np.sqrt(g*yx)

zsurf=yx+zpiso
zn=yn+zpiso
zc=yc+zpiso
zborde=h+zpiso

#GRÁFICAS:
plt.figure(1)
plt.plot(Q,y,'b')
plt.plot(Q_des,yn,'r*')
plt.grid(True,'both','both')
plt.title('Niveles canal rectangular')
plt.ylabel('y [m]')
plt.xlabel('Q [m^3/s]')

plt.figure(2)
plt.plot(Q,V,'b')
plt.plot(Q_des,Vn,'r*')
plt.grid(True,'both','both')
plt.title('Velocidades canal rectangular')
plt.ylabel('V [m/s]')
plt.xlabel('Q [m^3/s]')

plt.figure(3)
plt.plot(Q,Fr,'b')
plt.plot(Q_des,Fr_n,'r*')
plt.grid(True,'both','both')
plt.title('Nº de Froude canal rectangular')
plt.ylabel('Fr')
plt.xlabel('Q [m^3/s]')

plt.figure(4)
plt.plot(x2,ysurf,'b')
plt.plot(x2,y3,'r--')
plt.plot(x1,y_wall,'k')
plt.title('Sección transversal flujo normal')
plt.ylabel('y [m]')
plt.xlabel('x [m]')
plt.xlim([-1.1*b/2,1.1*b/2])
plt.ylim([-0.1*h,1.1*h])
plt.axis('equal')
plt.legend(['Nivel normal','Nivel crítico'])

plt.figure(5)
plt.subplot(311)
plt.plot(x,zborde,'k')
plt.plot(x,zsurf,'b')
plt.plot(x,zn,'g-.')
plt.plot(x,zc,'r--')
plt.plot(x,zpiso,'k')
plt.title('Nivel, velocidad y Fr a lo largo del canal rectangular')
plt.ylabel('z [m]')
plt.legend(['Borde del canal','Superficie agua','Nivel normal','Nivel crítico','Piso canal'])

plt.subplot(312)
plt.plot(x,Vx,'b')
plt.ylabel('V [m/s]')

plt.subplot(313)
plt.plot(x,Frx,'b')
plt.xlabel('x [m] - cadecastro.com')
plt.ylabel('Fr')