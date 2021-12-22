#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 19 15:08:01 2021
CANAL ABIERTO CIRCULAR
Autor: Carlos Armando De Castro (cadecastro.com)
"""
import numpy as np
import matplotlib.pyplot as plt
D=1.6#[m] Diámetro del tubo
S0=0.0011# Pendiente del ducto
n=0.013#Coef. rugosidad de Manning
Q_des=0.89#[m^3/s] Caudal deseado

#Análisis de flujo gradualmente variado:
y0=1.40#[m] Altura inicial
L=1000#[m] Longitud del canal
Nx=20000#Intervalos análisis FGV.

#Parámetros análisis:
N=1000#Número de pasos curvas normales
g=9.81#[m/s^2] Gravedad

#Parámetros solución numérica niveles:
yn1=D/2#Suposición inicial nivel normal
yc1=D/8#Suposición inicial nivel crítico
yc2=D/4#Suposición inicial nivel crítico
tol=1e-4#[m] Tolerancia absoluta nivel normal
it_max=1000#Iteraciones máximas
#_________________________________________________
#SOLUCIÓN:
pi=np.pi
r=D/2#[m] Radio del ducto
y=np.linspace(D/N,D,N)
#CURVAS EN FUNCIÓN DEL CAUDAL:
theta=np.arcsin((r-y)/r)
A=r*r*(pi/2-theta)-(r-y)*r*np.cos(theta)
P=r*(pi-2*theta)
Rh=A/P
V=np.power(Rh,2/3)*np.power(S0,0.5)/n
Q=A*V
h=A/(2*np.sqrt(r*r-(r-y)*(r-y)))
Fr=V/np.sqrt(g*h)
del theta
del A
del P
del Rh
del h
#PROFUNDIDAD NORMAL:
it=1
cdc=1
while cdc==1:
    theta1=np.arcsin((r-yn1)/r)
    P=r*(pi-2*theta1)
    A=np.power(n*Q_des*np.power(P,2/3)/np.sqrt(S0),3/5)
    yn=r-(r*r*(pi/2-theta1)-A)/(r*np.cos(theta1))
    residual=abs(yn-yn1)
    if residual<=tol or it>=it_max:
       cdc=0      
    it=it+1
    yn1=yn

#Valores de flujo normal:
theta=np.arcsin((r-yn)/r)
An=r*r*(pi/2-theta)-(r-yn)*r*np.cos(theta)
Vn=Q_des/An
h=An/(2*np.sqrt(r*r-(r-yn)*(r-yn)))
Fr_n=Vn/np.sqrt(g*h)
print('--------------------------')
print('VALORES FLUJO NORMAL:')
print('Nivel normal yn = ',np.format_float_positional(yn,precision=3),'m')
print('Velocidad normal = ',np.format_float_positional(Vn,precision=3),'m/s')
print('Número de Froude normal = ',np.format_float_positional(Fr_n,precision=3))
print('Iteraciones totales normal = ',it)
print('Residual absoluto normal = ',np.format_float_scientific(residual,precision=3),'m')
#Dibujo de la tubería a nivel normal:
x1=np.linspace(-r,r,N+1)
y_inf=-np.sqrt(r*r-x1*x1)
y_sup=np.sqrt(r*r-x1*x1)
ysurf=np.array([yn-r,yn-r])
x2=np.array([-np.sqrt(r*r-(r-yn)*(r-yn)),np.sqrt(r*r-(r-yn)*(r-yn))])

#PROFUNDIDAD CRÍTICA:
it=1
cdc=1
while cdc==1:
    A1=r*r*(pi/2-np.arcsin(1-yc1/r))-(r-yc1)*r*np.cos(np.arcsin(1-yc1/r))
    A2=r*r*(pi/2-np.arcsin(1-yc2/r))-(r-yc2)*r*np.cos(np.arcsin(1-yc2/r))
    h1=A1/(2*np.sqrt(r*r-(r-yc1)*(r-yc1)))
    h2=A2/(2*np.sqrt(r*r-(r-yc2)*(r-yc2)))
    f1=pi/2*r*r-r*r*np.arcsin(1-yc1/r)-(r-yc1)/r*np.sqrt(r*r-(r-yc1)*(r-yc1))-Q_des/np.sqrt(g*h1)
    f2=pi/2*r*r-r*r*np.arcsin(1-yc2/r)-(r-yc2)/r*np.sqrt(r*r-(r-yc2)*(r-yc2))-Q_des/np.sqrt(g*h2)
    yc=yc1-(yc1-yc2)*f1/(f1-f2)
    residual=abs(yc-yc1)
    if residual<=tol or it>=it_max:
       cdc=0
    it=it+1
    yc2=yc1
    yc1=yc
Vc=np.sqrt(g*yc)#[m/s]
print('--------------------------')
print('VALORES FLUJO CRÍTICO:')
print('Nivel crítico yc =',np.format_float_positional(yc,precision=3),'m')
print('Velocidad crítica Vc =',np.format_float_positional(Vc,precision=3),'m/s')
print('Iteraciones totales nivel crítico = ',it)
print('Residual absoluto nivel crítico = ',np.format_float_scientific(residual,precision=3),'m')
#Dibujo línea nivel crítico:
y3=np.array([yc-r,yc-r])
x3=np.array([-np.sqrt(r*r-(r-yc)*(r-yc)),np.sqrt(r*r-(r-yc)*(r-yc))])

#FLUJO GRADUALMENTE VARIADO:
dx=L/Nx#[m]
x=np.linspace(0,L,Nx+1)#[m]
yx=np.zeros(len(x))
zpiso=np.zeros(len(x))
Vx=np.zeros(len(x))
Frx=np.zeros(len(x))
del theta
if Fr_n>=1:
    yx[0]=y0
    for i in range(0,len(x)-1):
        theta=np.arcsin(1-yx[i]/r)
        A=r*r*(pi/2-theta)-(r-yx[i])*r*np.cos(theta)
        P=r*(pi-2*theta)
        Rh=A/P
        Vx[i]=Q_des/A
        h=A/(2*np.sqrt(r*r-(r-yx[i])*(r-yx[i])))
        Frx[i]=Vx[i]/np.sqrt(g*h)
        S=n*n*Vx[i]*Vx[i]/np.power(Rh,4/3)
        yx[i+1]=yx[i]+dx*(S0-S)/(1-Frx[i]*Frx[i])
        zpiso[i+1]=zpiso[i]-S0*dx
    theta=np.arcsin((r-yx[len(x)-1])/r)
    A=r*r*(pi/2-theta)-(r-yx[len(x)-1])*r*np.cos(theta)
    Vx[len(x)-1]=Q_des/A
    h=A/(2*np.sqrt(r*r-(r-yx[len(x)-1])*(r-yx[len(x)-1])))
    Frx[len(x)-1]=Vx[len(x)-1]/np.sqrt(g*h)
else:
    yx[len(x)-1]=y0
    for i in range(1,len(x)-1):
        theta=np.arcsin((r-yx[len(x)-i])/r)
        A=r*r*(pi/2-theta)-(r-yx[len(x)-i])*r*np.cos(theta)
        P=r*(pi-2*theta)
        Rh=A/P
        Vx[len(x)-i]=Q_des/A
        h=A/(2*np.sqrt(r*r-(r-yx[len(x)-i])*(r-yx[len(x)-i])))
        Frx[len(x)-i]=Vx[len(x)-i]/np.sqrt(g*h)
        S=n*n*Vx[len(x)-i]*Vx[len(x)-i]/np.power(Rh,4/3)
        yx[len(x)-i-1]=yx[len(x)-i]-dx*(S0-S)/(1-Frx[len(x)-i]*Frx[len(x)-i])
        zpiso[i+1]=zpiso[i]-S0*dx
    yx[0]=yx[1]-dx*(S0-S)/(1-Frx[1]*Frx[1])
del theta
del A
del h
theta=np.arcsin(1-yx/r)
A=r*r*(pi/2-theta)-(r-yx)*r*np.cos(theta)
Vx=Q_des/A
h=A/(2*np.sqrt(r*r-(r-yx)*(r-yx)))
Frx=Vx/np.sqrt(g*h)

zsurf=yx+zpiso
zn=yn+zpiso
zc=yc+zpiso
zborde=D+zpiso

#GRÁFICAS:
plt.figure(1)
plt.plot(Q,y,'b')
plt.plot(Q_des,yn,'r*')
plt.grid(True,'both','both')
plt.title('Niveles canal circular')
plt.ylabel('y [m]')
plt.xlabel('Q [m^3/s]')

plt.figure(2)
plt.plot(Q,V,'b')
plt.plot(Q_des,Vn,'r*')
plt.grid(True,'both','both')
plt.title('Velocidades canal circular')
plt.ylabel('V [m/s]')
plt.xlabel('Q [m^3/s]')

plt.figure(3)
plt.plot(Q,Fr,'b')
plt.plot(Q_des,Fr_n,'r*')
plt.grid(True,'both','both')
plt.title('Nº de Froude canal circular')
plt.ylabel('Fr')
plt.xlabel('Q [m^3/s]')

plt.figure(4)
plt.plot(x2,ysurf,'b')
plt.plot(x3,y3,'r--')
plt.plot(x1,y_sup,'k')
plt.plot(x1,y_inf,'k')
plt.title('Sección transversal flujo normal')
plt.ylabel('y [m]')
plt.xlabel('x [m]')
plt.axis('square')
plt.xlim([-1.05*r,r*1.05])
plt.ylim([-1.05*r,r*1.05])
plt.legend(['Nivel normal','Nivel crítico'])

plt.figure(5)
plt.subplot(311)
plt.plot(x,zborde,'k')
plt.plot(x,zsurf,'b')
plt.plot(x,zn,'g-.')
plt.plot(x,zc,'r--')
plt.plot(x,zpiso,'k')
plt.title('Nivel, velocidad y Fr a lo largo del canal circular')
plt.ylabel('z [m]')
plt.legend(['Techo canal','Superficie agua','Nivel normal','Nivel crítico','Piso canal'])
plt.subplot(312)
plt.plot(x,Vx,'b')
plt.ylabel('V [m/s]')
plt.subplot(313)
plt.plot(x,Frx,'b')
plt.xlabel('x [m] - cadecastro.com')
plt.ylabel('Fr')