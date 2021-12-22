#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 17:57:13 2021
CANAL ABIERTO TRAPEZOIDAL
Autor: Carlos Armando De Castro (cadecastro.com)
"""
import numpy as np
import matplotlib.pyplot as plt
#Datos canal circular:
b=2.0#[m] Base del canal
W=2.0#[m] Largo sección inclinada
theta=30#[grados] Inclinación pared lateral desde horizontal
S0=0.001050# Pendiente del canal
n=0.013#Coef. rugosidad de Manning
Q_des=2.0#[m^3/s] Caudal deseado

#Análisis de flujo gradualmente variado:
y0=0.95#[m] Altura inicial
L=1000#[m] Longitud del canal
Nx=20000#Intervalos análisis FGV.

#Parámetros análisis:
N=1000#Número de pasos curvas normales
g=9.81#[m/s^2] Gravedad

#SOLUCIÓN:
theta=theta*np.pi/180
h=W*np.sin(theta)
#Parámetros solución numérica niveles:
yn1=h/2#Suposición inicial nivel normal
yc1=h/4#Suposición inicial nivel crítico
tol=1e-4#[m] Tolerancia absoluta niveles
it_max=1000#Iteraciones máximas
y=np.linspace(h/N,h,N)
V=np.zeros(N)
Q=np.zeros(N)
Fr=np.zeros(N)

#CURVAS EN FUNCIÓN DEL CAUDAL:
A=b*y+y*y*1/np.tan(theta)
P=b+2*y/np.sin(theta)
Rh=A/P
V=np.power(Rh,2/3)*np.sqrt(S0)/n
Q=A*V
D=A/(2*(b/2+y*1/np.tan(theta)))
Fr=V/np.sqrt(g*D)
del A
del P
del Rh
del D
#PROFUNDIDAD NORMAL:
it=1
cdc=1
while cdc==1:
    P=b+2*yn1/np.sin(theta)
    A=np.power(n*Q_des*np.power(P,2/3)/np.sqrt(S0),3/5)
    yn=0.5*np.tan(theta)*(np.sqrt(b*b+4*A*1/np.tan(theta))-b)
    residual=abs(yn-yn1)
    if residual<=tol or it>=it_max:
       cdc=0
    it=it+1
    yn1=yn

#Valores de flujo normal:
An=b*yn+yn*yn*1/np.tan(theta)
Vn=Q_des/An
D=A/(2*(b/2+yn*1/np.tan(theta)))
Fr_n=Vn/np.sqrt(g*D)
print('--------------------------')
print('VALORES FLUJO NORMAL:')
print('Nivel normal yn =',np.format_float_positional(yn,3),'m')
print('Velocidad normal = ',np.format_float_positional(Vn,3),'m/s')
print('Número de Froude normal = ',np.format_float_positional(Fr_n,3))
print('Residual absoluto normal =',np.format_float_scientific(residual,3),'m')
#Dibujo del canal nivel normal:
x1=np.array([-b/2-W*np.cos(theta),-b/2,b/2,b/2+W*np.cos(theta)])
y_wall=np.array([h,0,0,h])
x2=np.array([-b/2-yn*1/np.tan(theta),b/2+yn*1/np.tan(theta)])
ysurf=np.array([yn,yn])

#PROFUNDIDAD CRÍTICA:
it=1
cdc=1
while cdc==1:
    yc=yc1-(yc1*yc1*1/np.tan(theta)+b*yc1-Q_des/np.sqrt(g*yc1))/(2*1/np.tan(theta)*yc1+b+0.5*Q_des*np.power(yc1,-1.5))/np.sqrt(g)
    residual=abs(yc-yc1)
    if residual<=tol or it>=it_max:
       cdc=0
    it=it+1
    yc1=yc
D=A/(2*(b/2+yc*1/np.tan(theta)))
Vc=np.sqrt(g*D)#[m/s]
print('--------------------------')
print('VALORES FLUJO CRÍTICO:')
print('Nivel crítico yc =',np.format_float_positional(yc,3),'m')
print('Velocidad crítica Vc [m/s] =',np.format_float_positional(Vc,3),'m/s')
print('Iteraciones totales nivel crítico =',it)
print('Residual absoluto nivel crítico =',np.format_float_scientific(residual,3),'m')
x3=np.array([-b/2-yc*1/np.tan(theta),b/2+yc*1/np.tan(theta)])
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
        A=b*yx[i]+yx[i]*yx[i]*1/np.tan(theta)
        P=b+2*yx[i]/np.sin(theta)
        Rh=A/P
        Vx[i]=Q_des/A
        D=A/(2*(b/2+yx[i]*1/np.tan(theta)))
        Frx[i]=Vx[i]/np.sqrt(g*D)
        S=n*n*Vx[i]*Vx[i]/np.power(Rh,4/3)
        yx[i+1]=yx[i]+dx*(S0-S)/(1-Frx[i]*Frx[i])
        zpiso[i+1]=zpiso[i]-S0*dx
else:
    yx[len(x)-1]=y0
    for i in range(1,len(x)-1):
        A=b*yx[len(x)-i]+1/np.tan(theta)*yx[len(x)-i]*yx[len(x)-i]
        P=b+2*yx[len(x)-i]/np.sin(theta)
        Rh=A/P
        Vx[len(x)-i]=Q_des/A
        D=A/(2*(b/2+yx[len(x)-i]*1/np.tan(theta)))
        Frx[len(x)-i]=Vx[len(x)-i]/np.sqrt(g*D)
        S=n*n*Vx[len(x)-i]*Vx[len(x)-i]/np.power(Rh,4/3)
        yx[len(x)-i-1]=yx[len(x)-i]-dx*(S0-S)/(1-Frx[len(x)-i]*Frx[len(x)-i])
        zpiso[i+1]=zpiso[i]-S0*dx
    yx[0]=yx[1]-dx*(S0-S)/(1-Frx[1]*Frx[1])
del A 
del D
A=b*yx+1/np.tan(theta)*yx*yx
Vx=Q_des/A
D=A/(2*(b/2+yx*1/np.tan(theta)))
Frx=Vx/np.sqrt(g*D)

zsurf=yx+zpiso
zn=yn+zpiso
zc=yc+zpiso
zborde=h+zpiso

#GRÁFICAS:
plt.figure(1)
plt.plot(Q,y,'b')
plt.plot(Q_des,yn,'r*')
plt.grid(True,'both','both')
plt.title('Niveles canal trapezoidal')
plt.ylabel('y [m]')
plt.xlabel('Q [m^3/s]')

plt.figure(2)
plt.plot(Q,V,'b')
plt.plot(Q_des,Vn,'r*')
plt.grid(True,'both','both')
plt.title('Velocidades canal trapezoidal')
plt.ylabel('V [m/s]')
plt.xlabel('Q [m^3/s]')

plt.figure(3)
plt.plot(Q,Fr,'b')
plt.plot(Q_des,Fr_n,'r*')
plt.grid(True,'both','both')
plt.title('Nº de Froude canal trapezoidal')
plt.ylabel('Fr')
plt.xlabel('Q [m^3/s]')

plt.figure(4)
plt.plot(x2,ysurf,'b')
plt.plot(x3,y3,'r--')
plt.plot(x1,y_wall,'k')
plt.title('Sección transversal flujo normal')
plt.ylabel('y [m]')
plt.xlabel('x [m]')
plt.xlim([-1.1*b/2-W*np.cos(theta), 1.1*b/2+W*np.cos(theta)])
plt.ylim([-0.1*h, 1.1*h])
plt.legend(['Nivel normal','Nivel crítico'])
plt.axis('equal')

plt.figure(5)
plt.subplot(311)
plt.plot(x,zborde,'k')
plt.plot(x,zsurf,'b')
plt.plot(x,zn,'g-.')
plt.plot(x,zc,'r--')
plt.plot(x,zpiso,'k')
plt.title('Nivel a lo largo del canal trapezoidal')
plt.ylabel('z [m]')
plt.legend(['Borde del canal','Superficie agua','Nivel normal','Nivel crítico','Piso canal'])

plt.subplot(312)
plt.plot(x,Vx,'b')
plt.ylabel('V [m/s]')

plt.subplot(313)
plt.plot(x,Frx,'b')
plt.xlabel('x [m] - cadecastro.com')
plt.ylabel('Fr')