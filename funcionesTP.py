from cmath import sqrt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#masas
m1=2
m2=1
m3=1

def init(pasos):
    p1=np.zeros((pasos,3))
    p2=np.zeros((pasos,3))
    p3=np.zeros((pasos,3))

    p1[0],p2[0],p3[0]=np.array([1.0,1.0,1.0]),np.array([-1,-1,1]),np.array([-1,1,-1])

    v1=np.zeros((pasos,3))
    v2=np.zeros((pasos,3))
    v3=np.zeros((pasos,3))

    v1[0],v2[0],v3[0]=np.array([0.5,-0.5,0.0]),np.array([-0.5,0.0,0.5]),np.array([0.0,0.5,-0.5])
    return p1,p2,p3,v1,v2,v3

def calcularAceleracion(m2,m3,r1,r2,r3,G):
    r1_2=r2-r1
    r1_3=r3-r1
    distancia1=np.linalg.norm(r1_2)
    distancia2=np.linalg.norm(r1_3)
    aceleracion1=(G*m2*r1_2)/(distancia1**3)
    aceleracion2=(G*m3*r1_3)/(distancia2**3)
    return aceleracion1+aceleracion2

def graficar(p1,p2,p3,titulo_ventana='Gráfico 3D'):
    figura = plt.figure()
    ax = figura.add_subplot(111, projection='3d')
    ax.plot(p1[:,0],p1[:,1],p1[:,2],color='red',label='Cuerpo 1')
    ax.plot(p2[:,0],p2[:,1],p2[:,2],color='blue',label='Cuerpo 2')
    ax.plot(p3[:,0],p3[:,1],p3[:,2],color='green',label='Cuerpo 3')
    ax.scatter(p1[-1,0],p1[-1,1],p1[-1,2],color='red',s=30,label='Fin Cuerpo 1')
    ax.scatter(p2[-1,0],p2[-1,1],p2[-1,2],color='blue',s=30,label='Fin Cuerpo 2')
    ax.scatter(p3[-1,0],p3[-1,1],p3[-1,2],color='green',s=30,label='Fin Cuerpo 3')
    # Establecer el título de la ventana
    figura.canvas.manager.set_window_title(titulo_ventana)
    # ax.legend()
    # plt.show()
    return figura,ax

#1-A-EULER----------------------------------------------------------------------------------------------------
def calcularEuler(c1, c2, Delta_t):
    y=c1+c2*Delta_t
    return y
def calcularEuler_Heunn(x_1,v_1,v_2,Delta_t):
    x_2=x_1+(v_1+v_2)*Delta_t/2
    return x_2

def ejecutarEuler(G,pasos,delta_t):
    # G=2
    # pasos=5000
    # delta_t=0.001
    titulo_ventana='Euler'

    p1,p2,p3,v1,v2,v3=init(pasos)

    for i in range(pasos-1):
        
        p1[i+1]=calcularEuler(p1[i],v1[i],delta_t)
        p2[i+1]=calcularEuler(p2[i],v2[i],delta_t)
        p3[i+1]=calcularEuler(p3[i],v3[i],delta_t)

        v1[i+1]=calcularEuler(v1[i],calcularAceleracion(m2,m3,p1[i],p2[i],p3[i],G),delta_t)
        v2[i+1]=calcularEuler(v2[i],calcularAceleracion(m1,m3,p2[i],p1[i],p3[i],G),delta_t)
        v3[i+1]=calcularEuler(v3[i],calcularAceleracion(m1,m2,p3[i],p1[i],p2[i],G),delta_t)

    return graficar(p1,p2,p3,titulo_ventana)

def ejecutarEuler_heunn(G,pasos,delta_t):
    titulo_ventana='Euler-Heunn'
    p1,p2,p3,v1,v2,v3=init(pasos)

    for i in range(pasos-1):
        
        p1_2=calcularEuler(p1[i],v1[i],delta_t)
        p2_2=calcularEuler(p2[i],v2[i],delta_t)
        p3_2=calcularEuler(p3[i],v3[i],delta_t)

        v1[i+1]=calcularEuler_Heunn(v1[i],calcularAceleracion(m2,m3,p1[i],p2[i],p3[i],G),calcularAceleracion(m2,m3,p1_2,p2_2,p3_2,G),delta_t)
        v2[i+1]=calcularEuler_Heunn(v2[i],calcularAceleracion(m1,m3,p2[i],p1[i],p3[i],G),calcularAceleracion(m1,m3,p2_2,p1_2,p3_2,G),delta_t)
        v3[i+1]=calcularEuler_Heunn(v3[i],calcularAceleracion(m1,m2,p3[i],p1[i],p2[i],G),calcularAceleracion(m1,m2,p3_2,p1_2,p2_2,G),delta_t)

        p1[i+1]=calcularEuler_Heunn(p1[i],v1[i],v1[i+1],delta_t)
        p2[i+1]=calcularEuler_Heunn(p2[i],v2[i],v2[i+1],delta_t)
        p3[i+1]=calcularEuler_Heunn(p3[i],v3[i],v3[i+1],delta_t)

    return graficar(p1,p2,p3,titulo_ventana)
#Fin__1-A-EULER----------------------------------------------------------------------------------------------------