from cmath import sqrt
from math import pi
from math import cos
from math import sin
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#masas
#----------------------------------------------------------------
ca=pi/180

#tierra
vp=30.29e3 #velocidad perihelio
pp=147.10e9 #posicion perihelio
a=282.94*ca #angulo
i=0.00*ca #inclinacion

#cuerpo5 Marte
vp2=26.614e3 #velocidad perihelio
pp2=206.67e9 #posicion perihelio 
a2=286.46*ca 
i2=1.85*ca


#----------------------------------------------------------------
m3 = 1.989e30  #masa del sol
m1 = 5.972e24  #masa de la tierra
m2 = 6.39e23 #masa de marte
def init(pasos):
    #G = 6.67430e-11  
    # m_sol = 1.989e30  
    # m_tierra = 5.972e24  
    # m_planetaX = 6.39e23

    et=np.zeros(pasos)

    p1=np.zeros((pasos,3))
    p2=np.zeros((pasos,3))
    p3=np.zeros((pasos,3))

    #estable
    # p1[0],p2[0],p3[0]=np.array([1.496e11,0.0,0.0]),np.array([2.279e11,0.0,0.0]),np.array([0.0,0.0,0.0])
    #inestable 1
    # p1[0],p2[0],p3[0]=np.array([2.2e11,0.0,0.0]),np.array([-2.2e11,0.0,0.0]),np.array([0.0,0.0,0.0])
    #paper
    p1[0],p2[0],p3[0]=np.array([pp*cos(a)*cos(i),pp*sin(a)*cos(i),pp*sin(i)]),np.array([pp2*cos(a2)*cos(i2),pp2*sin(a2)*cos(i2),pp2*sin(i2)]),np.array([0.0,0.0,0.0])

    v1=np.zeros((pasos,3))
    v2=np.zeros((pasos,3))
    v3=np.zeros((pasos,3))

    #estable
    # v1[0],v2[0],v3[0]=np.array([0.0,29780.0,0.0]),np.array([0.0,24077.0,0.0]),np.array([0.0,0.0,0.0])
    #inestable 1
    # v1[0],v2[0],v3[0]=np.array([0.0,25000.0,0.0]),np.array([0.0,25000.0,0.0]),np.array([0.0,0.0,5000])
    #paper
    v1[0],v2[0],v3[0]=np.array([-vp*sin(a),vp*cos(a),0.0]),np.array([-vp2*sin(a2),vp2*cos(a2),0.0]),np.array([0.0,0.0,0.0])
    return p1,p2,p3,v1,v2,v3,et

def calcularAceleracion(m2,m3,r1,r2,r3,G):
    r1_2=r2-r1
    r1_3=r3-r1
    distancia1=np.linalg.norm(r1_2)
    distancia2=np.linalg.norm(r1_3)
    aceleracion1=(G*m2*r1_2)/(distancia1**3)
    aceleracion2=(G*m3*r1_3)/(distancia2**3)
    return aceleracion1+aceleracion2

# def graficar(p1,p2,p3,titulo_ventana='Gráfico 3D'):
#     figura = plt.figure()
#     ax = figura.add_subplot(111, projection='3d')
#     ax.plot(p1[:,0],p1[:,1],p1[:,2],color='red',label='Cuerpo 1')
#     ax.plot(p2[:,0],p2[:,1],p2[:,2],color='blue',label='Cuerpo 2')
#     ax.plot(p3[:,0],p3[:,1],p3[:,2],color='green',label='Cuerpo 3')
#     ax.scatter(p1[-1,0],p1[-1,1],p1[-1,2],color='red',s=30,label='Fin Cuerpo 1')
#     ax.scatter(p2[-1,0],p2[-1,1],p2[-1,2],color='blue',s=30,label='Fin Cuerpo 2')
#     ax.scatter(p3[-1,0],p3[-1,1],p3[-1,2],color='green',s=30,label='Fin Cuerpo 3')
#     # Establecer el título de la ventana
#     figura.canvas.manager.set_window_title(titulo_ventana)
#     # ax.legend()
#     # plt.show()
#     return figura,ax
def graficar(p1,p2,p3,energia_total, titulo_ventana='Gráfico 3D'):
    figura = plt.figure(figsize=(12, 6))

    # Subplot 1: Gráfico 3D de las posiciones
    ax1 = figura.add_subplot(121, projection='3d')
    ax1.plot(p1[:,0], p1[:,1], p1[:,2], color='red', label='Cuerpo 1')
    ax1.plot(p2[:,0], p2[:,1], p2[:,2], color='blue', label='Cuerpo 2')
    ax1.plot(p3[:,0], p3[:,1], p3[:,2], color='green', label='Cuerpo 3')
    ax1.scatter(p1[-1,0], p1[-1,1], p1[-1,2], color='red', s=30, label='Fin Cuerpo 1')
    ax1.scatter(p2[-1,0], p2[-1,1], p2[-1,2], color='blue', s=30, label='Fin Cuerpo 2')
    ax1.scatter(p3[-1,0], p3[-1,1], p3[-1,2], color='green', s=30, label='Fin Cuerpo 3')
    ax1.set_title('Trayectorias de los cuerpos')
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('Z')
    ax1.legend()

    # Subplot 2: Gráfico de la energía total
    ax2 = figura.add_subplot(122)
    ax2.plot(energia_total, color='purple', label='Energía Total')
    ax2.set_title('Energía Total del Sistema')
    ax2.set_xlabel('Tiempo')
    ax2.set_ylabel('Energía')
    ax2.legend()
    # Establecer el título de la ventana
    figura.canvas.manager.set_window_title(titulo_ventana)

    return figura, (ax1, ax2)

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

    p1,p2,p3,v1,v2,v3,et=init(pasos)

    for i in range(pasos-1):
        et[i]=calcularEnergia_i(p1[i],p2[i],p3[i],v1[i],v2[i],v3[i],G)
        p1[i+1]=calcularEuler(p1[i],v1[i],delta_t)
        p2[i+1]=calcularEuler(p2[i],v2[i],delta_t)
        p3[i+1]=calcularEuler(p3[i],v3[i],delta_t)

        v1[i+1]=calcularEuler(v1[i],calcularAceleracion(m2,m3,p1[i],p2[i],p3[i],G),delta_t)
        v2[i+1]=calcularEuler(v2[i],calcularAceleracion(m1,m3,p2[i],p1[i],p3[i],G),delta_t)
        v3[i+1]=calcularEuler(v3[i],calcularAceleracion(m1,m2,p3[i],p1[i],p2[i],G),delta_t)
    # Calcular la energía total para el último paso
    et[pasos - 1] = calcularEnergia_i(p1[pasos - 1], p2[pasos - 1], p3[pasos - 1], v1[pasos - 1], v2[pasos - 1], v3[pasos - 1], G)

    print('------------------EULER-----------------')
    print('energia segun trapecio', ejecutar_trapecio(et, delta_t))
    print('energia segun simpson', ejecutar_simpson(et, delta_t))
    print('energia segun gauss', ejecutar_cuadratura_gauss(et, delta_t))
    # print("et",et)
    # et_filtrado = [e for e in et if e != 0]
    
    return graficar(p1,p2,p3,et,titulo_ventana)

def ejecutarEuler_heunn(G,pasos,delta_t):
    titulo_ventana='Euler-Heunn'
    p1,p2,p3,v1,v2,v3,et=init(pasos)

    for i in range(pasos-1):
        et[i]=calcularEnergia_i(p1[i],p2[i],p3[i],v1[i],v2[i],v3[i],G)
        p1_2=calcularEuler(p1[i],v1[i],delta_t)
        p2_2=calcularEuler(p2[i],v2[i],delta_t)
        p3_2=calcularEuler(p3[i],v3[i],delta_t)

        v1[i+1]=calcularEuler_Heunn(v1[i],calcularAceleracion(m2,m3,p1[i],p2[i],p3[i],G),calcularAceleracion(m2,m3,p1_2,p2_2,p3_2,G),delta_t)
        v2[i+1]=calcularEuler_Heunn(v2[i],calcularAceleracion(m1,m3,p2[i],p1[i],p3[i],G),calcularAceleracion(m1,m3,p2_2,p1_2,p3_2,G),delta_t)
        v3[i+1]=calcularEuler_Heunn(v3[i],calcularAceleracion(m1,m2,p3[i],p1[i],p2[i],G),calcularAceleracion(m1,m2,p3_2,p1_2,p2_2,G),delta_t)

        p1[i+1]=calcularEuler_Heunn(p1[i],v1[i],v1[i+1],delta_t)
        p2[i+1]=calcularEuler_Heunn(p2[i],v2[i],v2[i+1],delta_t)
        p3[i+1]=calcularEuler_Heunn(p3[i],v3[i],v3[i+1],delta_t)

    et[pasos - 1] = calcularEnergia_i(p1[pasos - 1], p2[pasos - 1], p3[pasos - 1], v1[pasos - 1], v2[pasos - 1], v3[pasos - 1], G)
    print('------------------EULER-HEUN-----------------')
    print('energia segun trapecio', ejecutar_trapecio(et, delta_t))
    print('energia segun simpson', ejecutar_simpson(et, delta_t))
    print('energia segun gauss', ejecutar_cuadratura_gauss(et, delta_t))
    return graficar(p1,p2,p3,et,titulo_ventana)
#Fin__1-A-EULER----------------------------------------------------------------------------------------------------

# def ejecutarRK4(G,pasos,delta_t):
#     titulo_ventana='RK4'
#     p1,p2,p3,v1,v2,v3,et=init(pasos)
    
#     for i in range(pasos-1):
#         et[i]=calcularEnergia_i(p1[i],p2[i],p3[i],v1[i],v2[i],v3[i],G)

#         a1=calcularAceleracion(m2,m3,p1[i],p2[i],p3[i],G)
#         a2=calcularAceleracion(m1,m3,p2[i],p1[i],p3[i],G)
#         a3=calcularAceleracion(m1,m2,p3[i],p1[i],p2[i],G)
#         if(i>0):
#             v1[i]=v1[i-1]+a1*delta_t
#             v2[i]=v2[i-1]+a2*delta_t
#             v3[i]=v3[i-1]+a3*delta_t

#         k1_1=v1[i]
#         k1_2=v2[i]
#         k1_3=v3[i]

#         k2_1=v1[i]+a1*delta_t*k1_1/2
#         k2_2=v2[i]+a2*delta_t*k1_2/2
#         k2_3=v3[i]+a3*delta_t*k1_3/2

#         k3_1=v1[i]+a1*delta_t*k2_1/2
#         k3_2=v2[i]+a2*delta_t*k2_2/2
#         k3_3=v3[i]+a3*delta_t*k2_3/2

#         k4_1=v1[i]+a1*delta_t*k3_1
#         k4_2=v2[i]+a2*delta_t*k3_2
#         k4_3=v3[i]+a3*delta_t*k3_3

#         p1[i+1]=p1[i]+(delta_t/6)*(k1_1+2*k2_1+2*k3_1+k4_1)
#         p2[i+1]=p2[i]+(delta_t/6)*(k1_2+2*k2_2+2*k3_2+k4_2)
#         p3[i+1]=p3[i]+(delta_t/6)*(k1_3+2*k2_3+2*k3_3+k4_3)

#     return graficar(p1,p2,p3,et,titulo_ventana)

def ejecutarRK4(G, pasos, delta_t):
    titulo_ventana = 'RK4'
    p1, p2, p3, v1, v2, v3, et = init(pasos)
    
    for i in range(pasos - 1):
        et[i] = calcularEnergia_i(p1[i], p2[i], p3[i], v1[i], v2[i], v3[i], G)

        # Calcular aceleraciones iniciales
        a1 = calcularAceleracion(m2, m3, p1[i], p2[i], p3[i], G)
        a2 = calcularAceleracion(m1, m3, p2[i], p1[i], p3[i], G)
        a3 = calcularAceleracion(m1, m2, p3[i], p1[i], p2[i], G)

        # k1
        k1_v1 = a1
        k1_v2 = a2
        k1_v3 = a3
        k1_p1 = v1[i]
        k1_p2 = v2[i]
        k1_p3 = v3[i]

        # k2
        k2_v1 = calcularAceleracion(m2, m3, p1[i] + k1_p1 * delta_t / 2, p2[i] + k1_p2 * delta_t / 2, p3[i] + k1_p3 * delta_t / 2, G)
        k2_v2 = calcularAceleracion(m1, m3, p2[i] + k1_p2 * delta_t / 2, p1[i] + k1_p1 * delta_t / 2, p3[i] + k1_p3 * delta_t / 2, G)
        k2_v3 = calcularAceleracion(m1, m2, p3[i] + k1_p3 * delta_t / 2, p1[i] + k1_p1 * delta_t / 2, p2[i] + k1_p2 * delta_t / 2, G)
        k2_p1 = v1[i] + k1_v1 * delta_t / 2
        k2_p2 = v2[i] + k1_v2 * delta_t / 2
        k2_p3 = v3[i] + k1_v3 * delta_t / 2

        # k3
        k3_v1 = calcularAceleracion(m2, m3, p1[i] + k2_p1 * delta_t / 2, p2[i] + k2_p2 * delta_t / 2, p3[i] + k2_p3 * delta_t / 2, G)
        k3_v2 = calcularAceleracion(m1, m3, p2[i] + k2_p2 * delta_t / 2, p1[i] + k2_p1 * delta_t / 2, p3[i] + k2_p3 * delta_t / 2, G)
        k3_v3 = calcularAceleracion(m1, m2, p3[i] + k2_p3 * delta_t / 2, p1[i] + k2_p1 * delta_t / 2, p2[i] + k2_p2 * delta_t / 2, G)
        k3_p1 = v1[i] + k2_v1 * delta_t / 2
        k3_p2 = v2[i] + k2_v2 * delta_t / 2
        k3_p3 = v3[i] + k2_v3 * delta_t / 2

        # k4
        k4_v1 = calcularAceleracion(m2, m3, p1[i] + k3_p1 * delta_t, p2[i] + k3_p2 * delta_t, p3[i] + k3_p3 * delta_t, G)
        k4_v2 = calcularAceleracion(m1, m3, p2[i] + k3_p2 * delta_t, p1[i] + k3_p1 * delta_t, p3[i] + k3_p3 * delta_t, G)
        k4_v3 = calcularAceleracion(m1, m2, p3[i] + k3_p3 * delta_t, p1[i] + k3_p1 * delta_t, p2[i] + k3_p2 * delta_t, G)
        k4_p1 = v1[i] + k3_v1 * delta_t
        k4_p2 = v2[i] + k3_v2 * delta_t
        k4_p3 = v3[i] + k3_v3 * delta_t

        # Actualizar posiciones y velocidades
        p1[i + 1] = p1[i] + (delta_t / 6) * (k1_p1 + 2 * k2_p1 + 2 * k3_p1 + k4_p1)
        p2[i + 1] = p2[i] + (delta_t / 6) * (k1_p2 + 2 * k2_p2 + 2 * k3_p2 + k4_p2)
        p3[i + 1] = p3[i] + (delta_t / 6) * (k1_p3 + 2 * k2_p3 + 2 * k3_p3 + k4_p3)
        v1[i + 1] = v1[i] + (delta_t / 6) * (k1_v1 + 2 * k2_v1 + 2 * k3_v1 + k4_v1)
        v2[i + 1] = v2[i] + (delta_t / 6) * (k1_v2 + 2 * k2_v2 + 2 * k3_v2 + k4_v2)
        v3[i + 1] = v3[i] + (delta_t / 6) * (k1_v3 + 2 * k2_v3 + 2 * k3_v3 + k4_v3)

    et[pasos - 1] = calcularEnergia_i(p1[pasos - 1], p2[pasos - 1], p3[pasos - 1], v1[pasos - 1], v2[pasos - 1], v3[pasos - 1], G)

    print('------------------RK4-----------------')
    print('energia segun trapecio', ejecutar_trapecio(et, delta_t))
    print('energia segun simpson', ejecutar_simpson(et, delta_t))
    print('energia segun gauss', ejecutar_cuadratura_gauss(et, delta_t))
 
    return graficar(p1, p2, p3, et, titulo_ventana)
    
def ejecutarRK5(G, pasos, delta_t):
    titulo_ventana = 'RK5'
    p1, p2, p3, v1, v2, v3, et = init(pasos)
    
    for i in range(pasos - 1):
        et[i] = calcularEnergia_i(p1[i], p2[i], p3[i], v1[i], v2[i], v3[i], G)
        # k1
        k1_1 = v1[i]
        k1_2 = v2[i]
        k1_3 = v3[i]

        # p1+1
        p1_2 = p1[i] + (delta_t / 4) * k1_1
        p2_2 = p2[i] + (delta_t / 4) * k1_2
        p3_2 = p3[i] + (delta_t / 4) * k1_3

        # k2
        a2_1 = calcularAceleracion(m2, m3, p1_2, p2_2, p3_2, G)
        a2_2 = calcularAceleracion(m1, m3, p2_2, p1_2, p3_2, G)
        a2_3 = calcularAceleracion(m1, m2, p3_2, p1_2, p2_2, G)
        k2_1 = v1[i] + (delta_t / 4) * a2_1
        k2_2 = v2[i] + (delta_t / 4) * a2_2
        k2_3 = v3[i] + (delta_t / 4) * a2_3

        # p2+1
        p1_3 = p1[i] + (delta_t / 8) * (k1_1 + k2_1)
        p2_3 = p2[i] + (delta_t / 8) * (k1_2 + k2_2)
        p3_3 = p3[i] + (delta_t / 8) * (k1_3 + k2_3)

        # k3
        a3_1 = calcularAceleracion(m2, m3, p1_3, p2_3, p3_3, G)
        a3_2 = calcularAceleracion(m1, m3, p2_3, p1_3, p3_3, G)
        a3_3 = calcularAceleracion(m1, m2, p3_3, p1_3, p2_3, G)
        k3_1 = v1[i] + (delta_t / 4) * a3_1
        k3_2 = v2[i] + (delta_t / 4) * a3_2
        k3_3 = v3[i] + (delta_t / 4) * a3_3

        # p3+1
        p1_4 = p1[i] + delta_t * (k3_1 - (k2_1 / 2))
        p2_4 = p2[i] + delta_t * (k3_2 - (k2_2 / 2))
        p3_4 = p3[i] + delta_t * (k3_3 - (k2_3 / 2))

        # k4
        a4_1 = calcularAceleracion(m2, m3, p1_4, p2_4, p3_4, G)
        a4_2 = calcularAceleracion(m1, m3, p2_4, p1_4, p3_4, G)
        a4_3 = calcularAceleracion(m1, m2, p3_4, p1_4, p2_4, G)
        k4_1 = v1[i] + (delta_t / 2) * a4_1
        k4_2 = v2[i] + (delta_t / 2) * a4_2
        k4_3 = v3[i] + (delta_t / 2) * a4_3

        # p4+1
        p1_5 = p1[i] + (delta_t) * (k1_1 * 3 / 16 + k4_1 * 9 / 16)
        p2_5 = p2[i] + (delta_t) * (k1_2 * 3 / 16 + k4_2 * 9 / 16)
        p3_5 = p3[i] + (delta_t) * (k1_3 * 3 / 16 + k4_3 * 9 / 16)

        # k5
        a5_1 = calcularAceleracion(m2, m3, p1_5, p2_5, p3_5, G)
        a5_2 = calcularAceleracion(m1, m3, p2_5, p1_5, p3_5, G)
        a5_3 = calcularAceleracion(m1, m2, p3_5, p1_5, p2_5, G)
        k5_1 = v1[i] + (delta_t * 3 / 4) * a5_1
        k5_2 = v2[i] + (delta_t * 3 / 4) * a5_2
        k5_3 = v3[i] + (delta_t * 3 / 4) * a5_3

        # p5+1
        p1_6 = p1[i] + (delta_t / 7) * (k1_1 * (-3) + 2 * k2_1 + 12 * k3_1 - k4_1 * 12 + k5_1 * 8)
        p2_6 = p2[i] + (delta_t / 7) * (k1_2 * (-3) + 2 * k2_2 + 12 * k3_2 - k4_2 * 12 + k5_2 * 8)
        p3_6 = p3[i] + (delta_t / 7) * (k1_3 * (-3) + 2 * k2_3 + 12 * k3_3 - k4_3 * 12 + k5_3 * 8)

        # k6
        a6_1 = calcularAceleracion(m2, m3, p1_6, p2_6, p3_6, G)
        a6_2 = calcularAceleracion(m1, m3, p2_6, p1_6, p3_6, G)
        a6_3 = calcularAceleracion(m1, m2, p3_6, p1_6, p2_6, G)
        k6_1 = v1[i] + delta_t * a6_1
        k6_2 = v2[i] + delta_t * a6_2
        k6_3 = v3[i] + delta_t * a6_3

        # p+1
        p1[i + 1] = p1[i] + (delta_t / 90) * (7 * k1_1 + 32 * k3_1 + 12 * k4_1 + 32 * k5_1 + 7 * k6_1)
        p2[i + 1] = p2[i] + (delta_t / 90) * (7 * k1_2 + 32 * k3_2 + 12 * k4_2 + 32 * k5_2 + 7 * k6_2)
        p3[i + 1] = p3[i] + (delta_t / 90) * (7 * k1_3 + 32 * k3_3 + 12 * k4_3 + 32 * k5_3 + 7 * k6_3)
        
        v1[i + 1] = v1[i] + (delta_t / 90) * (7 * a2_1 + 32 * a3_1 + 12 * a4_1 + 32 * a5_1 + 7 * a6_1)
        v2[i + 1] = v2[i] + (delta_t / 90) * (7 * a2_2 + 32 * a3_2 + 12 * a4_2 + 32 * a5_2 + 7 * a6_2)
        v3[i + 1] = v3[i] + (delta_t / 90) * (7 * a2_3 + 32 * a3_3 + 12 * a4_3 + 32 * a5_3 + 7 * a6_3)

    et[pasos - 1] = calcularEnergia_i(p1[pasos - 1], p2[pasos - 1], p3[pasos - 1], v1[pasos - 1], v2[pasos - 1], v3[pasos - 1], G)

    print('------------------RK5-----------------')
    print('energia segun trapecio', ejecutar_trapecio(et, delta_t))
    print('energia segun simpson', ejecutar_simpson(et, delta_t))
    print('energia segun gauss', ejecutar_cuadratura_gauss(et, delta_t))

    return graficar(p1, p2, p3, et, titulo_ventana)



def calcularECinetica(m_i,v_i):
    return (m_i/2)*np.linalg.norm(v_i)**2

def calcularEPotencial(m1,m2,m3,p1,p2,p3,G):
    r1_2=p2-p1
    r1_3=p3-p1
    r2_3 = p3 - p2
    distancia1=np.linalg.norm(r1_2)
    distancia2=np.linalg.norm(r1_3)
    distancia3 = np.linalg.norm(r2_3)
    return -G * m1 * m2 / distancia1 - G * m1 * m3 / distancia2 - G * m2 * m3 / distancia3 #suma de las energias potenciales de los cuerpos

def calcularEnergiaTotal(ep,ec):
    return ep+ec

def calcularEnergia_i(p1_i,p2_i,p3_i,v1_i,v2_i,v3_i,G):
    ep=calcularEPotencial(m1,m2,m3,p1_i,p2_i,p3_i,G)
    ec=calcularECinetica(m1,v1_i)+calcularECinetica(m2,v2_i)+calcularECinetica(m3,v3_i)
    et=calcularEnergiaTotal(ep,ec)
    return np.abs(et)

def ejecutar_trapecio(et, delta_t):
    acumulado = 0

    for i in range(1, len(et)):
        acumulado += abs((et[i] + et[i-1]) * delta_t / 2)

    return acumulado

def ejecutar_simpson(et, delta_t):
    acumulado = 0
  
    for i in range(1, len(et) - 1):
        acumulado += abs((et[i-1] + 4 * et[i] + et[i+1]) * delta_t / 6)

    return acumulado

def ejecutar_cuadratura_gauss(et, delta_t):
    acumulado = 0

    for i in range(1, len(et) - 1):
        acumulado += abs((et[i-1] + 4 * et[i] + et[i+1]) * delta_t / 6)

    return acumulado