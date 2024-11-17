import funcionesTP as func
import matplotlib.pyplot as plt

fig1,(ax1,ax11)=func.ejecutarEuler(6.67430e-11 ,31536,1000)

fig2,(ax2,ax22)=func.ejecutarEuler_heunn(6.67430e-11 ,31536,1000)

fig3,(ax3,ax33)=func.ejecutarRK4(6.67430e-11 ,31536,1000)

fig3,(ax3,ax33)=func.ejecutarRK5(6.67430e-11 ,31536,1000)

plt.show()