import funcionesTP as func
import matplotlib.pyplot as plt

fig1,ax1=func.ejecutarEuler(2.2,5000,0.001)

fig2,ax2=func.ejecutarEuler_heunn(2.2,5000,0.001)

fig3,ax3=func.ejecutarRK4(2.2,5000,0.001)

plt.show()