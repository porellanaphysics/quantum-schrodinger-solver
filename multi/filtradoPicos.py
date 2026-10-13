import numpy as np
import matplotlib.pyplot as plt

listaN = [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]
listaCiclos = [4, 6, 8, 10]

for i in listaN:
    for j in listaCiclos:
        
        # Leer datos del archivo .dat
        datos = np.loadtxt(f'resultados_N{i}_ciclos{j}.dat')
        x = datos[:, 0]
        y = datos[:, 1]

        # Parámetros del filtro
        ventana = 41  # Tamaño de la ventana (ajusta según necesites, mejor si es impar)

        # Aplicar media móvil sin bordes
        y_filtrado = np.convolve(y, np.ones(ventana)/ventana, mode='valid')

        # Ajustar el eje x para coincidir con los datos filtrados
        x_filtrado = x[(ventana-1)//2 : -(ventana-1)//2]

        # Crear array con los resultados
        resultados = np.column_stack((x_filtrado, y_filtrado))

        # Guardar en archivo CSV
        np.savetxt(f'datFiltro/filtro_N{i}_ciclos{j}.csv', resultados, 
                   header='tiempo,valor_filtrado', 
                   delimiter=',', 
                   fmt='%.8f', 
                   comments='')

        # Graficar
        plt.figure()
        plt.plot(x, y, label='Original', alpha=0.5)
        plt.plot(x_filtrado, y_filtrado, label=f'Media móvil', linewidth=2, color='red')
        plt.title('Resultado del Filtrado y Datos Guardados en CSV')
        plt.legend()
        plt.grid(True)
        #plt.show()
        plt.savefig(f'graficas/grafico_filtro_N{i}_ciclos{j}.png')

        print(f"¡Datos guardados exitosamente en 'datos_filtrados.csv'!")
        print(f"Se filtraron {ventana-1} puntos ({(ventana-1)/len(x)*100:.1f}% de los datos)")