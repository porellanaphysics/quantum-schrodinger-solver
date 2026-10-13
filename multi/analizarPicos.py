import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt



listaN = [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]
listaCiclos = [4, 6, 8, 10]

tolerancia_ruido = 0.001  # Ajusta este valor según sea necesario

for N in listaN:
    for ciclos in listaCiclos:
        archivo = f"datFiltro/filtro_N{N}_ciclos{ciclos}.csv"
        try:
            # Cargar datos del archivo
            datos = np.genfromtxt(archivo, delimiter=',', skip_header=1)
            tiempo = datos[:, 0]
            valores = datos[:, 1]
            # Encontrar picos positivos
            indices_picos, _ = find_peaks(valores, prominence=tolerancia_ruido, width=60)
            tiempos_picos = tiempo[indices_picos]
            
            # Guardar resultados
            archivo_salida = f"picos/picos_N{N}_ciclos{ciclos}.txt"
            np.savetxt(archivo_salida, tiempos_picos, header="Tiempos de picos positivos")
            print(f"Picos guardados en {archivo_salida}")
            
            # Crear y guardar la gráfica
            
            plt.figure()
            plt.plot(tiempo, valores, label="Datos")
            plt.plot(tiempo[indices_picos], valores[indices_picos], "x", label="Picos")
            plt.xlabel("Tiempo")
            plt.ylabel("Valores")
            plt.title(f"Gráfica de picos (N={N}, ciclos={ciclos})")
            plt.legend()
            archivo_grafica = f"grafPicos/grafica_N{N}_ciclos{ciclos}.png"
            plt.savefig(archivo_grafica)
            plt.close()
            print(f"Gráfica guardada en {archivo_grafica}")
        except Exception as e:
            print(f"Error procesando {archivo}: {e}")