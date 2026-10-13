import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

listaN = [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]
listaCiclos = [4, 6, 8, 10]

def leer_tiempos(file_path):
    tiempos = []
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            tiempos.append(float(line))
    return np.array(tiempos)


listaMedia = []
listaDesviacion = []

for i in listaN:
    for j in listaCiclos:
        # Ruta del archivo
        file_path = f'picos/picos_N{i}_ciclos{j}.txt'

        # Leer los tiempos y guardarlos en un array
        tiempos_array = leer_tiempos(file_path)
        print(tiempos_array)

        # Calcular la diferencia entre periodos sucesivos
        diferencias = np.diff(tiempos_array)

        # Imprimir las diferencias
        print("Diferencias entre periodos sucesivos:")
        print(diferencias)

        # Calcular la media y la desviación típica de las diferencias
        media_diferencias = np.mean(diferencias)
        desviacion_diferencias = np.std(diferencias)

        # Guardar los resultados en listas
        listaMedia.append(media_diferencias)
        listaDesviacion.append(desviacion_diferencias)

        print(f"Media de las diferencias: {media_diferencias}")
        print(f"Desviación típica de las diferencias: {desviacion_diferencias}")
        # Convertir las listas en arrays para facilitar el manejo

print(listaMedia)
# Calcular la media de las medias y la desviación típica de las medias
media_de_medias = np.mean(listaMedia)
desviacion_de_medias = np.std(listaMedia)

# Imprimir los resultados
print("----------------------------------------------------------")
print(f"Media de las medias: {media_de_medias}")
print(f"Desviación típica de las medias: {desviacion_de_medias}")
print("----------------------------------------------------------")
listaN_array = np.array(listaN)
listaCiclos_array = np.array(listaCiclos)
listaMedia_array = np.array(listaMedia).reshape(len(listaN), len(listaCiclos))
# Crear una malla para los valores de N y ciclos
N_mesh, Ciclos_mesh = np.meshgrid(listaN_array, listaCiclos_array, indexing='ij')
# Crear la figura y el eje 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# Graficar la superficie
surf = ax.plot_surface(N_mesh, Ciclos_mesh, listaMedia_array, cmap='viridis', edgecolor='k')
# Etiquetas de los ejes
ax.set_xlabel('N')
ax.set_ylabel('Ciclos')
ax.set_zlabel('Periodo (Media de las diferencias)')
# Título
ax.set_title('Periodo en función de N y Ciclos')
# Ajustar los límites del eje Z de 0 a 150
ax.set_zlim(0, 150)
# Barra de color
fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)
# Mostrar la gráfica
plt.show()
