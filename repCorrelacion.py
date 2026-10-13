import numpy as np
import matplotlib.pyplot as plt

# Leer los datos del archivo
try:
    with open('correlacion.dat', 'r') as file:
        datos = [float(line.strip()) for line in file]
except FileNotFoundError:
    print("Error: No se encontró el archivo 'correlacion.dat'")
    exit()

# Crear el gráfico
plt.figure(figsize=(10, 6))
plt.plot(datos, 'b-', linewidth=2, label='Correlación')
plt.title('Gráfico de Correlación', fontsize=16)
plt.xlabel('Índice', fontsize=14)
plt.ylabel('Valor de Correlación', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=12)

# Ajustar los límites del eje y para mejor visualización
plt.ylim(min(datos)-0.001, 1.0)

# Mostrar el gráfico
plt.tight_layout()
plt.show()