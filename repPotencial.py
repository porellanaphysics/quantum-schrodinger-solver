import matplotlib.pyplot as plt

# Leer los datos del archivo
with open('potencial.txt', 'r') as file:
    data = [int(line.strip()) for line in file if line.strip()]

# Crear la gráfica
plt.figure(figsize=(10, 6))
plt.plot(data, 'b-', linewidth=1.5)  # Línea azul continua
plt.title('Gráfica del Potencial', fontsize=14)
plt.xlabel('Índice', fontsize=12)
plt.ylabel('Valor del Potencial', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)  # Cuadrícula punteada
plt.tight_layout()

# Mostrar la gráfica
plt.show()