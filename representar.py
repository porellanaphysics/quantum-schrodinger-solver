import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.gridspec import GridSpec

# Cargar todos los datos
data_wave = np.loadtxt('onda_resultado.dat')  # Función de onda
data_potential = np.loadtxt('potencial.dat')  # Potencial
norm_values = np.loadtxt('norma.dat')  # Valores de la norma

# Verificar y preparar los datos del potencial
if data_potential.shape != data_wave.shape:
    if len(data_potential.shape) == 1:
        data_potential = np.tile(data_potential, (data_wave.shape[0], 1))
    else:
        raise ValueError("Dimensiones incompatibles entre onda y potencial")

# Calcular amplitud máxima inicial de la onda (primer frame)
initial_wave_max = np.max(np.abs(data_wave[0]))

# Escalar el potencial para visualización (1.5 veces la amplitud máxima inicial)
scaling_factor = 1.5 * initial_wave_max / np.max(data_potential) if np.max(data_potential) > 0 else 1
data_potential_display = data_potential * scaling_factor

# Crear el array de tiempos
times = np.arange(len(norm_values))

# Configuración de la figura
fig = plt.figure(figsize=(12, 8))
gs = GridSpec(2, 1, height_ratios=[1, 2])

# Gráfico de la norma
ax1 = fig.add_subplot(gs[0])
norm_line, = ax1.plot(times, norm_values, 'g-', lw=2)
ax1.set_xlabel('Paso de tiempo')
ax1.set_ylabel('Norma L²')
ax1.set_title('Conservación de la norma')
ax1.grid(True)
ax1.set_ylim(min(norm_values)*0.999, max(norm_values)*1.001)

# Gráfico de la animación
ax2 = fig.add_subplot(gs[1])
ax2.set_xlim(0, data_wave.shape[1] - 1)

# Ajustar límites Y considerando tanto la onda como el potencial escalado
combined_max = max(np.max(np.abs(data_wave)), np.max(data_potential_display))
ax2.set_ylim(-combined_max * 1.1, combined_max * 1.1)
ax2.set_xlabel('Posición')
ax2.set_ylabel('Amplitud')
ax2.set_title(f'Evolución de la función de onda y potencial (escalado {scaling_factor:.2f}x)')
ax2.grid(True)

# Objetos gráficos
line_wave, = ax2.plot([], [], 'b-', lw=1.5, label='Función de onda')
line_potential, = ax2.plot([], [], 'r--', lw=2, label='Potencial (escalado)')
ax2.legend()

# Punto móvil para la norma
norm_point, = ax1.plot([], [], 'ro', markersize=8)

# Inicialización
def init():
    line_wave.set_data([], [])
    line_potential.set_data([], [])
    norm_point.set_data([], [])
    return line_wave, line_potential, norm_point,

# Actualización
def update(frame):
    x = np.arange(data_wave.shape[1])
    line_wave.set_data(x, data_wave[frame])
    line_potential.set_data(x, data_potential_display[frame])
    norm_point.set_data([times[frame]], [norm_values[frame]])
    return line_wave, line_potential, norm_point,

# Crear animación
ani = animation.FuncAnimation(
    fig, update, 
    frames=min(len(times), data_wave.shape[0]),
    init_func=init, 
    blit=True, 
    interval=50
)

plt.tight_layout()
plt.show()