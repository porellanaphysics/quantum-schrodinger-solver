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

# Parámetro de velocidad (ajústalo aquí)
INTERVALO_MS = 2  # Tiempo entre frames en milisegundos (menos = más rápido)

# Configuración de la figura
fig = plt.figure(figsize=(12, 8))
gs = GridSpec(2, 1, height_ratios=[1, 2])

# Gráfico de la norma
ax1 = fig.add_subplot(gs[0])
norm_line, = ax1.plot(np.arange(len(norm_values)), norm_values, 'g-', lw=2)
ax1.set_xlim(0, len(norm_values)-1)
ax1.set_ylim(min(norm_values)*0.999, max(norm_values)*1.001)
ax1.set_xlabel('Paso de tiempo')
ax1.set_ylabel('Norma L²')
ax1.set_title('Conservación de la norma')
ax1.grid(True)



# Gráfico de la animación
ax2 = fig.add_subplot(gs[1])
ax2.set_xlim(0, data_wave.shape[1] - 1)
wave_max = np.max(np.abs(data_wave))
ax2.set_ylim(-wave_max * 1.1, wave_max * 1.1)
ax2.set_xlabel('Posición')
ax2.set_ylabel('Amplitud')
ax2.set_title('Evolución de la función de onda y potencial')
ax2.grid(True)

# Objetos gráficos
line_wave, = ax2.plot([], [], 'b-', lw=1.5, label='Función de onda')
line_potential, = ax2.plot([], [], 'r--', lw=2, label='Potencial')
norm_point, = ax1.plot([], [], 'ro', markersize=8)
time_text = ax1.text(0.02, 0.95, '', transform=ax1.transAxes)
ax2.legend()

# Inicialización
def init():
    line_wave.set_data([], [])
    line_potential.set_data([], [])
    norm_point.set_data([], [])
    time_text.set_text('Tiempo: 0')
    return line_wave, line_potential, norm_point, time_text

# Actualización
def update(frame):
    line_wave.set_data(np.arange(data_wave.shape[1]), data_wave[frame])
    line_potential.set_data(np.arange(data_potential.shape[1]), data_potential[frame])
    norm_point.set_data([frame], [norm_values[frame]])
    time_text.set_text(f'Tiempo: {frame}/{data_wave.shape[0]-1}')
    return line_wave, line_potential, norm_point, time_text

# Crear animación
ani = animation.FuncAnimation(
    fig, update, frames=data_wave.shape[0],
    init_func=init, blit=True,
    interval=INTERVALO_MS,
    repeat=False
)

# Guardar algunos fotogramas de la animación como imágenes PNG
SAVE_FRAMES = True  # Cambiar a False si no deseas guardar los fotogramas
FRAMES_TO_SAVE = [int(data_wave.shape[0] * i / 20) for i in range(5)]  # Primeros 5 intervalos del 5%

if SAVE_FRAMES:
    for frame in FRAMES_TO_SAVE:
        update(frame)  # Actualizar los datos para el fotograma deseado
        plt.savefig(f'grafOnda/frame_{frame:04d}.png')  # Guardar el fotograma como imagen PNG
plt.tight_layout()
plt.show()

