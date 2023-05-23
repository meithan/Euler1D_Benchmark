# Este programa grafica la salida del código que resuelve las ecuaciones
# de Euler en 1D, en formato ASCII.

# En el archivo de datos, en formato ASCII, cada fila debe tener 4 columnas
# separadas por espacios:
# coord x
# densidad
# velocidad
# presión

# ==============================================================================

import sys

import matplotlib.pyplot as plt
import numpy as np

# ==============================================================================

# Nombre del archivo de datos a leer (en formato ASCII)
# fname = "output_20.dat"

# Alternativamente, se puede especificar el nombre del archivo en la terminal,
# python plot_Euler1D.py output_20.dat
# descomentando la siguiente línea y comentando la anterior:
fname = sys.argv[1]

# Leemos los datos ASCII del archivo que se cargan como una matriz que
# tiene 4 columnas (posición x, densidad, velocidad, presión) y tantas
# filas como celdas en la simulación.
datos = np.loadtxt(fname)

# Para más comodido le damos definimos variables individuales para las
# columnas
x = datos[:, 0]
dens = datos[:, 1]
vel = datos[:, 2]
pres = datos[:, 3]

# También calculamos la temperatura; nótese que pres y dens son vectores
# (arreglos de numpy), así que esta operación se realiza sobre todos los
# elementos en todo el vector
temp = pres/dens

# Crear nueva figura de tamaño 10 x 8 "pulgadas"
plt.figure(figsize=(10,8))

# Con subplot(filas, columnas, posición) especificamos que queremos un arreglo
# de sub-gráficas con ese número de filas y columnas, y que actualmente queremos
# graficar en esa posición (la numeración empieza en 1 y corre de izquierda a
# derecha, y de arriba hacia abajo)
plt.subplot(2,2,1)
# El comando plot() recibe un vector con los valores de x y otro con los valores
# de y como los primeros dos argumentos. Los demás argumentos son opciones,
# como el estilo de los marcadores y el color.
plt.plot(x, dens, marker=".", color="blue")
# Título de esta (sub-)gráfica
plt.title("Densidad")
# Rango del eje x que queremos mostrar
plt.xlim(0, 1)
# Rango del eje y que queremos mostrar
plt.ylim(0, 1.05)
# Cuadrícula
plt.grid(ls=":")

# Los siguientes 3 grupos grafican diferentes variables en subplots distintos

plt.subplot(2,2,2)
plt.plot(x, vel, marker=".", color="green")
plt.title("Velocidad")
plt.xlim(0, 1)
plt.ylim(-0.05, 1.0)
plt.grid(ls=":")

plt.subplot(2,2,3)
plt.plot(x, pres, marker=".", color="darkorange")
plt.title("Presión")
plt.xlim(0, 1)
plt.ylim(0, 1.05)
plt.grid(ls=":")

plt.subplot(2,2,4)
plt.plot(x, temp, marker=".", color="red")
plt.title("Temperatura")
plt.xlim(0, 1)
plt.ylim(0.7, 1.2)
plt.grid(ls=":")

plt.tight_layout()

if "--save" in sys.argv:
  out_fname = "output.png"
  plt.savefig(out_fname)
  print("Wrote", out_fname)
else:
  plt.show()
