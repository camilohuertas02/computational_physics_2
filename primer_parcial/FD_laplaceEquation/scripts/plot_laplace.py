# scripts/plot_laplace.py

# Script para graficar la solución de la ecuación de Laplace con Matplotlib


import matplotlib

matplotlib.use('TkAgg')

import matplotlib.pyplot as plt

import numpy as np



# 📥 Leer datos desde el archivo generado por C++

x, y, V = [], [], []


with open("laplace.dat") as file:

    for line in file:

        if line.strip() == "" or line.startswith("#"):

            continue

        xi, yi, Vi = map(float, line.strip().split())

        x.append(xi)

        y.append(yi)

        V.append(Vi)


# 🔁 Convertir listas a arrays de NumPy

x = np.array(x)

y = np.array(y)

V = np.array(V)


# 📏 Obtener dimensiones únicas para reconstruir la malla

nx = len(np.unique(x))

ny = len(np.unique(y))


X = x.reshape((nx, ny))

Y = y.reshape((nx, ny))

Z = V.reshape((nx, ny))


# 🎨 Crear figura 3D

fig = plt.figure(figsize=(10, 7))

ax = fig.add_subplot(111, projection='3d')


# 🟠 Superficie coloreada

surf = ax.plot_surface(X, Y, Z, cmap="plasma", edgecolor='k', linewidth=0.4, alpha=0.9)


# 🎯 Curvas de nivel proyectadas en z=0

ax.contour(X, Y, Z, zdir='z', offset=0, cmap="plasma")


# 🧭 Etiquetas

ax.set_xlabel('X (m)')

ax.set_ylabel('Y (m)')

ax.set_zlabel('Voltaje (V)')

ax.set_title("Distribución de Voltaje - Ecuación de Laplace")


# 🎨 Barra de colores

fig.colorbar(surf, shrink=0.6, aspect=10, label='Voltaje (V)')


# 💾 Guardar imagen automáticamente

plt.tight_layout()

plt.savefig("resultado_python.png", dpi=300)


# 👁️ Mostrar en pantalla

plt.show() 
