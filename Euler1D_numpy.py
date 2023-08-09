# ============================================================================
# Euler1D.py
#
# por J.C. Toledo-Roy
# 20 / oct / 2020
#
# Este programa resuelve las ecuaciones de Euler de la dinámica de gases en 1D
# usando el método de Lax o el de Macormack.
#
# Incluye algunos tubos de choque, como el de Sod, como pruebas.
#
# Versión implementada exclusivamente con arreglos y operaciones de Numpy
#
# ==============================================================================

import math
import sys
import time

import numpy as np

# ==============================================================================
# CONSTANTES Y VARIABLES GLOBALES
# ==============================================================================

# Constantes con nombre (a usar como opciones más abajo) -- No modificar

# Solvers numéricos
SOLVER_LAX = 1
SOLVER_MACORMACK = 2

# Problemas (ICs) -- ver descripciones en la siguiente sección
SHOCK_TUBE1 = 1
SHOCK_TUBE2 = 2
SHOCK_TUBE3 = 3
SHOCK_TUBE3A = 4
SHOCK_TUBE4 = 5
SHOCK_TUBE5 = 6

# Tipos de condiciones de frontera
BC_FREEFLOW = 1
BC_INFLOW = 2
BC_REFLECTIVE = 3
BC_PERIODIC = 4

# ----------------------------------------------------------------------------
# CONSTANTES Y VARIABLES GLOBALES ESPECIFICADAS POR EL USUARIO

# Parámetros de la simulación
NEQ = 3             # Número de ecuaciones
NX = 5000            # Tamaño de la malla
XL = 0.0            # Coordenada física del extremo izquierdo
XR = 1.0            # Coordenada física del extremo derecho
TFIN = 0.20         # Tiempo final de integración
CFL = 0.9           # Parametro de Courant
DTOUT = TFIN/10     # Intervalo para escribir a disco
GAMMA = 1.4         # Razón de capacidades caloríficas

# Viscosidad artficial, sólo para Macormack
ETA = 0.1

# Solucionador numérico: SOLVER_LAX o SOLVER_MACORMACK
NUM_SOLVER = SOLVER_LAX
# NUM_SOLVER = SOLVER_MACORMACK

# Problema a resolver, una de las constantes siguientes:
# SHOCK_TUBE1: el tubo de choque clásico de Sod
# SHOCK_TUBE2: choque fuerte: IC con contraste de presión de 10^4
# SHOCK_TUBE3: discontinuidad de contact inmóvil
# SHOCK_TUBE3A: como anterior, pero DC se mueve lentamente
# SHOCK_TUBE4: dos ondas de rarefacción
# SHOCK_TUBE5: dos choques
PROBLEM = SHOCK_TUBE1

# Condiciones de frontera, usar las constantes siguientes:
# BC_FREEFLOW: salida libre (gradiente cero)
# BC_INFLOW: entrada; debe ser especificada en boundary()
# BC_REFLECTIVE: reflectiva
# BC_PERIODIC: periódica
BC_LEFT = BC_FREEFLOW
BC_RIGHT = BC_FREEFLOW

# Directorio donde escribir las salidas (usar "./" para dir actual)
# Debe terminar en una diagonal '/'
OUT_DIR = "./temp/"

# La "plantilla" para los nombres de archivos de salida, incluyendo
# el formato para el número de salida
# Por ejemplo, "output_%02i.txt" usará dos cifras (con un 0 si necesario)
OUT_FNAME = "output_%02i.txt"

do_output = False

# ============================================================================
# NO ES NECESARIO MODIFICAR NADA DEBAJO DE ESTA LÍNEA
# ============================================================================

# Constantes calculadas de las anteriores
DX = (XR-XL)/NX      # Espaciamiento de la malla

xcoords = XL + np.arange(NX+2) * DX

# ----------------------------------------------------------------------------
# VARIABLES GLOBALES

U = np.empty((NEQ,NX+2))          # Variables conservadas actuales
UP = np.empty((NEQ,NX+2))         # Variables conservadas "avanzadas"
F = np.empty((NEQ,NX+2))          # Flujos físicos
UT = np.empty((NEQ,NX+2))         # Us temporales (sólo usadas en Macormack)
PRIM = np.empty((NEQ,NX+2))       # Variables primitivas

# Variables globales de la simulación
dt = 0             # Paso de tiempo
t = 0              # Tiempo actual
it = 0             # Iteración actual
nout = 0           # Número de la siguiente salida
tout = 0           # Tiempo para el siguiente output
start_time = 0     # Tiempo de inicio

# ==============================================================================
# CONDICIONES INICIALES
# ==============================================================================

# ------------------------------------------------------------------------------
# Condición inicial para imponer tubos de choque genéricos
def shocktube_IC(U, x0, rhoL, uL, pL, rhoR, uR, pR):
  
  U1L = rhoL
  U2L = rhoL*uL
  U3L = 0.5*rhoL*uL*uL + pL/(GAMMA-1)
  U1R = rhoR
  U2R = rhoR*uR
  U3R = 0.5*rhoR*uR*uR + pR/(GAMMA-1)

  left = xcoords <= x0
  U[0,:] = np.where(left, U1L, U1R)
  U[1,:] = np.where(left, U2L, U2R)
  U[2,:] = np.where(left, U3L, U3R)

# ==============================================================================
# Impone las condiciones iniciales
# Debe inicializar los valores de U en todo el dominio
# Nótese que también debemos llenar las celdas fantasma
def initflow(U):

  # Inicializar variables conservadas U de acuerdo al problema elegido
  if (PROBLEM == SHOCK_TUBE1):
    # Tubo de choque clásico de Sod
    shocktube_IC(U, 0.5, 1.0, 0.0, 1.0, 0.125, 0.0, 0.1)
  elif (PROBLEM == SHOCK_TUBE2):
    # Choque fuerte
    shocktube_IC(U, 0.5, 1.0, 0.0, 1000.0, 1.0, 0.0, 0.1)
  elif (PROBLEM == SHOCK_TUBE3):
    # Discontinuidad de contacto estacionaria
    shocktube_IC(U, 0.5, 1.4, 0.0, 1.0, 1.0, 0.0, 1.0)
  elif (PROBLEM == SHOCK_TUBE3A):
    # Discontinuidad de contacto que se mueve lentamente
    shocktube_IC(U, 0.5, 1.4, 0.1, 1.0, 1.0, 0.1, 1.0)
  elif (PROBLEM == SHOCK_TUBE4):
    # Dos ondas de rarefacción
    shocktube_IC(U, 0.5, 1.0, -1.0, 0.4, 1.0, 1.0, 0.4)
  elif (PROBLEM == SHOCK_TUBE5):
    # Dos choques
    shocktube_IC(U, 0.5, 1.0, 1.0, 0.4, 1.0, -1.0, 0.4)
  else:
    print("Initial condition not set in initflow")
    sys.exit(1)

# ==============================================================================
# Inicializaciones de variables globales del código
def initmain():

  global t, it, nout, tout

  t = 0.0
  it = 0
  tout = 0.0
  nout = 0

# ==============================================================================

# Aplicar condiciones de frontera a celdas fantasma
# El arreglo pasado es al que aplicaremos las BCs
def boundary(U):
  
  # BC a la izquierda
  if (BC_LEFT == BC_FREEFLOW):
    U[:,0] = U[:,1]
  elif (BC_LEFT == BC_INFLOW):
    U[0,0] = 0.0
    U[1,0] = 0.0
    U[2,0] = 0.0
  elif (BC_LEFT == BC_REFLECTIVE):
    U[0,0] =  U[0,1]
    U[1,0] = -U[1,1]
    U[2,0] =  U[2,1]
  elif (BC_LEFT == BC_PERIODIC):
    U[:,0] = U[:,NX]

  # BC a la derecha
  if (BC_RIGHT == BC_FREEFLOW):
    U[:,NX+1] = U[:,NX]
  elif (BC_RIGHT == BC_INFLOW):
    U[0,NX+1] = 0.0
    U[1,NX+1] = 0.0
    U[2,NX+1] = 0.0
  elif (BC_RIGHT == BC_REFLECTIVE):
    U[0,NX+1] =  U[0,NX]
    U[1,NX+1] = -U[1,NX]
    U[2,NX+1] =  U[2,NX]
  elif (BC_RIGHT == BC_PERIODIC):
    U[:,NX+1] = U[:,1]

# ==============================================================================

# Calcula las primitivas a partir de las U pasadas
# Esto incluye las celdas fantasma
def flow2prim(U, PRIM):

  PRIM[0,:] = U[0,:]
  PRIM[1,:] = U[1,:]/U[0,:]
  PRIM[2,:] = (GAMMA-1)*(U[2,:] - U[1,:]**2/(2*U[0,:]))

# ==============================================================================

# Calcular los flujos físicos F -- Ecuaciones de Euler 1D
# No olvidar calcular F en las celdas fantasma!
def fluxes(PRIM, F):
  
  F[0,:] = PRIM[0,:] * PRIM[1,:]
  F[1,:] = PRIM[0,:] * PRIM[1,:]**2 + PRIM[2,:]
  F[2,:] = PRIM[1,:] * (0.5*PRIM[0,:]*PRIM[1,:]**2 + GAMMA/(GAMMA-1)*PRIM[2,:])

# ==============================================================================

# Calcula el paso de tiempo resultante de la condición CFL
def timestep(PRIM):
  
  # Determinamos la velocidad máxima en la malla
  cs = np.sqrt(GAMMA*PRIM[2,:]/PRIM[0,:])
  u_plus_cs = np.abs(PRIM[1,:]) + cs
  max_speed = u_plus_cs.max()      

  # Condición de estabilidad CFL
  dt = CFL * DX / max_speed

  return dt

# ==============================================================================
# SOLVERS NUMÉRICOS

# Método de Lax-Friedrichs
def Lax(U, F, UP):
  
  UP[:,1:NX+1] = (U[:,2:NX+2] + U[:,0:NX])/2 - dt/(2*DX) * (F[:,2:NX+2] - F[:,0:NX])

# -----------------------------------------------------------------------------

# Método de Macormack
def Macormack(U, F, UP):
  
  # Paso predictor: actualizamos las UT con flujos hacia adelante
  UT[:,1:NX+1] = U[:,1:NX+1] - dt/DX * (F[:,2:NX+2] - F[:,1:NX+1])

  # Aplicamos las BCs a las UT
  boundary(UT)

  # Actualizar las primitivas de nuevo usando las UT esta vez
  calc_prims(UT, PRIM)

  # Re-calculamos los flujos F pero usando las primitivas actualizadas
  fluxes(PRIM, F)

  # Paso corrector: obtenemos las UP usando U, UT y F actualizados
  UP[:,1:NX+1] = (U[:,1:NX+1] + UT[:,1:NX+1])/2 - dt/(2*DX) * (F[:,1:NX+1] - F[:,0:NX])

# -----------------------------------------------------------------------------
# "Wrapper" para el solver numérico: se selecciona el apropiado
def solver(U, F, UP):
  
  if (NUM_SOLVER == SOLVER_LAX):
    Lax(U, F, UP)
  elif (NUM_SOLVER == SOLVER_MACORMACK):
    Macormack(U, F, UP)

# ==============================================================================
# STEPPING

# Volcar las UPs sobre las Us "finales" del paso de tiempo -- sin viscosidad
def step_simple(U, UP):
  
  U[:,:] = UP[:,:]

# ------------------------------------------------------------------------------

# Volcar las UPs sobre las Us "finales" del paso de tiempo, aplicando
# viscosidad artificial donde haya máximos o mínimos locales
def step_viscosity(U, UP):

  raise(Exception(NotImplementedError()))

  for i in range(1,NX+1):
    for e in range(NEQ):
      # Aplicamos la viscosidad sólo donde hay mínimos/máximos locales
      # En las demás celdas simplemente copiamos las UP sobre las U
      # ETA debe ser estrictamente menor que 1/2
      if ((U[i+1][e]-U[i][e])*(U[i][e]-U[i-1][e]) < 0):
        U[i][e] = UP[i][e] + ETA*(UP[i+1][e] + UP[i-1][e] - 2*UP[i][e])
      else:
        U[i][e] = UP[i][e]

  # Lo de arriba no toca las dos celdas fantasma, pero éstas también
  # deben ser copiadas (sin viscosidad)
  for e in range(NEQ):
    U[e,0] = UP[0][e]
    U[NX+1][e] = UP[NX+1][e]

# ------------------------------------------------------------------------------
# "Wrapper" para el stepping, dependiendo del solver usado
def step(U, UP):

  global t, it
    
  if (NUM_SOLVER == SOLVER_LAX):
    step_simple(U, UP)
  elif (NUM_SOLVER == SOLVER_MACORMACK):
    step_viscosity(U, UP)

  # Avanzar el estado de la simulación
  t = t + dt
  it = it + 1

# ==============================================================================

# Escribe a disco el estado de la simulación
def output(PRIM):

  global nout, tout

  if do_output:
  
    # Generar el nombre del archivo de salida
    fname = OUT_FNAME % nout

    # Generar la ruta del archivo de salida (inc directorio)
    fpath = OUT_DIR
    if OUT_DIR != "" and not OUT_DIR.endswith("/"):
      fpath = OUT_DIR + "/"
    fpath = fpath + fname

      # Abrir el archivo para escritura
    fout = open(fpath, "w")

    # Escribir los valores de U al archivo (sólo celdas físicas)
    for i in range(1, NX+1):
      fout.write("{} {} {} {}\n".format(xcoords[i], PRIM[0,i], PRIM[1,i], PRIM[2,i]))

    # Cerrar archivo
    fout.close()

    print("Se escribió salida", nout)

  # Avanzar variables de output
  nout = nout + 1
  tout = DTOUT * nout

# ==============================================================================
# Programa principal
# ==============================================================================

# Condición inicial e inicializaciones
initmain()

# Condición inicial e inicializaciones
initflow(U)

# Llenar celdas fantasma
boundary(U)

# Calcular las primitivas
flow2prim(U, PRIM)

# Escribir condición inicial a disco
output(PRIM)

# Bucle principal
clock_start = time.time()
while (t < TFIN):

  # Calcular el paso de tiempo
  dt = timestep(PRIM)

  # Actualizar flujos físicos
  fluxes(PRIM, F)

  # Aplicar el método numérico para calcular las UP
  solver(U, F, UP)

  # Aplicar condiciones de frontera a las UP recién calculadas
  boundary(UP)

  # Avanzar las U, con viscosidad para Macormack
  step(U, UP)

  # Actualizar las primitivas PRIM usando las nuevas U
  flow2prim(U, PRIM)

  # Escribir a disco
  if (t >= tout):
    output(PRIM)

# Imprimir tiempo transcurrido
elapsed = time.time() - clock_start
# print("Se calcularion {} iteraciones {:.3f} s".format(it,elapsed))
print(elapsed)

# ==============================================================================
