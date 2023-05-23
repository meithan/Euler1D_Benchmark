# ============================================================================
# Euler1D.jl
#
# By J.C. Toledo-Roy
# 18 / may / 2023
#
# Programa minimalista que resuelve el tubo de choque de Sod para las
# ecuaciones de Euler en 1D usando el método de Lax-Friedichs
#
# ==============================================================================

# ==============================================================================
# CONSTANTES Y VARIABLES GLOBALES
# ==============================================================================

# ----------------------------------------------------------------------------
# CONSTANTES Y VARIABLES GLOBALES ESPECIFICADAS POR EL USUARIO

# Parámetros de la simulación
const NEQ = 3             # Número de ecuaciones
const NX = 5000            # Tamaño de la malla
const XL = 0.0            # Coordenada física del extremo izquierdo
const XR = 1.0           # Coordenada física del extremo derecho
const TFIN = 0.20         # Tiempo final de integración
const CFL = 0.9           # Parametro de Courant
const DTOUT = TFIN/10     # Intervalo para escribir a disco
const GAMMA = 1.4         # Razón de capacidades caloríficas

# Directorio donde escribir las salidas (usar "./" para dir actual)
# Debe terminar en una diagonal '/'
const OUT_DIR = "./temp/"

const do_output = false

# ============================================================================
# NO ES NECESARIO MODIFICAR NADA DEBAJO DE ESTA LÍNEA
# ============================================================================

# Constantes calculadas de las anteriores
DX::Float64 = (XR-XL)/NX      # Espaciamiento de la malla

# ----------------------------------------------------------------------------
# VARIABLES GLOBALES

U = zeros(Float64, (NX+2, NEQ))      # Variables conservadas actuales
UP = zeros(Float64, (NX+2, NEQ))     # Variables conservadas "avanzadas"
F = zeros(Float64, (NX+2, NEQ))      # Flujos físicos
UT = zeros(Float64, (NX+2, NEQ))     # Us temporales (sólo usadas en Macormack)
PRIM = zeros(Float64, (NX+2, NEQ))   # Variables primitivas

# Variables globales de la simulación
dt::Float64 = 0             # Paso de tiempo
t::Float64 = 0              # Tiempo actual
it::Int64 = 0               # Iteración actual
nout::Int64 = 0             # Número de la siguiente salida
tout::Float64 = 0           # Tiempo para el siguiente output
start_time::Float64 = 0     # Tiempo de inicio

# ==============================================================================
# CONDICIONES INICIALES
# ==============================================================================

# Devuelve la coordenada X de la celda i
function xcoord(i)
  XL + i*DX
end

# ------------------------------------------------------------------------------
# Condición inicial
function initflow(U)
  
  x0 = 0.5
  rhoL = 1.0
  uL = 0.0
  pL = 1.0
  rhoR = 0.125
  uR = 0.0
  pR = 0.1

  U1L = rhoL
  U2L = rhoL*uL
  U3L = 0.5*rhoL*uL*uL + pL/(GAMMA-1)
  U1R = rhoR
  U2R = rhoR*uR
  U3R = 0.5*rhoR*uR*uR + pR/(GAMMA-1)

  for i in 1:NX+2
    x = xcoord(i)
    if (x <= x0)
      U[i,1] = U1L
      U[i,2] = U2L
      U[i,3] = U3L
    else
      U[i,1] = U1R
      U[i,2] = U2R
      U[i,3] = U3R
    end
  end

end

# ==============================================================================
# Inicializaciones de variables globales del código
function initmain()

  global t = 0.0
  global it = 0
  global tout = 0.0
  global nout = 0

end

# ==============================================================================

# Aplicar condiciones de frontera a celdas fantasma
# El arreglo pasado es al que aplicaremos las BCs
function boundary(U)
  
  # BC a la izquierda
  for e in 1:NEQ
    U[1,e] = U[2,e]
  end
  
  # BC a la derecha
  for e in 1:NEQ
    U[NX+2,e] = U[NX+1,e]
  end

end

# ==============================================================================

# Calcula las primitivas a partir de las U pasadas
# Esto incluye las celdas fantasma
function flow2prim(U, PRIM)

  for i in 1:NX+2
    PRIM[i,1] = U[i,1]
    PRIM[i,2] = U[i,2]/U[i,1]
    PRIM[i,3] = (GAMMA-1)*(U[i,3] - U[i,2]*U[i,2]/(2*U[i,1]))
  end

end

# ==============================================================================

# Calcular los flujos físicos F -- Ecuaciones de Euler 1D
# No olvidar calcular F en las celdas fantasma!
function fluxes(PRIM, F)
  
  for i in 1:NX+2
    F[i,1] = PRIM[i,1] * PRIM[i,2]
    F[i,2] = PRIM[i,1] * PRIM[i,2]*PRIM[i,2] + PRIM[i,3]
    F[i,3] = PRIM[i,2]*(0.5*PRIM[i,1]*PRIM[i,2]*PRIM[i,2] + GAMMA/(GAMMA-1)*PRIM[i,3])
  end

end

# ==============================================================================

# Calcula el paso de tiempo resultante de la condición CFL
function timestep(PRIM, dt)
  
  # Determinamos la velocidad máxima en la malla
  max_speed = 0.0
  for i in 2:NX+1
    cs = sqrt(GAMMA*PRIM[i,3]/PRIM[i,1])
    u_plus_cs = abs(PRIM[i,2]) + cs
    if (u_plus_cs > max_speed)
      max_speed = u_plus_cs
    end
  end

  # Condición de estabilidad CFL
  _dt = CFL * DX / max_speed

  return _dt

end

# ==============================================================================
# SOLVERS NUMÉRICOS

# Método de Lax-Friedrichs
function Lax(U, F, UP)
  
  for i in 2:NX+1
    for e in 1:NEQ
      UP[i,e] = (U[i+1,e] + U[i-1,e])/2 - dt/(2*DX) * (F[i+1,e] - F[i-1,e])
    end
  end

end

# ==============================================================================
# STEPPING

# Volcar las UPs sobre las Us "finales" del paso de tiempo -- sin viscosidad
function step(U, UP)
  
  for i in 1:NX+1
    for e in 1:NEQ
      U[i,e] = UP[i,e]
    end
  end

  # Avanzar el estado de la simulación
  global t = t + dt
  global it = it + 1

end

# ==============================================================================

# Escribe a disco el estado de la simulación
function output(PRIM)

  if do_output
  
    # Generar el nombre del archivo de salida
    fname = "output_" * string(nout) * ".txt"

    # Generar la ruta del archivo de salida (inc directorio)
    fpath = OUT_DIR * fname
    
    # Abrir el archivo para escritura
    open(fpath, "w") do file

        # Escribir los valores de PRIM al archivo (sólo celdas físicas)
        for i in 2:NX+1
          x = XL + i*DX
          line = string(x) * " " * string(PRIM[i,1]) * " " * string(PRIM[i,2]) * " " * string(PRIM[i,3]) * "\n"
          write(file, line)
        
        end 

        println("Se escribió salida ", nout)

    end

  end

  # Avanzar variables de output
  global nout = nout + 1
  global tout = DTOUT * nout

end

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
while t < TFIN

  # Calcular el paso de tiempo
  global dt = timestep(PRIM, dt)

  # Actualizar flujos físicos
  fluxes(PRIM, F)

  # Aplicar el método numérico para calcular las UP
  Lax(U, F, UP)

  # Aplicar condiciones de frontera a las UP recién calculadas
  boundary(UP)

  # Avanzar las U, con viscosidad para Macormack
  step(U, UP)

  # Actualizar las primitivas PRIM usando las nuevas U
  flow2prim(U, PRIM)

  # Escribir a disco
  if (t >= tout)
    output(PRIM)
  end

end

# Imprimir tiempo transcurrido


# ==============================================================================
