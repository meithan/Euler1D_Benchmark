/*------------------------------------------------------------------------------
Euler1D.cpp

por J.C. Toledo-Roy
20 Oct 2020

Resuelve las ecuaciones de Euler de la dinámica de gases en 1D
usando el método de Lax o el de Macormack.

Incluye algunos tubos de choque, como el de Sod, como pruebas.

------------------------------------------------------------------------------*/

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <time.h>
using namespace std;

/*============================================================================*/
// Constantes con nombre (a usar como opciones más abajo) -- No modificar

// Solucionadores numéricos
const int SOLVER_LAX = 1;
const int SOLVER_MACORMACK = 2;

// Problemas (ICs) -- ver descripciones en la siguiente sección
const int SHOCK_TUBE_SOD = 1;
const int SHOCK_TUBE1 = 1;
const int SHOCK_TUBE2 = 2;
const int SHOCK_TUBE3 = 3;
const int SHOCK_TUBE3A = 4;
const int SHOCK_TUBE4 = 5;
const int SHOCK_TUBE5 = 6;

// Tipos de condiciones de frontera
const int BC_FREEFLOW = 1;
const int BC_INFLOW = 2;
const int BC_REFLECTIVE = 3;
const int BC_PERIODIC = 4;

/*============================================================================*/
// CONSTANTES Y VARIABLES GLOBALES ESPECIFICADAS POR EL USUARIO

// Parámetros de la simulación
const int    NEQ = 3;          // Número de ecuaciones
const int    NX = 5000;         // Tamaño de la malla
const double XL = 0.0;         // Coordenada física del extremo izquierdo
const double XR = 1.0;         // Coordenada física del extremo derecho
const double TFIN = 0.20;       // Tiempo final de integración
const double CFL = 0.9;        // Parametro de Courant
const double DTOUT = TFIN/10;     // Intervalo para escribir a disco
const double GAMMA = 1.4;      // Razón de capacidades caloríficas

// Viscosidad artficial; sólo para Macormack
const double ETA = 0.05;

// Solucionador numérico: SOLVER_LAX o SOLVER_MACORMACK
const int NUM_SOLVER = SOLVER_LAX;
// const int NUM_SOLVER = SOLVER_MACORMACK;

// Problema a resolver: una de las constantes siguientes:
// SHOCK_TUBE1: el tubo de choque clásico de Sod
// SHOCK_TUBE2: choque fuerte: IC con contraste de presión de 10^4
// SHOCK_TUBE3: discontinuidad de contact inmóvil
// SHOCK_TUBE3A: como anterior, pero DC se mueve lentamente
// SHOCK_TUBE4: dos ondas de rarefacción
// SHOCK_TUBE5: dos choques
const int PROBLEM = SHOCK_TUBE_SOD;

// Condiciones de frontera
// BC_FREEFLOW: salida libre (gradiente cero)
// BC_INFLOW: entrada; debe ser especificada en boundary()
// BC_REFLECTIVE: reflectiva
// BC_PERIODIC: periódica
const int BC_LEFT = BC_FREEFLOW;
const int BC_RIGHT = BC_FREEFLOW;

// Directorio donde escribir las salidas (usar "./" para dir actual)
// Debe terminar en una diagonal '/'
const char* OUT_DIR = "./temp/";

const bool do_output = false;

/*============================================================================*/
/* NO ES NECESARIO MODIFICAR NADA DEBAJO DE ESTA LÍNEA                        */
/*============================================================================*/

// Constantes calculadas de las anteriores

const double DX = (XR-XL)/NX;      // Espaciamiento de la malla

/*============================================================================*/
// VARIABLES GLOBALES

double    U[NX+2][NEQ];      // Variables conservadas actuales
double   UP[NX+2][NEQ];      // Variables conservadas "avanzadas"
double    F[NX+2][NEQ];      // Flujos físicos
double   UT[NX+2][NEQ];      // Us temporales; sólo usadas en Macormack
double PRIM[NX+2][NEQ];      // Variables primitivas

double dt;         // Paso de tiempo
double t;          // Tiempo actual
int it;            // Iteración actual
clock_t start;     // Tiempo de inicio
double tout;       // Tiempo de la siguiente salida
int nout;          // Número de la siguiente salida

/*============================================================================*/
// CONDICIONES INICIALES

// Devuelve la coordenada X de la celda i
double xcoord(int i) {
  return XL + i*DX;
}

// -----------------------------------------------------------------------------

// Condición inicial para un tubo de choque genéricos
void shocktube_IC(double U[NX+2][NEQ], double x0, double rhoL, double uL, double pL, double rhoR, double uR, double pR) {

  const double U1L = rhoL;
  const double U2L = rhoL*uL;
  const double U3L = 0.5*rhoL*uL*uL + pL/(GAMMA-1);
  const double U1R = rhoR;
  const double U2R = rhoR*uR;
  const double U3R = 0.5*rhoR*uR*uR + pR/(GAMMA-1);
  double x;

  for (int i = 1; i <= NX; i++) {
    x = xcoord(i);
    if (x <= x0) {
      U[i][0] = U1L;
      U[i][1] = U2L;
      U[i][2] = U3L;
    } else {
      U[i][0] = U1R;
      U[i][1] = U2R;
      U[i][2] = U3R;
    }
  }

  // Descomentar esto para eliminar el artefacto de "escalones"
  // en el método de Lax
  // for (int e = 0; e <= 2; e++) {
  //   UU[e][NX/2] = 0.5 * (UU[e][NX/2+1] + UU[e][NX/2-1]);
  // }

}

// -----------------------------------------------------------------------------

// Impone las condiciones iniciales
void initflow(double UU[NX+2][NEQ]) {

  // Inicializar variables conservadas U de acuerdo al problema elegido
  if (PROBLEM == SHOCK_TUBE1) {
    // Tubo de choque clásico de Sod
    shocktube_IC(UU, 0.5, 1.0, 0.0, 1.0, 0.125, 0.0, 0.1);
  } else if (PROBLEM == SHOCK_TUBE2) {
    // Choque fuerte
    shocktube_IC(UU, 0.5, 1.0, 0.0, 1000.0, 1.0, 0.0, 0.1);
  } else if (PROBLEM == SHOCK_TUBE3) {
    // Discontinuidad de contacto estacionaria
    shocktube_IC(UU, 0.5, 1.4, 0.0, 1.0, 1.0, 0.0, 1.0);
  } else if (PROBLEM == SHOCK_TUBE3A) {
    // Discontinuidad de contacto que se mueve lentamente
    shocktube_IC(UU, 0.5, 1.4, 0.1, 1.0, 1.0, 0.1, 1.0);
  } else if (PROBLEM == SHOCK_TUBE4) {
    // Dos ondas de rarefacción
    shocktube_IC(UU, 0.5, 1.0, -1.0, 0.4, 1.0, 1.0, 0.4);
  } else if (PROBLEM == SHOCK_TUBE5) {
    // Dos choques
    shocktube_IC(UU, 0.5, 1.0, 1.0, 0.4, 1.0, -1.0, 0.4);
  } else {
    printf("Initial condition not set in initflow\n");
    exit(1);
  }

}

// -----------------------------------------------------------------------------

// Inicializaciones generales, como variables globales
void initmain() {

  t = 0;
  it = 0;
  nout = 0;
  tout = 0;

}

/*============================================================================*/

// Aplica condiciones de frontera a celdas fantasma del arreglo pasado
void boundary(double U[NX+2][NEQ]) {

  // BC a la izquierda
  if (BC_LEFT == BC_FREEFLOW) {
    for (int e = 0; e < NEQ; e++) {
      U[0][e] = U[1][e];
    }
  } else if (BC_LEFT == BC_INFLOW) {
    U[0][0] = 0.0;
    U[0][1] = 0.0;
    U[0][2] = 0.0;
  } else if (BC_LEFT == BC_REFLECTIVE) {
    U[0][0] = U[1][0];
    U[0][1] = -U[1][1];
    U[0][2] = U[1][2];
  } else if (BC_LEFT == BC_PERIODIC) {
    for (int e = 0; e < NEQ; e++) {
      U[0][e] = U[NX][e];
    }
  }

  // BC a la derecha
  if (BC_RIGHT == BC_FREEFLOW) {
    for (int e = 0; e < NEQ; e++) {
      U[NX+1][e] = U[NX][e];
    }
  } else if (BC_RIGHT == BC_INFLOW) {
    U[NX+1][0] = 0.0;
    U[NX+1][1] = 0.0;
    U[NX+1][2] = 0.0;
  } else if (BC_RIGHT == BC_REFLECTIVE) {
    U[NX+1][0] = U[NX][0];
    U[NX+1][1] = -U[NX][1];
    U[NX+1][2] = U[NX][2];
  } else if (BC_RIGHT == BC_PERIODIC) {
    for (int e = 0; e < NEQ; e++) {
      U[NX+1][e] = U[0][e];
    }
  }

}
/*============================================================================*/

// Calcula las primitivas a partir de las U pasadas
// Esto incluye las celdas fantasma
void calc_prims(double U[NX+2][NEQ], double PRIM[NX+2][NEQ]) {

  for (int i = 0; i <= NX+1; i++) {
    PRIM[i][0] = U[i][0];
    PRIM[i][1] = U[i][1]/U[i][0];
    PRIM[i][2] = (GAMMA-1)*(U[i][2] - U[i][1]*U[i][1]/(2*U[i][0]));
  }

}

/*============================================================================*/

// Calcular los flujos físicos F -- Ecuaciones de Euler 1D
// No olvidar calcular F en las celdas fantasma!
void fluxes(double PRIM[NX+2][NEQ], double F[NX+2][NEQ]) {

  for (int i = 0; i <= NX+1; i++) {
    F[i][0] = PRIM[i][0] * PRIM[i][1];
    F[i][1] = PRIM[i][0] * PRIM[i][1]*PRIM[i][1] + PRIM[i][2];
    F[i][2] = PRIM[i][1]*(0.5*PRIM[i][0]*PRIM[i][1]*PRIM[i][1] + GAMMA/(GAMMA-1) * PRIM[i][2]);
  }

}

/*============================================================================*/

// Calcula la velocidad del sonido dada la presión y densidad
double sound_speed(double pres, double rho) {
  return sqrt(GAMMA*pres/rho);
}

/*============================================================================*/

// Calcula el paso de tiempo resultante de la condición CFL
double timestep(double P[NX+2][NEQ]) {

  // Determinamos la velocidad máxima en la malla
  double cs, u_plus_cs;
  double max_speed = 0.0;
  for (int i = 1; i <= NX; i++) {
    cs = sound_speed(P[i][2], P[i][0]);
    u_plus_cs = fabs(P[i][1]) + cs;
    if (u_plus_cs > max_speed) max_speed = u_plus_cs;
  }

  // Condición de estabilidad CFL
  return CFL * DX / max_speed;

}

/*============================================================================*/
// SOLVERS NUMÉRICOS

// -----------------------------------------------------------------------------
// Método de Lax-Friedrichs
void Lax(double U[NX+2][NEQ], double F[NX+2][NEQ], double UP[NX+2][NEQ]) {

  // Aplicar método de Lax
  for (int i = 1; i <= NX; i++) {
    for (int e = 0; e < NEQ; e++) {
      UP[i][e] = (U[i+1][e] + U[i-1][e])/2.0
                 - dt/(2*DX) * (F[i+1][e] - F[i-1][e]);
    }
  }

}

// -----------------------------------------------------------------------------
// Método de Macormack
void Macormack(double U[NX+2][NEQ], double F[NX+2][NEQ], double UP[NX+2][NEQ]) {

  // Paso predictor: actualizamos las UT con flujos hacia adelante
  for (int i = 1; i <= NX; i++) {
    for (int e = 0; e < NEQ; e++) {
      UT[i][e] = U[i][e] - dt/DX * (F[i+1][e] - F[i][e]);
    }
  }

  // Aplicamos las BCs a las UT
  boundary(UT);

  // Actualizar las primitivas de nuevo usando las UT esta vez
  calc_prims(UT, PRIM);

  // Re-calculamos los flujos F pero usando las primitivas actualizadas
  fluxes(PRIM, F);

  // Paso corrector: obtenemos las UP usando U, UT y F actualizados
  for (int i = 1; i <= NX; i++) {
    for (int e = 0; e < NEQ; e++) {
      UP[i][e] = (U[i][e] + UT[i][e])/2 - dt/(2*DX) * (F[i][e] - F[i-1][e]);
    }
  }

}

// -----------------------------------------------------------------------------
// "Wrapper" para el solver numérico: se selecciona el apropiado
void solver(double U[NX+2][NEQ], double F[NX+2][NEQ], double UP[NX+2][NEQ]) {
  
  if (NUM_SOLVER == SOLVER_LAX) {
    Lax(U, F, UP);
  } else if (NUM_SOLVER == SOLVER_MACORMACK) {
    Macormack(U, F, UP);
  }

}

/*============================================================================*/
// STEPPING

// Volcar las UPs sobre las Us "finales" del paso de tiempo -- sin viscosidad
// Incluye las celdas fantasma
void step_simple(double U[NX+2][NEQ], double UP[NX+2][NEQ]) {

  // Esto incluye las celdas fantasma
  for (int i = 0; i <= NX+1; i++) {
    for (int e = 0; e < NEQ; e++) {
      U[i][e] = UP[i][e];
    }
  }

}

// -----------------------------------------------------------------------------
// Volcar las UPs sobre las Us "finales" del paso de tiempo, aplicando
// viscosidad artificial donde haya máximos o mínimos locales
// Incluye las celdas fantasma
void step_viscosity(double U[NX+2][NEQ], double UP[NX+2][NEQ]) {

  for (int i = 1; i <= NX; i++) {
    for (int e = 0; e < NEQ; e++) {
      // Aplicamos la viscosidad sólo donde hay mínimos/máximos locales
      // En las demás celdas simplemente copiamos las UP sobre las U
      // ETA debe ser estrictamente menor que 1/2
      if ((U[i+1][e]-U[i][e])*(U[i][e]-U[i-1][e]) < 0) {
        U[i][e] = UP[i][e] + ETA*(UP[i+1][e] + UP[i-1][e] - 2*UP[i][e]);
      } else {
        U[i][e] = UP[i][e];
      }
    }
  }

  // Lo de arriba no toca las dos celdas fantasma, pero éstas también
  // deben ser copiadas (sin viscosidad)
  for (int e = 0; e < NEQ; e++) {
    U[0][e] = UP[0][e];
    U[NX+1][e] = UP[NX+1][e];
  }

}

// -----------------------------------------------------------------------------
// "Wrapper" para el stepping, dependiendo del solver usado
void step(double U[NX+2][NEQ], double UP[NX+2][NEQ]) {
  
  if (NUM_SOLVER == SOLVER_LAX || ETA == 0) {
    step_simple(U, UP);
  } else if (NUM_SOLVER == SOLVER_MACORMACK) {
    step_viscosity(U, UP);
  }

}

/*============================================================================*/

// Escribe a disco el estado de la simulación
void output(double PRIM[NX+2][NEQ]) {

  if (do_output) {

    // Generar el nombre del archivo de salida
    char fname[80];
    sprintf(fname, "%s/output_%02i.txt", OUT_DIR, nout);

    // Abrir el archivo
    fstream fout(fname, ios::out);

    // Escribir la posición, densidad, velocidad y presión a disco, separadas
    // por espacios
    double x;
    for (int i = 0; i <= NX; i++) {
      x = xcoord(i);
      fout << x << " " << PRIM[i][0] << " " << PRIM[i][1] << " " << PRIM[i][2] << endl;
    }

    // Cerrar archivo
    fout.close();

    printf("Salida %i, t=%.7f, it=%i, dt=%.7f\n", nout, t, it, dt);
  
  }

  // Avanzar variables de output
  nout = nout + 1;
  tout = nout * DTOUT;

}

/*============================================================================*/
// PROGRAMA PRINCIPAL

int main() {

  // Inicializaciones generales
  initmain();

  // Condición inicial
  initflow(U);

  // Llenar celdas fantasma
  boundary(U);

  // Actualizar las primitivas
  calc_prims(U, PRIM);

  // Escribir condición inicial a disco
  output(PRIM);

  // Tiempo de inicio de la simulación
  start = clock();

  // BUCLE PRINCIPAL
  while (t < TFIN) {

    // Calcular el paso de tiempo dt
    dt = timestep(PRIM);

    // Actualizar los flujos F
    fluxes(PRIM, F);

    // Aplicar el método numérico para calcular las UP
    solver(U, F, UP);

    // Aplicar condiciones de frontera a las UP recién calculadas
    boundary(UP);

    // Avanzar las U, con viscosidad para Macormack
    step(U, UP);

    // Actualizar las primitivas PRIM usando las nuevas U
    calc_prims(U, PRIM);

    // Avanzar el estado de la simulación
    t = t + dt;
    it = it + 1;

    // Escribir a disco, si aplicable
    if (t >= tout) {
      output(PRIM);
    }

  }

  // Ejecución terminada!
  printf("\nSe calcularon %i iteraciones en %.3f s\n\n", it, (double)(clock() - start)/CLOCKS_PER_SEC);

}

/*============================================================================*/
