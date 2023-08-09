/*==============================================================================
Euler1D.java

por J.C. Toledo-Roy
9 jun 2023

Versión en Java del código original.

Resuelve las ecuaciones de Euler de la dinámica de gases en 1D
usando el método de Lax o el de Macormack.

==============================================================================*/

import java.lang.Math;
import java.io.File;
import java.io.FileWriter;
import java.io.BufferedWriter;

public class Euler1D {

  /*==========================================================================*/
  // Constantes con nombre (a usar como opciones más abajo) -- No modificar

  // Solucionadores numéricos
  static final int SOLVER_LAX = 1;
  static final int SOLVER_MACORMACK = 2;

  // Problemas (ICs) -- ver descripciones en la siguiente sección
  static final int SHOCK_TUBE_SOD = 1;
  static final int SHOCK_TUBE1 = 1;
  static final int SHOCK_TUBE2 = 2;
  static final int SHOCK_TUBE3 = 3;
  static final int SHOCK_TUBE3A = 4;
  static final int SHOCK_TUBE4 = 5;
  static final int SHOCK_TUBE5 = 6;

  // Tipos de condiciones de frontera
  static final int BC_FREEFLOW = 1;
  static final int BC_INFLOW = 2;
  static final int BC_REFLECTIVE = 3;
  static final int BC_PERIODIC = 4;

  /*==========================================================================*/
  // CONSTANTES Y VARIABLES GLOBALES ESPECIFICADAS POR EL USUARIO

  // Parámetros de la simulación
  static final int    NEQ = 3;          // Número de ecuaciones
  static final int    NX = 5000;         // Tamaño de la malla
  static final double XL = 0.0;         // Coordenada física del extremo izquierdo
  static final double XR = 1.0;         // Coordenada física del extremo derecho
  static final double TFIN = 0.20;       // Tiempo final de integración
  static final double CFL = 0.9;        // Parametro de Courant
  static final double DTOUT = TFIN/10;     // Intervalo para escribir a disco
  static final double GAMMA = 1.4;      // Razón de capacidades caloríficas

  // Viscosidad artficial; sólo para Macormack
  static final double ETA = 0.05;

  // Solucionador numérico: SOLVER_LAX o SOLVER_MACORMACK
  static final int NUM_SOLVER = SOLVER_LAX;
  // final int NUM_SOLVER = SOLVER_MACORMACK;

  // Problema a resolver: una de las finalantes siguientes:
  // SHOCK_TUBE1: el tubo de choque clásico de Sod
  // SHOCK_TUBE2: choque fuerte: IC con contraste de presión de 10^4
  // SHOCK_TUBE3: discontinuidad de contact inmóvil
  // SHOCK_TUBE3A: como anterior, pero DC se mueve lentamente
  // SHOCK_TUBE4: dos ondas de rarefacción
  // SHOCK_TUBE5: dos choques
  static final int PROBLEM = SHOCK_TUBE_SOD;

  // Condiciones de frontera
  // BC_FREEFLOW: salida libre (gradiente cero)
  // BC_INFLOW: entrada; debe ser especificada en boundary()
  // BC_REFLECTIVE: reflectiva
  // BC_PERIODIC: periódica
  static final int BC_LEFT = BC_FREEFLOW;
  static final int BC_RIGHT = BC_FREEFLOW;

  // Directorio donde escribir las salidas (usar "./" para dir actual)
  // Debe terminar en una diagonal '/'
  static final String OUT_DIR = "./temp/";

  static final boolean do_output = false;

  /*==========================================================================*/
  /* NO ES NECESARIO MODIFICAR NADA DEBAJO DE ESTA LÍNEA                        */
  /*==========================================================================*/

  // finalantes calculadas de las anteriores

  static final double DX = (XR-XL)/NX;      // Espaciamiento de la malla

  /*==========================================================================*/
  // VARIABLES GLOBALES

  static double[][] U = new double[NX+2][NEQ];      // Variables conservadas actuales
  static double[][] UP = new double[NX+2][NEQ];      // Variables conservadas "avanzadas"
  static double[][] F = new double[NX+2][NEQ];      // Flujos físicos
  static double[][] UT = new double[NX+2][NEQ];      // Us temporales; sólo usadas en Macormack
  static double[][] PRIM = new double[NX+2][NEQ];      // Variables primitivas

  static double dt;         // Paso de tiempo
  static double t;          // Tiempo actual
  static int it;            // Iteración actual
  static double tout;       // Tiempo de la siguiente salida
  static int nout;          // Número de la siguiente salida

  /*==========================================================================*/
  // CONDICIONES INICIALES

  // Devuelve la coordenada X de la celda i
  static double xcoord(int i) {
    return XL + i*DX;
  }

  // ---------------------------------------------------------------------------

  // Condición inicial para un tubo de choque genéricos
  static void shocktube_IC(double[][] UU, double x0, double rhoL, double uL, double pL, double rhoR, double uR, double pR) {

    final double U1L = rhoL;
    final double U2L = rhoL*uL;
    final double U3L = 0.5*rhoL*uL*uL + pL/(GAMMA-1);
    final double U1R = rhoR;
    final double U2R = rhoR*uR;
    final double U3R = 0.5*rhoR*uR*uR + pR/(GAMMA-1);
    double x;

    for (int i = 1; i <= NX; i++) {
      x = xcoord(i);
      if (x <= x0) {
        UU[i][0] = U1L;
        UU[i][1] = U2L;
        UU[i][2] = U3L;
      } else {
        UU[i][0] = U1R;
        UU[i][1] = U2R;
        UU[i][2] = U3R;
      }
    }

    // Descomentar esto para eliminar el artefacto de "escalones"
    // en el método de Lax
    // for (int e = 0; e <= 2; e++) {
    //   UU[e][NX/2] = 0.5 * (UU[e][NX/2+1] + UU[e][NX/2-1]);
    // }

  }

  // ---------------------------------------------------------------------------

  // Impone las condiciones iniciales
  static void initflow(double[][] UU) {

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
      System.out.printf("Initial condition not set in initflow\n");
      System.exit(1);
    }

  }

  // ---------------------------------------------------------------------------

  // Inicializaciones generales, como variables globales
  static void initmain() {

    t = 0;
    it = 0;
    nout = 0;
    tout = 0;

  }

  /*==========================================================================*/

  // Aplica condiciones de frontera a celdas fantasma del arreglo pasado
  static void boundary(double[][] UU) {

    // BC a la izquierda
    if (BC_LEFT == BC_FREEFLOW) {
      for (int e = 0; e < NEQ; e++) {
        UU[0][e] = UU[1][e];
      }
    } else if (BC_LEFT == BC_INFLOW) {
      UU[0][0] = 0.0;
      UU[0][1] = 0.0;
      UU[0][2] = 0.0;
    } else if (BC_LEFT == BC_REFLECTIVE) {
      UU[0][0] = UU[1][0];
      UU[0][1] = -UU[1][1];
      UU[0][2] = UU[1][2];
    } else if (BC_LEFT == BC_PERIODIC) {
      for (int e = 0; e < NEQ; e++) {
        UU[0][e] = UU[NX][e];
      }
    }

    // BC a la derecha
    if (BC_RIGHT == BC_FREEFLOW) {
      for (int e = 0; e < NEQ; e++) {
        UU[NX+1][e] = UU[NX][e];
      }
    } else if (BC_RIGHT == BC_INFLOW) {
      UU[NX+1][0] = 0.0;
      UU[NX+1][1] = 0.0;
      UU[NX+1][2] = 0.0;
    } else if (BC_RIGHT == BC_REFLECTIVE) {
      UU[NX+1][0] = UU[NX][0];
      UU[NX+1][1] = -UU[NX][1];
      UU[NX+1][2] = UU[NX][2];
    } else if (BC_RIGHT == BC_PERIODIC) {
      for (int e = 0; e < NEQ; e++) {
        UU[NX+1][e] = UU[0][e];
      }
    }

  }

  /*==========================================================================*/

  // Calcula las primitivas a partir de las U pasadas
  // Esto incluye las celdas fantasma
  static void calc_prims(double[][] UU, double[][] PP) {

    for (int i = 0; i <= NX+1; i++) {
      PP[i][0] = UU[i][0];
      PP[i][1] = UU[i][1]/UU[i][0];
      PP[i][2] = (GAMMA-1)*(UU[i][2] - UU[i][1]*UU[i][1]/(2*UU[i][0]));
    }

  }

  /*==========================================================================*/

  // Calcular los flujos físicos F -- Ecuaciones de Euler 1D
  // No olvidar calcular F en las celdas fantasma!
  static void fluxes(double[][] PP, double[][] FF) {

    for (int i = 0; i <= NX+1; i++) {
      FF[i][0] = PP[i][0] * PP[i][1];
      FF[i][1] = PP[i][0] * PP[i][1]*PP[i][1] + PP[i][2];
      FF[i][2] = PP[i][1]*(0.5*PP[i][0]*PP[i][1]*PP[i][1] + GAMMA/(GAMMA-1) * PP[i][2]);
    }

  }

  /*==========================================================================*/

  // Calcula la velocidad del sonido dada la presión y densidad
  static double sound_speed(double pres, double rho) {
    return Math.sqrt(GAMMA*pres/rho);
  }

  /*==========================================================================*/

  // Calcula el paso de tiempo resultante de la condición CFL
  static double timestep(double[][] PP) {

    // Determinamos la velocidad máxima en la malla
    double cs, u_plus_cs;
    double max_speed = 0.0;
    for (int i = 1; i <= NX; i++) {
      cs = sound_speed(PP[i][2], PP[i][0]);
      u_plus_cs = Math.abs(PP[i][1]) + cs;
      if (u_plus_cs > max_speed) max_speed = u_plus_cs;
    }

    // Condición de estabilidad CFL
    return CFL * DX / max_speed;

  }

  /*==========================================================================*/
  // SOLVERS NUMÉRICOS

  // ---------------------------------------------------------------------------
  // Método de Lax-Friedrichs
  static void Lax(double[][] UU, double[][] FF, double[][] UUP) {

    // Aplicar método de Lax
    for (int i = 1; i <= NX; i++) {
      for (int e = 0; e < NEQ; e++) {
        UUP[i][e] = (UU[i+1][e] + UU[i-1][e])/2.0
        - dt/(2*DX) * (FF[i+1][e] - FF[i-1][e]);
      }
    }

  }

  // ---------------------------------------------------------------------------
  // Método de Macormack
  static void Macormack(double[][] UU, double[][] PP, double[][] FF, double[][] UUP) {

    // Paso predictor: actualizamos las UT con flujos hacia adelante
    for (int i = 1; i <= NX; i++) {
      for (int e = 0; e < NEQ; e++) {
        UT[i][e] = UU[i][e] - dt/DX * (FF[i+1][e] - FF[i][e]);
      }
    }

    // Aplicamos las BCs a las UT
    boundary(UT);

    // Actualizar las primitivas de nuevo usando las UT esta vez
    calc_prims(UT, PP);

    // Re-calculamos los flujos F pero usando las primitivas actualizadas
    fluxes(PP, F);

    // Paso corrector: obtenemos las UP usando U, UT y F actualizados
    for (int i = 1; i <= NX; i++) {
      for (int e = 0; e < NEQ; e++) {
        UUP[i][e] = (UU[i][e] + UT[i][e])/2 - dt/(2*DX) * (FF[i][e] - FF[i-1][e]);
      }
    }

  }

  // ---------------------------------------------------------------------------
  // "Wrapper" para el solver numérico: se selecciona el apropiado
  static void solver(double[][] UU, double[][] PP, double[][] FF, double[][] UUP) {
    
    if (NUM_SOLVER == SOLVER_LAX) {
      Lax(UU, FF, UUP);
    } else if (NUM_SOLVER == SOLVER_MACORMACK) {
      Macormack(UU, PP, FF, UUP);
    }

  }

  /*==========================================================================*/
  // STEPPING

  // Volcar las UPs sobre las Us "finales" del paso de tiempo -- sin viscosidad
  // Incluye las celdas fantasma
  static void step_simple(double[][] UU, double[][] UUP) {

    // Esto incluye las celdas fantasma
    for (int i = 0; i <= NX+1; i++) {
      for (int e = 0; e < NEQ; e++) {
        UU[i][e] = UUP[i][e];
      }
    }

  }

  // ---------------------------------------------------------------------------
  // Volcar las UPs sobre las Us "finales" del paso de tiempo, aplicando
  // viscosidad artificial donde haya máximos o mínimos locales
  // Incluye las celdas fantasma
  static void step_viscosity(double[][] UU, double[][] UUP) {

    for (int i = 1; i <= NX; i++) {
      for (int e = 0; e < NEQ; e++) {
        // Aplicamos la viscosidad sólo donde hay mínimos/máximos locales
        // En las demás celdas simplemente copiamos las UP sobre las U
        // ETA debe ser estrictamente menor que 1/2
        if ((UU[i+1][e]-UU[i][e])*(UU[i][e]-UU[i-1][e]) < 0) {
          UU[i][e] = UUP[i][e] + ETA*(UUP[i+1][e] + UUP[i-1][e] - 2*UUP[i][e]);
        } else {
          UU[i][e] = UUP[i][e];
        }
      }
    }

    // Lo de arriba no toca las dos celdas fantasma, pero éstas también
    // deben ser copiadas (sin viscosidad)
    for (int e = 0; e < NEQ; e++) {
      UU[0][e] = UUP[0][e];
      UU[NX+1][e] = UUP[NX+1][e];
    }

  }

  // ---------------------------------------------------------------------------
  // "Wrapper" para el stepping, dependiendo del solver usado
  static void step(double[][] UU, double[][] UUP) {
    
    if (NUM_SOLVER == SOLVER_LAX || ETA == 0) {
      step_simple(UU, UUP);
    } else if (NUM_SOLVER == SOLVER_MACORMACK) {
      step_viscosity(UU, UUP);
    }

  }

  /*==========================================================================*/

  // Escribe a disco el estado de la simulación
  static void output(double[][] PP) {

    if (do_output) {

      // Generar el nombre del archivo de salida
      String fname = String.format("%s/output_%02d.txt", OUT_DIR, nout);      

      // Abrir el archivo
      try {
        File file = new File(fname);
        BufferedWriter writer = new BufferedWriter(new FileWriter(file));

        // Escribir la posición, densidad, velocidad y presión a disco, separadas
        // por espacios
        double x;
        for (int i = 0; i <= NX; i++) {
          x = xcoord(i);
          writer.write(String.format("%f %f %f %f\n", x, PP[i][0], PP[i][1], PP[i][2]));
        }

        // Cerrar archivo
        writer.close();

        System.out.printf("Salida %d, t=%.7f, it=%d, dt=%.7f\n", nout, t, it, dt);

      } catch (Exception e) {
        e.getStackTrace();
        System.out.printf("Couldn'f write to disk\n");
        System.exit(1);
      }
 
    }

    // Avanzar variables de output
    nout = nout + 1;
    tout = nout * DTOUT;

  }

  /*============================================================================*/
  // PROGRAMA PRINCIPAL

  public static void main(String[] args) {

    long start = System.nanoTime();

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

    // BUCLE PRINCIPAL
    while (t < TFIN) {

      // Calcular el paso de tiempo dt
      dt = timestep(PRIM);

      // Actualizar los flujos F
      fluxes(PRIM, F);

      // Aplicar el método numérico para calcular las UP
      solver(U, PRIM, F, UP);

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

    double elapsed = (double)(System.nanoTime()-start)/1.0e9;
    System.out.printf("%.3f", elapsed);

  }

/*============================================================================*/

}