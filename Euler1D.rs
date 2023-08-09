#![allow(warnings)]
/*------------------------------------------------------------------------------
Euler1D.rs

por J.C. Toledo-Roy
19 May 2023

Versión en Rust del código original. Compilar con -C opt-level=3 para optimizar.

Resuelve las ecuaciones de Euler de la dinámica de gases en 1D
usando el método de Lax o el de Macormack.

Incluye algunos tubos de choque, como el de Sod, como pruebas.

------------------------------------------------------------------------------*/

use std::fs::File;
use std::io::Write;
use std::io::BufWriter;

/*============================================================================*/
// Constantes con nombre (a usar como opciones más abajo) -- No modificar

// Solucionadores numéricos
const SOLVER_LAX: u8 = 1;
const SOLVER_MACORMACK: u8 = 2;

// Problemas (ICs) -- ver descripciones en la siguiente sección
const SHOCK_TUBE_SOD: u8 = 1;
const SHOCK_TUBE1: u8 = 1;
const SHOCK_TUBE2: u8 = 2;
const SHOCK_TUBE3: u8 = 3;
const SHOCK_TUBE3A: u8 = 4;
const SHOCK_TUBE4: u8 = 5;
const SHOCK_TUBE5: u8 = 6;

// Tipos de condiciones de frontera
const BC_FREEFLOW: u8 = 1;
const BC_INFLOW: u8 = 2;
const BC_REFLECTIVE: u8 = 3;
const BC_PERIODIC: u8 = 4;

/*============================================================================*/
// CONSTANTES Y VARIABLES GLOBALES ESPECIFICADAS POR EL USUARIO

// Parámetros de la simulación
const NEQ: usize = 3;          // Número de ecuaciones
const NX: usize = 5000;         // Tamaño de la malla
const XL: f64 = 0.0;         // Coordenada física del extremo izquierdo
const XR: f64 = 1.0;         // Coordenada física del extremo derecho
const TFIN: f64 = 0.20;       // Tiempo final de integración
const CFL: f64 = 0.9;        // Parametro de Courant
const DTOUT: f64 = TFIN/10.0;     // Intervalo para escribir a disco
const GAMMA: f64 = 1.4;      // Razón de capacidades caloríficas

// Viscosidad artficial; sólo para Macormack
const ETA: f64 = 0.05;

// Solucionador numérico: SOLVER_LAX o SOLVER_MACORMACK
const NUM_SOLVER: u8 = SOLVER_LAX;
// const int NUM_SOLVER = SOLVER_MACORMACK;

// Problema a resolver: una de las constantes siguientes:
// SHOCK_TUBE1: el tubo de choque clásico de Sod
// SHOCK_TUBE2: choque fuerte: IC con contraste de presión de 10^4
// SHOCK_TUBE3: discontinuidad de contact inmóvil
// SHOCK_TUBE3A: como anterior, pero DC se mueve lentamente
// SHOCK_TUBE4: dos ondas de rarefacción
// SHOCK_TUBE5: dos choques
const PROBLEM: u8 = SHOCK_TUBE_SOD;

// Condiciones de frontera
// BC_FREEFLOW: salida libre (gradiente cero)
// BC_INFLOW: entrada; debe ser especificada en boundary()
// BC_REFLECTIVE: reflectiva
// BC_PERIODIC: periódica
const BC_LEFT: u8 = BC_FREEFLOW;
const BC_RIGHT: u8 = BC_FREEFLOW;

// Directorio donde escribir las salidas (usar "./" para dir actual)
// Debe terminar en una diagonal '/'
const OUT_DIR: &str = "./temp/";

const do_output: bool = false;

/*============================================================================*/
/* NO ES NECESARIO MODIFICAR NADA DEBAJO DE ESTA LÍNEA                        */
/*============================================================================*/

// Constantes calculadas de las anteriores

const DX: f64 = (XR-XL)/(NX as f64);      // Espaciamiento de la malla

/*============================================================================*/
// CONDICIONES INICIALES

// Devuelve la coordenada X de la celda i
fn xcoord(i: usize) -> f64 {
  return XL + (i as f64)*DX;
}

// -----------------------------------------------------------------------------

// Condición inicial para un tubo de choque genéricos
fn shocktube_IC(U: &mut Vec<Vec<f64>>, x0: f64, rhoL: f64, uL: f64, pL: f64, rhoR: f64, uR: f64, pR: f64) {

  let U1L = rhoL;
  let U2L = rhoL*uL;
  let U3L = 0.5*rhoL*uL*uL + pL/(GAMMA-1.0);
  let U1R = rhoR;
  let U2R = rhoR*uR;
  let U3R = 0.5*rhoR*uR*uR + pR/(GAMMA-1.0);
  let mut x;
  let i: i32;

  for i in 0..=NX+1 {
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
fn initflow(U: &mut Vec<Vec<f64>>) {

  // Inicializar variables conservadas U de acuerdo al problema elegido
  if (PROBLEM == SHOCK_TUBE1) {
    // Tubo de choque clásico de Sod
    shocktube_IC(U, 0.5, 1.0, 0.0, 1.0, 0.125, 0.0, 0.1);
  } else if (PROBLEM == SHOCK_TUBE2) {
    // Choque fuerte
    shocktube_IC(U, 0.5, 1.0, 0.0, 1000.0, 1.0, 0.0, 0.1);
  } else if (PROBLEM == SHOCK_TUBE3) {
    // Discontinuidad de contacto estacionaria
    shocktube_IC(U, 0.5, 1.4, 0.0, 1.0, 1.0, 0.0, 1.0);
  } else if (PROBLEM == SHOCK_TUBE3A) {
    // Discontinuidad de contacto que se mueve lentamente
    shocktube_IC(U, 0.5, 1.4, 0.1, 1.0, 1.0, 0.1, 1.0);
  } else if (PROBLEM == SHOCK_TUBE4) {
    // Dos ondas de rarefacción
    shocktube_IC(U, 0.5, 1.0, -1.0, 0.4, 1.0, 1.0, 0.4);
  } else if (PROBLEM == SHOCK_TUBE5) {
    // Dos choques
    shocktube_IC(U, 0.5, 1.0, 1.0, 0.4, 1.0, -1.0, 0.4);
  } else {
    println!("Initial condition not set in initflow\n");
    std::process::exit;
  }

}

// -----------------------------------------------------------------------------

// Inicializaciones generales, como variables globales
fn initmain(t: &mut f64, it: &mut i32, nout: &mut i32, tout: &mut f64) {

  *t = 0.0;
  *it = 0;
  *nout = 0;
  *tout = 0.0;

}

/*============================================================================*/

// Aplica condiciones de frontera a celdas fantasma del arreglo pasado
fn boundary(U: &mut Vec<Vec<f64>>) {

  // BC a la izquierda
  if (BC_LEFT == BC_FREEFLOW) {
    for e in 0..NEQ {
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
    for e in 0..NEQ {
      U[0][e] = U[NX][e];
    }
  }

  // BC a la derecha
  if (BC_RIGHT == BC_FREEFLOW) {
    for e in 0..NEQ {
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
    for e in 0..NEQ {
      U[NX+1][e] = U[0][e];
    }
  }

}

/*============================================================================*/

// Calcula las primitivas a partir de las U pasadas
// Esto incluye las celdas fantasma
fn calc_prims(U: &Vec<Vec<f64>>, PRIM: &mut Vec<Vec<f64>>) {

  for i in 0..=NX+1 {
    PRIM[i][0] = U[i][0];
    PRIM[i][1] = U[i][1]/U[i][0];
    PRIM[i][2] = (GAMMA-1.0)*(U[i][2] - U[i][1]*U[i][1]/(2.0*U[i][0]));
  }

}

/*============================================================================*/

// Calcula la velocidad del sonido dada la presión y densidad
fn sound_speed(pres: f64, rho: f64) -> f64 {
  return f64::sqrt(GAMMA*pres/rho);
}

/*============================================================================*/

// Calcula el paso de tiempo resultante de la condición CFL
fn timestep(P: &Vec<Vec<f64>>) -> f64 {

  // Determinamos la velocidad máxima en la malla
  let mut max_speed: f64 = 0.0;
  for i in 1..=NX {
    let cs = sound_speed(P[i][2], P[i][0]);
    let u_plus_cs = f64::abs(P[i][1]) + cs;
    if (u_plus_cs > max_speed) {
      max_speed = u_plus_cs;
    }
  }

  // Condición de estabilidad CFL
  return CFL * DX / max_speed;

}

/*============================================================================*/

// Calcular los flujos físicos F -- Ecuaciones de Euler 1D
// No olvidar calcular F en las celdas fantasma!
fn fluxes(PRIM: &Vec<Vec<f64>>, F: &mut Vec<Vec<f64>>) {

  for i in 0..=NX+1 {
    F[i][0] = PRIM[i][0] * PRIM[i][1];
    F[i][1] = PRIM[i][0] * PRIM[i][1]*PRIM[i][1] + PRIM[i][2];
    F[i][2] = PRIM[i][1] * (0.5*PRIM[i][0]*PRIM[i][1]*PRIM[i][1] + GAMMA/(GAMMA-1.0) * PRIM[i][2]);
  }

}


/*============================================================================*/
// SOLVERS NUMÉRICOS

// -----------------------------------------------------------------------------

// Método de Lax-Friedrichs
fn Lax(dt: f64, U: &Vec<Vec<f64>>, F: &Vec<Vec<f64>>, UP: &mut Vec<Vec<f64>>) {

  // Aplicar método de Lax
  for i in 1..=NX {
    for e in 0..NEQ {
      UP[i][e] = (U[i+1][e] + U[i-1][e])/2.0 - dt/(2.0*DX) * (F[i+1][e] - F[i-1][e]);
    }
  }

}

// -----------------------------------------------------------------------------

// Método de Macormack
fn Macormack(dt: f64, U: &Vec<Vec<f64>>, F: &mut Vec<Vec<f64>>, PRIM: &mut Vec<Vec<f64>>, UT: &mut Vec<Vec<f64>>, UP: &mut Vec<Vec<f64>>) {

  // Paso predictor: actualizamos las UT con flujos hacia adelante
  for i in 1..=NX {
    for e in 0..NEQ {
      UT[i][e] = U[i][e] - dt/DX * (F[i+1][e] - F[i][e]);
    }
  }

  // Aplicamos las BCs a las UT
  boundary(UT);

  // Actualizar las primitivas de nuevo usando las UT esta vez
  calc_prims(&UT, PRIM);

  // Re-calculamos los flujos F pero usando las primitivas actualizadas
  fluxes(&PRIM, F);

  // Paso corrector: obtenemos las UP usando U, UT y F actualizados
  for i in 1..=NX {
    for e in 0..NEQ {
      UP[i][e] = (U[i][e] + UT[i][e])/2.0 - dt/(2.0*DX) * (F[i][e] - F[i-1][e]);
    }
  }

}

/*============================================================================*/
// STEPPING

// Volcar las UPs sobre las Us "finales" del paso de tiempo -- sin viscosidad
// Incluye las celdas fantasma
fn step_simple(U: &mut Vec<Vec<f64>>, UP: &Vec<Vec<f64>>) {

  // Esto incluye las celdas fantasma
  for i in 1..=NX {
    for e in 0..NEQ {
      U[i][e] = UP[i][e];
    }
  }

}

// -----------------------------------------------------------------------------
// Volcar las UPs sobre las Us "finales" del paso de tiempo, aplicando
// viscosidad artificial donde haya máximos o mínimos locales
// Incluye las celdas fantasma
fn step_viscosity(U: &mut Vec<Vec<f64>>, UP: &Vec<Vec<f64>>) {

  for i in 1..=NX {
    for e in 0..NEQ {
      // Aplicamos la viscosidad sólo donde hay mínimos/máximos locales
      // En las demás celdas simplemente copiamos las UP sobre las U
      // ETA debe ser estrictamente menor que 1/2
      if ((U[i+1][e]-U[i][e])*(U[i][e]-U[i-1][e]) < 0.0) {
        U[i][e] = UP[i][e] + ETA*(UP[i+1][e] + UP[i-1][e] - 2.0*UP[i][e]);
      } else {
        U[i][e] = UP[i][e];
      }
    }
  }

  // Lo de arriba no toca las dos celdas fantasma, pero éstas también
  // deben ser copiadas (sin viscosidad)
  for e in 0..NEQ {
    U[0][e] = UP[0][e];
    U[NX+1][e] = UP[NX+1][e];
  }

}


/*============================================================================*/

// Escribe a disco el estado de la simulación
fn output(PRIM: &Vec<Vec<f64>>, nout: &mut i32, tout: &mut f64, t: f64, it: i32, dt: f64) {

  if (do_output) {

    // Generar el nombre del archivo de salida
    let fname: String = format!("{}/output_{:02}.txt", OUT_DIR, nout);

    // Abrir el archivo
    let mut file = BufWriter::new(File::create(fname).expect("Unable to create file"));

    // Escribir la posición, densidad, velocidad y presión a disco, separadas
    // por espacios

    for i in 1..=NX { 
      let x = xcoord(i);
      write!(file, "{} {} {} {}\n", x, PRIM[i][0], PRIM[i][1], PRIM[i][2]);
    }

    // Cerrar archivo
    drop(file);

    println!("Salida {}, t={:.7}, it={}, dt={:.7}", nout, t, it, dt);

  }

  // Avanzar variables de output
  *nout = *nout + 1;
  *tout = (*nout as f64) * DTOUT;

}

/*============================================================================*/
// PROGRAMA PRINCIPAL

fn main() {

  use std::time::Instant;
  let start = Instant::now();

  // VARIABLES GLOBALES

  let mut U = vec![vec![0.0; NEQ]; NX+2];;
  let mut UP = vec![vec![0.0; NEQ]; NX+2];;
  let mut F = vec![vec![0.0; NEQ]; NX+2];;
  let mut UT = vec![vec![0.0; NEQ]; NX+2];;
  let mut PRIM = vec![vec![0.0; NEQ]; NX+2];;

  let mut dt: f64 = 0.0;         // Paso de tiempo
  let mut t: f64 = 0.0;          // Tiempo actual
  let mut it: i32 = 0;            // Iteración actual
  let mut tout: f64 = 0.0;       // Tiempo de la siguiente salida
  let mut nout: i32 = 0;          // Número de la siguiente salida

  // Inicializaciones generales
  initmain(&mut t, &mut it, &mut nout, &mut tout);

  // Condición inicial
  initflow(&mut U);

  // Llenar celdas fantasma
  boundary(&mut U);

  // Actualizar las primitivas
  calc_prims(&U, &mut PRIM);

  // Escribir condición inicial a disco
  output(&PRIM, &mut nout, &mut tout, t, it, dt);

  // BUCLE PRINCIPAL
  while (t < TFIN) {
    
    // Calcular el paso de tiempo dt
    dt = timestep(&PRIM);

    // Actualizar los flujos F
    fluxes(&PRIM, &mut F);

    // Aplicar el método numérico para calcular las UP
    if (NUM_SOLVER == SOLVER_LAX) {
      Lax(dt, &U, &F, &mut UP);
    } else if (NUM_SOLVER == SOLVER_MACORMACK) {
      Macormack(dt, &U, &mut F, &mut PRIM, &mut UT, &mut UP);
    }

    // Aplicar condiciones de frontera a las UP recién calculadas
    boundary(&mut UP);

    // Avanzar las U, con viscosidad para Macormack
    if (NUM_SOLVER == SOLVER_LAX || ETA == 0.0) {
      step_simple(&mut U, &UP);
    } else if (NUM_SOLVER == SOLVER_MACORMACK) {
      step_viscosity(&mut U, &UP);
    }

    // Actualizar las primitivas PRIM usando las nuevas U
    calc_prims(&U, &mut PRIM);

    // Avanzar el estado de la simulación
    t = t + dt;
    it = it + 1;
    // println!("t={}, t={}, dt={}", it, t, dt);
  
    // Escribir a disco, si aplicable
    if (t >= tout) {
      output(&PRIM, &mut nout, &mut tout, t, it, dt);
    }

  }

  println!("{:.3?}", start.elapsed().as_secs_f64());

}

/*============================================================================*/
