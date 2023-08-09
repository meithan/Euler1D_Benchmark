! ==============================================================================
! Euler1D.f90
!
! por J.C. Toledo-Roy
! 20 / oct / 2020
!
! Este programa resuelve las ecuaciones de Euler de la dinámica de gases en 1D
! usando el método de Lax o el de Macormack.
!
! Incluye algunos tubos de choque, como el de Sod, como pruebas.
!
! ==============================================================================

! ==============================================================================
! CONSTANTES Y VARIABLES GLOBALES
! ==============================================================================
module named_constants
  implicit none
  ! Constantes con nombre (a usar como opciones más abajo) -- No modificar

  ! Solvers numéricos
  integer, parameter :: SOLVER_LAX = 1
  integer, parameter :: SOLVER_MACORMACK = 2

  ! Problemas (ICs) -- ver descripciones en la siguiente sección
  integer, parameter ::  SHOCK_TUBE1 = 1
  integer, parameter ::  SHOCK_TUBE2 = 2
  integer, parameter ::  SHOCK_TUBE3 = 3
  integer, parameter ::  SHOCK_TUBE3A = 4
  integer, parameter ::  SHOCK_TUBE4 = 5
  integer, parameter ::  SHOCK_TUBE5 = 6

  ! Tipos de condiciones de frontera
  integer, parameter ::  BC_FREEFLOW = 1
  integer, parameter ::  BC_INFLOW = 2
  integer, parameter ::  BC_REFLECTIVE = 3
  integer, parameter ::  BC_PERIODIC = 4
end module

module globals
  use named_constants
  implicit none
  ! ----------------------------------------------------------------------------
  ! CONSTANTES Y VARIABLES GLOBALES ESPECIFICADAS POR EL USUARIO

  ! Parámetros de la simulación
  integer, parameter :: NEQ = 3         ! Número de ecuaciones
  integer, parameter :: NX = 5000         ! Tamaño de la malla
  real, parameter :: XL = 0.0         ! Coordenada física del extremo izquierdo
  real, parameter :: XR = 1.0         ! Coordenada física del extremo derecho
  real, parameter :: TFIN = 0.20       ! Tiempo final de integración
  real, parameter :: CFL = 0.9        ! Parametro de Courant
  real, parameter :: DTOUT = TFIN/10     ! Intervalo para escribir a disco
  real, parameter :: GAMMA = 1.4      ! Razón de capacidades caloríficas

  ! Viscosidad artficial, sólo para Macormack
  real, parameter :: ETA = 0.05

  ! Solucionador numérico: SOLVER_LAX o SOLVER_MACORMACK
  integer, parameter :: NUM_SOLVER = SOLVER_LAX
  ! integer, parameter ::  NUM__SOLVER = SOLVER_MACORMACK

  ! Problema a resolver, una de las constantes siguientes:
  ! SHOCK_TUBE1: el tubo de choque clásico de Sod
  ! SHOCK_TUBE2: choque fuerte: IC con contraste de presión de 10^4
  ! SHOCK_TUBE3: discontinuidad de contact inmóvil
  ! SHOCK_TUBE3A: como anterior, pero DC se mueve lentamente
  ! SHOCK_TUBE4: dos ondas de rarefacción
  ! SHOCK_TUBE5: dos choques
  integer, parameter :: PROBLEM = SHOCK_TUBE1

  ! Condiciones de frontera, usar las constantes siguientes:
  ! BC_FREEFLOW: salida libre (gradiente cero)
  ! BC_INFLOW: entrada; debe ser especificada en boundary()
  ! BC_REFLECTIVE: reflectiva
  ! BC_PERIODIC: periódica
  integer, parameter :: BC_LEFT = BC_FREEFLOW
  integer, parameter :: BC_RIGHT = BC_FREEFLOW

  ! Directorio donde escribir las salidas (usar "./" para dir actual)
  ! Debe terminar en una diagonal '/'
  character(len=128), parameter :: OUT_DIR = "./temp/"

  logical, parameter :: do_output = .false.

  ! ============================================================================
  ! NO ES NECESARIO MODIFICAR NADA DEBAJO DE ESTA LÍNEA
  ! ============================================================================

  ! Constantes calculadas de las anteriores
  real, parameter ::  DX = (XR-XL)/NX      ! Espaciamiento de la malla

  ! ----------------------------------------------------------------------------
  ! VARIABLES GLOBALES

  real    U(NEQ,0:NX+1)      ! Variables conservadas actuales
  real   UP(NEQ,0:NX+1)      ! Variables conservadas "avanzadas"
  real    F(NEQ,0:NX+1)      ! Flujos físicos
  real   UT(NEQ,0:NX+1)      ! Us temporales (sólo usadas en Macormack)
  real PRIM(NEQ,0:NX+1)      ! Variables primitivas

  real :: time       ! Tiempo actual
  integer :: it      ! Iteración actual
  real :: dt         ! Paso de tiempo
  integer :: nout    ! Número de la siguiente salida
  real :: tout       ! Tiempo de la siguiente salida

  ! Para medir tiempo de ejecución
  integer :: clock_start, clock_count, clock_rate, clock_max

end module globals

! ==============================================================================
! SUBRUTINAS
! ==============================================================================

! Devuelve la coordenada X de la celda i
subroutine xcoord(i, x)
  
  use globals, only: XL, DX
  implicit none
  integer, intent(in) :: i
  real, intent(out) :: x
  x = XL + i*DX
end subroutine xcoord

! ==============================================================================
! Condición inicial para imponer tubos de choque genéricos
subroutine shocktube_IC(U, x0, rhoL, uL, pL, rhoR, uR, pR)
  
  use globals, only: NEQ, NX, GAMMA
  implicit none
  real, intent(inout) :: U(NEQ, 0:NX+1)
  real, intent(in) :: x0, rhoL, uL, pL, rhoR, uR, pR

  real :: U1L, U2L, U3L, U1R, U2R, U3R
  integer :: i
  real :: x

  U1L = rhoL
  U2L = rhoL*uL
  U3L = 0.5*rhoL*uL*uL + pL/(GAMMA-1)
  U1R = rhoR
  U2R = rhoR*uR
  U3R = 0.5*rhoR*uR*uR + pR/(GAMMA-1)

  do i=0,NX+1
    call xcoord(i, x)
    if (x <= x0) then
      U(1,i) = U1L
      U(2,i) = U2L
      U(3,i) = U3L
    else
      U(1,i) = U1R
      U(2,i) = U2R
      U(3,i) = U3R
    end if
  end do

  ! Para eliminar el artefacto de "escalones" del método de Lax
  ! do e=1,3
  !   U(e,NX/2) = 0.5 * (U(e,NX/2+1) + U(e,NX/2-1))
  ! end do

end subroutine shocktube_IC

! ==============================================================================
! Impone las condiciones iniciales
! Debe inicializar los valores de U en todo el dominio
! Nótese que también debemos llenar las celdas fantasma
subroutine initflow(U)
  
  use named_constants
  use globals, only: NEQ, NX, time, it, tout, nout, PROBLEM
  implicit none
  real, intent(inout) :: U(NEQ, 0:NX+1)

  ! Inicializar variables conservadas U de acuerdo al problema elegido
  if (PROBLEM == SHOCK_TUBE1) then
    ! Tubo de choque clásico de Sod
    call shocktube_IC(U, 0.5, 1.0, 0.0, 1.0, 0.125, 0.0, 0.1)
  
  else if (PROBLEM == SHOCK_TUBE2) then
    ! Choque fuerte
    call shocktube_IC(U, 0.5, 1.0, 0.0, 1000.0, 1.0, 0.0, 0.1)
  
  else if (PROBLEM == SHOCK_TUBE3) then
    ! Discontinuidad de contacto estacionaria
    call shocktube_IC(U, 0.5, 1.4, 0.0, 1.0, 1.0, 0.0, 1.0)
  
  else if (PROBLEM == SHOCK_TUBE3A) then
    ! Discontinuidad de contacto que se mueve lentamente
    call shocktube_IC(U, 0.5, 1.4, 0.1, 1.0, 1.0, 0.1, 1.0)
  
  else if (PROBLEM == SHOCK_TUBE4) then
    ! Dos ondas de rarefacción
    call shocktube_IC(U, 0.5, 1.0, -1.0, 0.4, 1.0, 1.0, 0.4)
  
  else if (PROBLEM == SHOCK_TUBE5) then
    ! Dos choques
    call shocktube_IC(U, 0.5, 1.0, 1.0, 0.4, 1.0, -1.0, 0.4)
  
  else
    print*, "Initial condition not set in initflow\n"
    stop
  
  end if

  ! Inicializar otras variables
  time = 0.0
  it = 0
  tout = 0.0
  nout = 0

end subroutine initflow

! ==============================================================================

! Aplicar condiciones de frontera a celdas fantasma
! El arreglo pasado es al que aplicaremos las BCs
subroutine boundary(U)
  
  use named_constants
  use globals, only: NEQ, NX, BC_LEFT, BC_RIGHT
  implicit none
  real, intent(inout) :: U(NEQ, 0:NX+1)

  integer :: e

  ! BC a la izquierda
  if (BC_LEFT == BC_FREEFLOW) then
    do e=1,NEQ
      U(e,0) = U(e,1)
    end do
  else if (BC_LEFT == BC_INFLOW) then
    U(1,0) = 0.0
    U(2,0) = 0.0
    U(3,0) = 0.0
  else if (BC_LEFT == BC_REFLECTIVE) then
    U(1,0) =  U(1,1)
    U(2,0) = -U(2,1)
    U(3,0) =  U(3,1)
  else if (BC_LEFT == BC_PERIODIC) then
    do e=1,NEQ
      U(e,0) = U(e,NX)
    end do
  end if

  ! BC a la derecha
  if (BC_RIGHT == BC_FREEFLOW) then
    do e=1,NEQ
      U(e,NX+1) = U(e,NX)
    end do
  else if (BC_RIGHT == BC_INFLOW) then
    U(1,NX+1) = 0.0
    U(2,NX+1) = 0.0
    U(3,NX+1) = 0.0
  else if (BC_RIGHT == BC_REFLECTIVE) then
    U(1,NX+1) =  U(1,NX)
    U(2,NX+1) = -U(2,NX)
    U(3,NX+1) =  U(3,NX)
  else if (BC_RIGHT == BC_PERIODIC) then
    do e=1,NEQ
      U(e,NX+1) = U(e,0)
    end do
  end if

end subroutine boundary

! ==============================================================================

! Calcula las primitivas a partir de las U pasadas
! Esto incluye las celdas fantasma
subroutine flow2prim(U, PRIM)
  
  use globals, only: NX, NEQ, GAMMA
  implicit none
  real, intent(in) :: U(NEQ, 0:NX+1)
  real, intent(out) :: PRIM(NEQ, 0:NX+1)

  integer :: i

  do i=0, NX+1
    PRIM(1,i) = U(1,i)
    PRIM(2,i) = U(2,i)/U(1,i)
    PRIM(3,i) = (GAMMA-1)*(U(3,i) - U(2,i)*U(2,i)/(2*U(1,i)))
  end do

end subroutine flow2prim

! ==============================================================================

! Calcular los flujos físicos F -- Ecuaciones de Euler 1D
! No olvidar calcular F en las celdas fantasma!
subroutine fluxes(PRIM, F)
  
  use globals, only: NX, NEQ, GAMMA
  implicit none
  real, intent(in) :: PRIM(NEQ, 0:NX+1)
  real, intent(out) :: F(NEQ, 0:NX+1)

  integer :: i

  do i=0,NX+1
    F(1,i) = PRIM(1,i) * PRIM(2,i)
    F(2,i) = PRIM(1,i) * PRIM(2,i)*PRIM(2,i) + PRIM(3,i)
    F(3,i) = PRIM(2,i)*(0.5*PRIM(1,i)*PRIM(2,i)*PRIM(2,i) + GAMMA/(GAMMA-1)*PRIM(3,i))
  end do

end subroutine fluxes

! ==============================================================================

! Calcula el paso de tiempo resultante de la condición CFL
subroutine timestep(PRIM, dt)
  
  use globals, only: NEQ, NX, GAMMA, CFL, DX
  implicit none
  real, intent(in) :: PRIM(NEQ, 0:NX+1)
  real, intent(out) :: dt

  real :: cs, u_plus_cs, max_speed
  integer :: i

  ! Determinamos la velocidad máxima en la malla
  max_speed = 0.0;
  do i=1,NX
    cs = sqrt(GAMMA*PRIM(3,i)/PRIM(1,i))
    u_plus_cs = abs(PRIM(2,i)) + cs
    if (u_plus_cs > max_speed) then
      max_speed = u_plus_cs
    end if
  end do

  ! Condición de estabilidad CFL
  dt = CFL * DX / max_speed;

end subroutine timestep

! ==============================================================================
! SOLVERS NUMÉRICOS

! Método de Lax-Friedrichs
subroutine Lax(U, F, UP)
  
  use globals, only: NEQ, NX, DX, dt
  implicit none
  real, intent(in) :: U(NEQ, 0:NX+1)
  real, intent(in) :: F(NEQ, 0:NX+1)
  real, intent(out) :: UP(NEQ, 0:NX+1)

  integer :: i, e

  do i=1,NX
    do e=1,NEQ
      UP(e,i) = (U(e,i+1) + U(e,i-1))/2.0 &
                 - dt/(2*DX) * (F(e,i+1) - F(e,i-1))
    end do
  end do

end subroutine Lax

! -----------------------------------------------------------------------------

! Método de Macormack
subroutine Macormack(U, F, UP)
  
  use globals, only: NEQ, NX, DX, dt, UT, PRIM
  implicit none
  real, intent(in) :: U(NEQ, 0:NX+1)
  real, intent(inout) :: F(NEQ, 0:NX+1)
  real, intent(out) :: UP(NEQ, 0:NX+1)

  integer :: i, e

  ! Paso predictor: actualizamos las UT con flujos hacia adelante
  do i=1,NX
    do e=1,NEQ
      UT(e,i) = U(e,i) - dt/DX * (F(e,i+1) - F(e,i))
    end do
  end do

  ! Aplicamos las BCs a las UT
  call boundary(UT)

  ! Actualizar las primitivas de nuevo usando las UT esta vez
  call flow2prim(UT, PRIM)

  ! Re-calculamos los flujos F pero usando las primitivas actualizadas
  call fluxes(PRIM, F)

  ! Paso corrector: obtenemos las UP usando U, UT y F actualizados
  do i=1,NX
    do e=1,NEQ
      UP(e,i) = (U(e,i) + UT(e,i))/2 - dt/(2*DX) * (F(e,i) - F(e,i-1))
    end do
  end do

end subroutine Macormack

! -----------------------------------------------------------------------------
subroutine solver(U, F, UP)
  
  use named_constants
  use globals, only: NEQ, NX, NUM_SOLVER
  implicit none
  real, intent(in) :: U(NEQ, 0:NX+1)
  real, intent(inout) :: F(NEQ, 0:NX+1)
  real, intent(out) :: UP(NEQ, 0:NX+1)

  if (NUM_SOLVER == SOLVER_LAX) then
    call Lax(U, F, UP)
  else if (NUM_SOLVER == SOLVER_MACORMACK) then
    call Macormack(U, F, UP)
  end if

end subroutine solver

! ==============================================================================
! STEPPING

! Volcar las UPs sobre las Us finales del paso de tiempo -- sin viscosidad
subroutine step_simple(U, UP)
  
  use globals, only: NEQ, NX
  implicit none
  real, intent(out) :: U(NEQ, 0:NX+1)
  real, intent(in) :: UP(NEQ, 0:NX+1)

  integer :: i, e

  do i=0,NX+1
    do e=1,NEQ
      U(e,i) = UP(e,i)
    end do
  end do

end subroutine step_simple

! ------------------------------------------------------------------------------

! Volcar las UPs sobre las Us "finales" del paso de tiempo, aplicando
! viscosidad artificial donde haya máximos o mínimos locales
subroutine step_viscosity(U, UP)
  
  use globals, only: NEQ, NX, ETA
  implicit none
  real, intent(out) :: U(NEQ, 0:NX+1)
  real, intent(in) :: UP(NEQ, 0:NX+1)

  integer :: i, e

  do i=1,NX
    do e=1,NEQ
      ! Aplicamos la viscosidad sólo donde hay mínimos/máximos locales
      ! En las demás celdas simplemente copiamos las UP sobre las U
      ! ETA debe ser estrictamente menor que 1/2
      if ((U(e,i+1)-U(e,i))*(U(e,i)-U(e,i-1)) < 0) then
        U(e,i) = UP(e,i) + ETA*(UP(e,i+1) + UP(e,i-1) - 2*UP(e,i))
      else
        U(e,i) = UP(e,i)
      end if
    end do
  end do

  ! Lo de arriba no toca las dos celdas fantasma, pero éstas también
  ! deben ser copiadas (sin viscosidad)
  do e=1,NEQ
    U(e,0) = UP(e,0)
    U(e,NX+1) = UP(e,NX+1)
  end do

end subroutine step_viscosity

! ------------------------------------------------------------------------------
! "Wrapper" para el stepping, dependiendo del solver usado
subroutine step(U, UP)
  
  use named_constants
  use globals, only: NEQ, NX, ETA, NUM_SOLVER
  implicit none
  real, intent(out) :: U(NEQ, 0:NX+1)
  real, intent(in) :: UP(NEQ, 0:NX+1)
  
  if (NUM_SOLVER == SOLVER_LAX) then
    call step_simple(U, UP)
  else if (NUM_SOLVER == SOLVER_MACORMACK) then
    call step_viscosity(U, UP)
  end if

end subroutine step


! ==============================================================================

! Escribe a disco el estado de la simulación
subroutine output(PRIM)
  use globals, only: NX, NEQ, XL, DX, nout, tout, DTOUT, OUT_DIR, do_output
  implicit none
  real, intent(in) :: PRIM(NEQ,0:NX+1)

  character(len=256) :: fname
  real :: x
  integer :: i

  if (do_output) then

    ! Generar el nombre del archivo de salida
    write(fname, "(a,a,i2.2,a)") TRIM(OUT_DIR), "output_", nout, ".txt"

    ! Abrir el archivo
    open(unit=10, file=fname, status='unknown')

    ! Escribir la posición, densidad, velocidad y presión a disco,
    ! separadas por espacios
    do i=1,NX
      x = XL + i*DX
      write(10,*) x, " ", PRIM(1,i), " ", PRIM(2,i), " ", PRIM(3,i)
    end do

    ! Cerrar archivo
    close(10)

    write(*,'(a,a)') "Se escribió ", trim(fname)

  end if

  ! Avanzar variables para el siguiente output
  nout = nout + 1
  tout = nout * DTOUT

end subroutine

! ==============================================================================
! Programa principal
! ==============================================================================

program Euler1D

  use named_constants
  use globals
  implicit none

  call system_clock(clock_start, clock_rate, clock_max)

  ! Condición inicial e inicializaciones
  call initflow(U)

  ! Llenar celdas fantasma
  call boundary(U)

  ! Calcular las primitivas
  call flow2prim(U, PRIM)

  ! Escribir condición inicial a disco
  call output(PRIM)

  ! Bucle principal
  do while (time < TFIN)

    ! Calcular el paso de tiempo
    call timestep(PRIM, dt)

    ! Actualizar flujos físicos
    call fluxes(PRIM, F)

    ! Aplicar el método numérico para calcular las UP
    call solver(U, F, UP)

    ! Aplicar condiciones de frontera a las UP recién calculadas
    call boundary(UP)

    ! Avanzar las U, con viscosidad para Macormack
    call step(U, UP)

    ! Actualizar las primitivas PRIM usando las nuevas U
    call flow2prim(U, PRIM)

    ! Avanzar el estado de la simulación
    time = time + dt
    it = it + 1

    ! Escribir a disco
    if (time >= tout) then
      call output(PRIM)
    end if

  end do

  ! Imprimir tiempo transcurrido
  call system_clock(clock_count, clock_rate, clock_max)
  write(*,'(f8.5)') (clock_count-clock_start)/real(clock_rate)
  !write(*,'(a,i5,a,f8.5,a)') "Se calcularon ", it, " iteraciones en ", (clock_count-clock_start)/real(clock_rate), " s"

end program Euler1D

! ==============================================================================
