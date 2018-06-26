! Run using: gfortran -I/usr/local/include fft.f90 Grid_Definition.f90
! Channel_Solvers_BC.f90 Channel_IC.f90 Channel_Program.f90 -lfftw3 -lm -o channel

PROGRAM Channel_Program
! Program to simulate the scalar transport equation with the assumption
! that velocity is known (so one way coupling), in a channel flow,
! later the momentum equations will be added. It calls the following modules :
! ------------------------------------------------------------------------------
! 1. Fourier_Spectral: To initialise the fft algorithms, generate wavenumbers
!                              in x and z and destroy fft plans
! 1 (a) FFT_Initialise
! Input : Nx, Nz
! Output : plan_fwd, plan_bkd
! 1 (b) create_wavenumbers_1D
! Input : NX, NZ, Lx, Lz
! Output : kx, kz
! 1 (c) FFT_destroy
! Input : plan_fwd, plan_bkd
! ------------------------------------------------------------------------------
! 2. Grid_Definition: To generate the grid based on the number of points
! (Nx, Ny and Nz), length (Lx, Ly and Lz) and stretching factor in y (stretch_y)
! Input : Nx, Ny, Nz, Lx, Ly, Lz, stretch_y
! Output : GX(x base grid), GY(y base grid), GYF(y fractional grid),
!          GZ(z base grid), DY(Delta y base grid), DYF(Delta y fractional grid)
! ------------------------------------------------------------------------------
! 3. Channel_IC : Module containing initial conditions for
!                 velocity, hydrostatic temperature and fluctuation temperature
! 3 (a) Initial_Conditions_velocity
! Input :
! Buck velocity and noise amplitude in velocity: U_bulk and Kick_ini_vel
! Boundary condition type at the walls (Dirichlet (1) and Neumann (2)):
! U_BC_Lower, U_BC_Upper, V_BC_Lower, V_BC_Upper, W_BC_Lower, W_BC_Upper
! - Value of the lower wall velocity (Dirichlet) or velocity gradient (Neumann)-
! U_wall_lower, V_wall_upper, W_wall_lower, U_wall_upper, V_wall_lower,
! W_wall_upper,
! Grid and wavenumber Information: NX, NY, NZ, Lx, Ly, Lz, kx, kz, DY, DYF, GYF
! FFT plans: plan_fwd, plan_bkd
! Output :  Velocities U, V, W
! 3 (b) Initial_Conditions_Background_Temperature
! Input :
! Grid Information: NX, NY, NZ, GX, GYF, GZ, Lx, Ly, Lz
! Boundary condition type: THB_BC_TYPE_X, THB_BC_TYPE_Z (1- Constant,
! 2- Sin variation, 3- Cos variation),
! THB_BC_TYPE_Y (1 - Linear profile and 2 - Tanh profile),
! THB_BC_Lower, THB_BC_Upper (1- Dirichlet, 2- Neumann): type used to generate
! y direction profile
! Boundary condition at walls: THB_wall_lower, THB_wall_upper
! (value (Dirichlet) or gradient (Neumann))
! Hydro_Background: put 0 to remove background temperature
! Output : Background temperature THB
! 3 (c) Initial_Conditions_Temperature
! Input :
! Grid Information: NX, NY, NZ, Ly, GYF
! Initial Condition types: TH_IC_TYPE_Y (initial disturbance about the
! hydrostatic in y direction) (1- no disturbance, 2- sin and 3-cos), Dist_amp
! gives the amplitude of the disturbance and noise amplitude Kick_ini_temp_fluct
! Boundary condition types: TH_BC_Lower, TH_BC_Upper (1- Dirichlet, 2- Neumann)
! Output: Initial temperature fluctuation: TH and
!         Start and End indices to update TH in evolution: K_Start, K_End
! ------------------------------------------------------------------------------
! 4. Channel_Solvers_BC: To get the time step values information and run time
!                        stepping
! 4 (a) courant (get the time step value)
! Input :
! Grid Information: NX, NY, NZ, Lx, Ly, Lz, DYF
! CFL , Richardson, Prandtl, Reynolds number: CFL, Ri, Pr, Re
! Velocities and lower wall background temperature: U, V, W, THB_wall_lower
! Output : time step: delta_t
! 4 (b) RK_SOLVER (to perform the time stepping)
! Input :
! Grid and wavenumber Information: NX, NY, NZ, kx, kz, DY, DYF
! Boundary condition type for scalar: TH_BC_Lower, TH_BC_Upper (1- Dirichlet,
!                                                                   2- Neumann)
! RK algorithm parameters: gamma, zeta, alpha
! time step: delta_t
! fft plans: plan_bkd, plan_fwd
! Richardson, Prandtl, Reynolds number: Ri, Pr, Re
! Velocities and Temperature from previous time steps: U, V, W, TH
! Background temperature: THB
! Start and End indices to update TH in evolution: K_Start, K_End
! Output :
! Temperature fluctuation: TH (later we will also output U, V and W)
! ------------------------------------------------------------------------------
USE Fourier_Spectral
USE Grid_Definition
USE Channel_IC
USE Channel_Solvers_BC
USE, INTRINSIC :: iso_c_binding
IMPLICIT NONE
include 'fftw3.f03'
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14), NX=32, NY=64, NZ=32
INTEGER, PARAMETER :: U_BC_Lower=1, U_BC_Upper=1, V_BC_Lower=1, V_BC_Upper=1, &
                      W_BC_Lower=1, W_BC_Upper=1, THB_BC_TYPE_X=1, &
                      THB_BC_TYPE_Y=1, THB_BC_TYPE_Z=1, THB_BC_Lower=1, &
                      THB_BC_Upper=1, Hydro_Background=1,&
                      TH_IC_TYPE_Y=1, TH_BC_Lower=1, TH_BC_Upper=1
REAL(KIND=DP), PARAMETER :: pi=4.0_DP*ATAN(1.0_DP), U_bulk=1.0_DP, &
  Kick_ini_vel=0.01_DP, Lx=2.0_DP*pi, Ly=5, Lz=2.0_DP*pi, Stretch_y=1.75_DP,&
  n1=0.0_DP, n2=0.0_DP, n3=0.0_DP, Kick_ini_temp_fluct=0.01_DP, &
    Dist_amp=0.1_DP,     U_wall_lower=0.0_DP,    U_wall_upper=0.0_DP, &
    V_wall_lower=0.0_DP, V_wall_upper=0.0_DP,    W_wall_lower=0.0_DP, &
    W_wall_upper=0.0_DP, THB_wall_lower=5.0_DP, THB_wall_upper=10.0_DP, &
    CFL=0.25_DP, Ri=0.25_DP, Pr=0.7_DP, Re=1000.0_DP, Kick_Dist_amp_P=0.0_DP
REAL(KIND=DP) :: delta_t, time, time_final=10.0_DP, T_ref=290.0_DP,&
                 Delta_T_Dim=270.0_DP
REAL(KIND=DP), DIMENSION(1:3) :: gamma, zeta, alpha
COMPLEX(C_DOUBLE_COMPLEX),PARAMETER :: ii=(0.d0, 1.d0)
REAL(KIND=DP), DIMENSION(1:NX) :: GX
REAL(KIND=DP), DIMENSION(0:NY+1) :: GY, GYF
REAL(KIND=DP), DIMENSION(1:NZ) :: GZ, kz
REAL(KIND=DP), DIMENSION(0:NY) :: DY, DYF
REAL(KIND=DP), DIMENSION(1:NX,1:NZ,0:NY+1) :: U, V, W, P, THB, TH, mu, mu_dbl_breve
REAL(KIND=DP), DIMENSION(1:Nx/2+1) :: kx
type(C_PTR) :: plan_fwd, plan_bkd
INTEGER :: I,J,K, K_Start, K_End
! ---------- Find a way to remove these from here and hard code them -----------
gamma(1) = 8.0_DP/15.0_DP
gamma(2) = 5.0_DP/12.0_DP
gamma(3) = 3.0_DP/4.0_DP
zeta(1)  = 0.0_DP
zeta(2)  = -17.0_DP/60.0_DP
zeta(3)  = -5.0_DP/12.0_DP
alpha(1) = 8.0_DP/15.0_DP
alpha(2) = 2.0_DP/15.0_DP
alpha(3) = 1.0_DP/3.0_DP
! --------------- Initialise fft plans and generate wavenumbers ----------------
CALL FFT_Initialise(NX, NZ, plan_fwd, plan_bkd)
CALL create_wavenumbers_1D( NX, NZ, Lx, Lz, kx, kz)
! ------------------------ Generate Grid and Initialise ------------------------
CALL xyz_grid(Lx, Ly, Lz, Stretch_y, NX, NY, NZ, GX, GY, GYF, GZ, DY, DYF)
CALL Initial_Conditions_velocity (U_BC_Lower, U_BC_Upper, V_BC_Lower, &
V_BC_Upper, W_BC_Lower, W_BC_Upper,U_wall_lower, V_wall_lower, W_wall_lower, &
U_wall_upper, V_wall_upper, W_wall_upper, NX, NY, NZ, Lx, Ly, Lz, kx, kz, DY, &
DYF, plan_fwd, plan_bkd, U_bulk, Kick_ini_vel, GYF, U, V, W)
CALL Initial_Conditions_Background_Temperature( NX, NY, NZ, THB_BC_TYPE_X,&
 THB_BC_TYPE_Y, THB_BC_TYPE_Z, THB_BC_Lower, THB_BC_Upper, Hydro_Background, &
 Lx, Ly, Lz, THB_wall_lower, THB_wall_upper, GX, GYF, GZ, THB)
CALL Initial_Conditions_Temperature(NX, NY, NZ, TH_IC_TYPE_Y, &
    TH_BC_Lower, TH_BC_Upper,GYF,Kick_ini_temp_fluct, Dist_amp, K_Start, K_End, Ly, TH)
CALL Viscosity_Temperature(NX, NY, NZ, TH, THB, T_ref,Delta_T_Dim, DY, DYF,K_start, K_end, mu, mu_dbl_breve)
!P=0.0_DP
CALL Initial_Conditions_Pressure(NX, NY, NZ, Kick_Dist_amp_P, P)
! ------------------------------------------------------------------------------
OPEN(unit=1, file='temperature_profile.dat', status='replace')
OPEN(unit=2, file='velocity_profile.dat', status='replace')
OPEN(unit=3, file='y.dat', status='replace')
WRITE(3,*) GYF
time = 0.0_DP
DO WHILE (time .le. time_final)
CALL courant( NX, NY, NZ, Lx, Ly, Lz, DYF, CFL, Ri, Pr, Re, U, V, W, &
                    THB_wall_lower, delta_t )
CALL RK_SOLVER ( K_start, K_end, NX, NY, NZ, TH_BC_Lower, TH_BC_Upper, &
U_BC_Lower,U_BC_Upper, V_BC_Lower,V_BC_Upper, W_BC_Lower, W_BC_Upper, &
kx, kz, gamma, zeta, alpha, delta_t, Lx, Lz, n1, n2, n3, plan_bkd, plan_fwd, DY, DYF, Pr, Re, Ri, &
U, V, W, P, TH, THB,mu, mu_dbl_breve, T_ref,Delta_T_Dim, U_wall_Lower,U_wall_Upper, &
 V_wall_Lower,V_wall_Upper, W_wall_Lower,W_wall_Upper )
  time=time+delta_t
  WRITE(1,*) TH(10,10,:)+THB(10,10,:)
  WRITE(2,*) U(10,10,:)
  print*, delta_t, time, maxval(TH), maxval(U), maxval(V), maxval(W)
END DO
CLOSE(1)
CLOSE(2)
CLOSE(3)
CALL FFT_destroy(plan_fwd, plan_bkd)
END PROGRAM Channel_Program
