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
! THB_BC_Lower, THB_BC_Upper (1- Dirichlet, 2- Neumann):: type used to generate
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
                      TH_IC_TYPE_Y=2, TH_BC_Lower=1, TH_BC_Upper=1, &
                      max_iter=40
REAL(KIND=DP), PARAMETER :: pi=4.0_DP*ATAN(1.0_DP), U_bulk=1.0_DP, &
  Kick_ini_vel=0.0_DP, Lx=2.0_DP*pi, Ly=5, Lz=2.0_DP*pi, Stretch_y=1.75_DP,&
  n1=0.0_DP, n2=0.0_DP, n3=0.0_DP, Kick_ini_temp_fluct=0.0_DP, &
    Dist_amp=1_DP,     U_wall_lower=0.0_DP,    U_wall_upper=0.0_DP, &
    V_wall_lower=0.0_DP, V_wall_upper=0.0_DP,    W_wall_lower=0.0_DP, &
    W_wall_upper=0.0_DP, THB_wall_lower=15.0_DP, THB_wall_upper=10.0_DP, &
    CFL=0.25_DP, Ri=0.25_DP, Pr=7_DP, Re=50.0_DP, Kick_Dist_amp_P=0.0_DP, &
    J, Eo, ETau, vel_norm_sq, TH_norm_sq
REAL(KIND=DP) :: delta_t, time, time_final=10.0_DP, T_ref=290.0_DP,&
                 Delta_T_Dim=270.0_DP
REAL(KIND=DP), DIMENSION(1:3) :: gamma, zeta, alpha
COMPLEX(C_DOUBLE_COMPLEX),PARAMETER :: ii=(0.d0, 1.d0)
REAL(KIND=DP), DIMENSION(1:NX) :: GX
REAL(KIND=DP), DIMENSION(0:NY+1) :: GY, GYF
REAL(KIND=DP), DIMENSION(1:NZ) :: GZ, kz
REAL(KIND=DP), DIMENSION(0:NY) :: DY, DYF
REAL(KIND=DP), DIMENSION(1:NX,1:NZ,0:NY+1) :: U, V, W, P, THB, TH, mu, mu_dbl_breve &
                                              Uo, Vo, Wo, THo, TH_sq
REAL(KIND=DP), DIMENSION(1:Nx/2+1) :: kx
type(C_PTR) :: plan_fwd, plan_bkd
INTEGER :: I,J,K, K_Start, K_End
REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: time_ary, temp_time

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
J_old=0
J=1
DO WHILE (J .gt. J_old)
! --------------------------- Initial Energy Start ---------------------------
CALL Vector_Volume_Integral(U,V,W,U,V,W,Nx,Ny,N_y,Nz,DY, DYF,Lx,Ly,Lz,vel_norm_sq)
TH_sq=TH*TH
Integrate_Volume(TH_sq,Nx,Ny,N_y,Nz,DY_F,Lx,Ly,Lz,TH_norm_sq)
Eo=(1.0_DP/2.0_DP)*(vel_norm_sq+(Ri/(T_ref**2))*TH_norm_sq)
! --------------------------- Initial Energy End ---------------------------
! --------------------------- Direct Solver Start ------------------------------
ALLOCATE(time_ary(1))
time_ary=time

u_fluc=U-U_bar
v_fluc=V-V_bar
w_fluc=W-W_bar
TH_fluc=TH-TH_bar

CALL Dissipation_Calculation( plan_bkd, plan_fwd, kx, kz, DY, DYF, &
 u_fluc, v_fluc, w_fluc, TH_fluc, mu,Dissipation)
 ALLOCATE(diss_ary(1))
 diss_ary=Dissipation

iter=0
! Direct Solution
DO WHILE ((time .lt. time_final) .OR. (iter .ge.2))
iter=iter+1
CALL courant( NX, NY, NZ, Lx, Ly, Lz, DYF, CFL, Ri, Pr, Re, U, V, W, &
                    THB_wall_lower, delta_t )

CALL RK_SOLVER ( K_start, K_end, NX, NY, NZ, TH_BC_Lower, TH_BC_Upper, &
U_BC_Lower,U_BC_Upper, V_BC_Lower,V_BC_Upper, W_BC_Lower, W_BC_Upper, &
kx, kz, gamma, zeta, alpha, delta_t, Lx, Lz, n1, n2, n3, plan_bkd, plan_fwd, DY, DYF, Pr, Re, Ri, &
U, V, W, P, TH, THB,mu, mu_dbl_breve, T_ref,Delta_T_Dim, U_wall_Lower,U_wall_Upper, &
 V_wall_Lower,V_wall_Upper, W_wall_Lower,W_wall_Upper )
  time=time+delta_t

  ALLOCATE(temp_time(size(time_ary)+1 ))
  temp_time(1:size(time_ary))=time_ary
  temp_time(size(time_ary)+1)=time
  DEALLOCATE(time_ary)
  ALLOCATE(time_ary(size(temp_time)))
  time_ary=temp_time
  DEALLOCATE(temp_time)

  u_fluc=U-U_bar
  v_fluc=V-V_bar
  w_fluc=W-W_bar
  TH_fluc=TH-TH_bar

  CALL Dissipation_Calculation( plan_bkd, plan_fwd, kx, kz, DY, DYF, &
   u_fluc, v_fluc, w_fluc, TH_fluc, mu,Dissipation)

   ALLOCATE(temp_diss(size(diss_ary)+1 ))
   temp_diss(1:size(diss_ary))=diss_ary
   temp_diss(size(diss_ary)+1)=Dissipation
   DEALLOCATE(diss_ary)
   ALLOCATE(diss_ary(size(temp_diss)))
   diss_ary=temp_diss
   DEALLOCATE(temp_diss)
END DO
delta_t=time-time_final
CALL RK_SOLVER ( K_start, K_end, NX, NY, NZ, TH_BC_Lower, TH_BC_Upper, &
U_BC_Lower,U_BC_Upper, V_BC_Lower,V_BC_Upper, W_BC_Lower, W_BC_Upper, &
kx, kz, gamma, zeta, alpha, delta_t, Lx, Lz, n1, n2, n3, plan_bkd, plan_fwd, DY, DYF, Pr, Re, Ri, &
U, V, W, P, TH, THB,mu, mu_dbl_breve, T_ref,Delta_T_Dim, U_wall_Lower,U_wall_Upper, &
 V_wall_Lower,V_wall_Upper, W_wall_Lower,W_wall_Upper )


 ALLOCATE(temp_time(size(time_ary)+1 ))
 temp_time(1:size(time_ary))=time_ary
 temp_time(size(time_ary)+1)=time_final
 DEALLOCATE(time_ary)
 ALLOCATE(time_ary(size(temp_time)))
 time_ary=temp_time
 DEALLOCATE(temp_time)


 u_fluc=U-U_bar
 v_fluc=V-V_bar
 w_fluc=W-W_bar
 TH_fluc=TH-TH_bar

 CALL Dissipation_Calculation( plan_bkd, plan_fwd, kx, kz, DY, DYF, &
  u_fluc, v_fluc, w_fluc, TH_fluc, mu,Dissipation)

  ALLOCATE(temp_diss(size(diss_ary)+1 ))
  temp_diss(1:size(diss_ary))=diss_ary
  temp_diss(size(diss_ary)+1)=Dissipation
  DEALLOCATE(diss_ary)
  ALLOCATE(diss_ary(size(temp_diss)))
  diss_ary=temp_diss
  DEALLOCATE(diss_ary)


 ! --------------------------- Initial Energy Start ---------------------------
 CALL Vector_Volume_Integral(U,V,W,U,V,W,Nx,Ny,N_y,Nz,DY, DYF,Lx,Ly,Lz,vel_norm_sq)
 TH_sq=TH*TH
 Integrate_Volume(TH_sq,Nx,Ny,N_y,Nz,DY_F,Lx,Ly,Lz,TH_norm_sq)
 ETau=(1.0_DP/2.0_DP)*(vel_norm_sq+(Ri/(T_ref**2))*TH_norm_sq)
 ! --------------------------- Initial Energy End ---------------------------

J_old=J
 ! ------------------- Obtain The Objective Function: Start ---------------------
 J=A1*ETau/Eo
 A2_diss=0.0_DP
   DO K = 2: size(time_ary)
 J=J+(diss_ary(K-1)+diss_ary(K))/2.0_DP*(time_ary(K)-time_ary(K-1))
 END DO
 ! ------------------- Obtain The Objective Function: End -----------------------

 ! ----------------------------- Direct Solver End -----------------------------

IF (J>J_old)
 ! ------------------- Set Initial Adjoint Fields: Start -----------------------
! Obtain the `initial' condition for the adjoint variables
v1=(A1/Eo* time_final)*U
v2=(A1/Eo* time_final)*V
v3=(A1/Eo* time_final)*W
stau=(A1*Ri* time_final)/(Eo*T_ref**2) * T
CALL Initial_Conditions_Pressure(NX, NY, NZ, Kick_Dist_amp_Q, Q)
! --------------------- Set Initial Adjoint Fields: End ------------------------

! -------------------------- Adjoint Solver Start ------------------------------
delta_t=time_final-time_ary(size(time_ary)-1)
DO I=1:size(time_ary)-1
! Adjoint Solution
SUBROUTINE RK_Solver_Back( K_start, K_end, NX, NY, NZ, v1_BC_Lower, v1_BC_Upper, &
v2_BC_Lower, v2_BC_Upper, v3_BC_Lower, v3_BC_Upper, stau_BC_Lower, stau_BC_Upper, &
U_total_4, V_total_4, W_total_4,  U_bar_4, V_bar_4, W_bar_4, TH_bar_4, TH_total_4, &
THB, DY, DYF, gamma, zeta, alpha, delta_t, Lx, Lz, n1, n2, n3, Pr, Re, Ri, T_ref, &
A2, v1_wall_Lower, v1_wall_Upper, v2_wall_Lower, v2_wall_Upper, v3_wall_Lower, &
v3_wall_Upper, stau_wall_Lower, stau_wall_Upper, plan_bkd, plan_fwd, kx, kz, &
v1, v2, v3, Q, stau )
IF (I .lt. size(time_ary)-1) THEN
delta_t=time_ary(size(time_ary)-I)-time_ary(size(time_ary)-I-1)
END IF
END DO

! -------------------------- Adjoint Solver End ------------------------------

! ------------------- Update Initial Real Field: Start -----------------------
update_initial_real_fields(Nx, Ny, Nz, Lx, Ly, Lz, Ri, T_ref, eps, &
                              DY, DYF,  v1, v2, v3, stau, U, V, W, TH)
! ------------------- Update Initial Real Field: End -----------------------

END IF
IF (iter .ge. max_iter) THEN
STOP
END IF
END DO
CALL FFT_destroy(plan_fwd, plan_bkd)
END PROGRAM Channel_Program
