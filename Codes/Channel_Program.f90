! Run using: gfortran -I/usr/local/include fft.f90 Grid_Definition.f90  Convergence_Check.f90 Channel_Solvers_BC.f90 Adjoint_Solvers.f90 Channel_IC.f90 Channel_Program.f90 -lfftw3 -lm -o channel
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
USE Adjoint_Solvers
USE Convergence_Check
USE, INTRINSIC :: iso_c_binding
IMPLICIT NONE
include 'fftw3.f03'

INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14), NX=32, NY=64, NZ=32
INTEGER, PARAMETER :: U_BC_Lower=1, U_BC_Upper=1, V_BC_Lower=1, V_BC_Upper=1, &
                      W_BC_Lower=1, W_BC_Upper=1, THB_BC_TYPE_X=1, &
                      THB_BC_TYPE_Y=1, THB_BC_TYPE_Z=1, THB_BC_Lower=1, &
                      THB_BC_Upper=1, Hydro_Background=1,&
                      TH_IC_TYPE_Y=2, TH_BC_Lower=1, TH_BC_Upper=1, &
                      max_iter=40, v1_BC_Lower=1, v1_BC_Upper=1, &
                     v2_BC_Lower=1, v2_BC_Upper=1, v3_BC_Lower=1, v3_BC_Upper=1, &
                     stau_BC_Lower=1, stau_BC_Upper=1


REAL(KIND=DP), PARAMETER :: pi=4.0_DP*ATAN(1.0_DP), U_bulk=1.0_DP, &
  Kick_ini_vel=0.1_DP, Lx=2.0_DP*pi, Ly=5, Lz=2.0_DP*pi, Stretch_y=1.75_DP,&
  n1=0.0_DP, n2=0.0_DP, n3=0.0_DP, Kick_ini_temp_fluct=0.0_DP, &
    Dist_amp=1_DP,     U_wall_lower=0.0_DP,    U_wall_upper=0.0_DP, &
    V_wall_lower=0.0_DP, V_wall_upper=0.0_DP,    W_wall_lower=0.0_DP, &
    W_wall_upper=0.0_DP, THB_wall_lower=15.0_DP, THB_wall_upper=10.0_DP, &
    CFL=0.25_DP, Ri=0.25_DP, Pr=7.0_DP, Re=50.0_DP, Kick_Dist_amp_P=0.0_DP, &
    v1_wall_Lower=0.0_DP, v1_wall_Upper=0.0_DP, v2_wall_Lower=0.0_DP, &
    v2_wall_Upper=0.0_DP, v3_wall_Lower=0.0_DP, &
    v3_wall_Upper=0.0_DP, stau_wall_Lower=0.0_DP, stau_wall_Upper=0.0_DP, &
    Kick_Dist_amp_Q=0.0_DP
REAL(KIND=DP) :: delta_t, time, time_final=0.01_DP, T_ref=290.0_DP,&
                 Delta_T_Dim=270.0_DP, A1=0.5_DP, A2=0.5_DP,Dissipation,&
    		 JJ, J_old, Eo, ETau, vel_norm_sq, TH_norm_sq, eps

REAL(KIND=DP), DIMENSION(1:3) :: gamma, zeta, alpha
COMPLEX(C_DOUBLE_COMPLEX),PARAMETER :: ii=(0.d0, 1.d0)
REAL(KIND=DP), DIMENSION(1:NX) :: GX
REAL(KIND=DP), DIMENSION(0:NY+1) :: GY, GYF
REAL(KIND=DP), DIMENSION(1:NZ) :: GZ, kz
REAL(KIND=DP), DIMENSION(0:NY) :: DY, DYF
REAL(KIND=DP), DIMENSION(1:NX,1:NZ,0:NY+1) :: U, V, W, P, THB, TH, mu, mu_dbl_breve, &
                                              Uo, Vo, Wo, THo, TH_sq, v1, v2, v3, Q, stau, &
                                              u_fluc, v_fluc, w_fluc, TH_fluc, &
                                              U_bar, V_bar, W_bar, TH_bar
COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(1:NX/2+1,1:NZ,0:NY+1) ::    F_Vel_func, F_Adj_Vel_func
REAL(KIND=DP), DIMENSION(1:NX,1:NZ,0:NY+1,2) ::    U_total_4_next, V_total_4_next, W_total_4_next, TH_total_4_next
REAL(KIND=DP), DIMENSION(1:NX,1:NZ,0:NY+1,3) ::    U_store,  V_store, W_store, TH_store
REAL(KIND=DP), DIMENSION(1:NX,1:NZ,0:NY+1,4) :: U_total_4, V_total_4, W_total_4,&
                            U_bar_4, V_bar_4, W_bar_4, TH_bar_4, TH_total_4
REAL(KIND=DP), DIMENSION(1:Nx/2+1) :: kx
type(C_PTR) :: plan_fwd, plan_bkd
INTEGER :: I,J,K, K_Start, K_End, iter, iter_big, iter_adj, iter_n, i1, i2, i3, i4
REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: time_ary, temp_time, Eo_check, Eo_check_temp, &
                                        diss_ary, temp_diss

CHARACTER::filename_1*17, filename_2*18, filename_3*23, filename_4*24
CHARACTER(LEN=:), allocatable :: filepath
filepath='/Users/arjunsharma/Documents/NextCloud/Non_Linear_Adjoint_Optimisation/Codes/IO_Files_Store/'
!CHARACTER(LEN=18) ::filename_2
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
    TH_BC_Lower, TH_BC_Upper,GYF,Kick_ini_temp_fluct, Dist_amp, K_Start, K_End,&
    Ly, TH)

CALL Viscosity_Temperature(NX, NY, NZ, TH, THB, T_ref,Delta_T_Dim, DY, DYF, &
K_start, K_end, mu, mu_dbl_breve)

CALL Poisson_RHS(Re, Ri, n1, n2, n3,NX, NY, NZ,  DY, DYF, U, V, W, TH, mu, kx, kz, &
K_start, K_end, plan_bkd, plan_fwd, F_Vel_func)

CALL Initial_Conditions_Pressure(NX, NY, NZ, kx, kz, DY, DYF, plan_bkd, &
                                        F_vel_func, P)


! --------Laminar Flow --------

CALL Initial_Conditions_velocity (U_BC_Lower, U_BC_Upper, V_BC_Lower, &
V_BC_Upper, W_BC_Lower, W_BC_Upper,U_wall_lower, V_wall_lower, W_wall_lower, &
U_wall_upper, V_wall_upper, W_wall_upper, NX, NY, NZ, Lx, Ly, Lz, kx, kz, DY, &
DYF, plan_fwd, plan_bkd, U_bulk, 0.0_DP, GYF, U_bar, V_bar, W_bar)

CALL Initial_Conditions_Temperature(NX, NY, NZ, TH_IC_TYPE_Y, &
    TH_BC_Lower, TH_BC_Upper,GYF,Kick_ini_temp_fluct, 0.0_DP, K_Start, K_End,&
    Ly, TH_bar)
! -------------------
! ------------------------- Obtain Initial Energy --------------------------
CALL Vector_Volume_Integral(U,V,W,U,V,W,Nx,Ny,Nz,DY, DYF,Lx,Ly,Lz, &
                            vel_norm_sq)
TH_sq=TH*TH
CALL Integrate_Volume(TH_sq,Nx,Ny,Ny,Nz,DY,Lx,Ly,Lz,TH_norm_sq)
Eo=(1.0_DP/2.0_DP)*(vel_norm_sq+(Ri/(T_ref**2))*TH_norm_sq)

J_old=0.0_DP
JJ=1.0_DP

iter_big=0


! Biggest Loop to check convergence of Lagrangian
DO WHILE (JJ .gt. J_old)
iter_big=iter_big+1
! Write the initial conditions to appropriate files
write(filename_1,'(A,I3.3,A)') 'U_total_0_',iter_big,'.bin'
open(unit=11,file=filepath//filename_1,form='unformatted',status='replace')
write(11)(((U(i1,i2,i3), &
     i1=1,NX),i2=1,NZ),i3=0,NY+1)
close(11)
write(filename_1,'(A,I3.3,A)') 'V_total_0_',iter_big,'.bin'
open(unit=11,file=filepath//filename_1,form='unformatted',status='replace')
write(11)(((V(i1,i2,i3), &
     i1=1,NX),i2=1,NZ),i3=0,NY+1)
close(11)
write(filename_1,'(A,I3.3,A)') 'W_total_0_',iter_big,'.bin'
open(unit=11,file=filepath//filename_1,form='unformatted',status='replace')
write(11)(((V(i1,i2,i3), &
     i1=1,NX),i2=1,NZ),i3=0,NY+1)
close(11)
write(filename_2,'(A,I3.3,A)') 'TH_total_0_',iter_big,'.bin'
open(unit=11,file=filepath//filename_2,form='unformatted',status='replace')
write(11)(((TH(i1,i2,i3), &
     i1=1,NX),i2=1,NZ),i3=0,NY+1)
close(11)

! --------------------------- Direct Solver Start ------------------------------
ALLOCATE(time_ary(1))
time_ary=time

u_fluc=U-U_bar
v_fluc=V-V_bar
w_fluc=W-W_bar
TH_fluc=TH-TH_bar

CALL Dissipation_Calculation(NX, NY, NZ,Lx, Ly, Lz, Ri, Pr, T_ref, &
plan_bkd, plan_fwd, kx, kz, DY, DYF, u_fluc, v_fluc, w_fluc, TH_fluc, mu, Dissipation)
 ALLOCATE(diss_ary(1))
 diss_ary=Dissipation
print *, Dissipation
iter=0
! Direct Solution
DO WHILE ((time .lt. time_final) .OR. (iter .ge.2)) ! Start RK solver (fwd) loop
iter=iter+1
CALL courant( NX, NY, NZ, Lx, Ly, Lz, DYF, CFL, Ri, Pr, Re, U, V, W, &
                    THB_wall_lower, delta_t )
CALL RK_SOLVER ( K_start, K_end, NX, NY, NZ, TH_BC_Lower, TH_BC_Upper, &
U_BC_Lower,U_BC_Upper, V_BC_Lower,V_BC_Upper, W_BC_Lower, W_BC_Upper, &
kx, kz, gamma, zeta, alpha, delta_t, Lx, Lz, n1, n2, n3, plan_bkd, plan_fwd,&
DY, DYF, Pr, Re, Ri, &
U, V, W, P, TH, THB,mu, mu_dbl_breve, T_ref,Delta_T_Dim, U_wall_Lower, &
 U_wall_Upper,V_wall_Lower,V_wall_Upper, W_wall_Lower,W_wall_Upper,&
 U_store,  V_store, W_store, TH_store  )

 ! writing binary files to store the total U, V, W and TH at each interation
 write(filename_3,'(A,I3.3,A,I7.7,A)') 'U_total_',iter_big,'_',iter,'.bin'
open(unit=11,file=filepath//filename_1,form='unformatted',status='replace')
 write(11)((((U_Store(i1,i2,i3,i4), &
      i1=1,NX),i2=1,NZ),i3=0,NY+1),i4=1,3)
 close(11)
 write(filename_3,'(A,I3.3,A,I7.7,A)') 'V_total_',iter_big,'_',iter,'.bin'
open(unit=11,file=filepath//filename_3,form='unformatted',status='replace')
 write(11)((((V_Store(i1,i2,i3,i4), &
      i1=1,NX),i2=1,NZ),i3=0,NY+1),i4=1,3)
 close(11)
 write(filename_3,'(A,I3.3,A,I7.7,A)') 'W_total_',iter_big,'_',iter,'.bin'
open(unit=11,file=filepath//filename_3,form='unformatted',status='replace')
 write(11)((((W_Store(i1,i2,i3,i4), &
      i1=1,NX),i2=1,NZ),i3=0,NY+1),i4=1,3)
 close(11)
 write(filename_4,'(A,I3.3,A,I7.7,A)') 'TH_total_',iter_big,'_',iter,'.bin'
open(unit=11,file=filepath//filename_4,form='unformatted',status='replace')
 write(11)((((TH_Store(i1,i2,i3,i4), &
      i1=1,NX),i2=1,NZ),i3=0,NY+1),i4=1,3)
 close(11)

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

CALL Dissipation_Calculation(NX, NY, NZ,Lx, Ly, Lz, Ri, Pr, T_ref, &
 plan_bkd, plan_fwd, kx, kz, DY, DYF, u_fluc, v_fluc, w_fluc, TH_fluc, mu,Dissipation)

ALLOCATE(temp_diss(size(diss_ary)+1 ))
temp_diss(1:size(diss_ary))=diss_ary
temp_diss(size(diss_ary)+1)=Dissipation
DEALLOCATE(diss_ary)
ALLOCATE(diss_ary(size(temp_diss)))
diss_ary=temp_diss
DEALLOCATE(temp_diss)
END DO ! Ending RK Solver (fwd) loop

delta_t=time-time_final ! Define the last iteration to finish exactly at t= time_final


CALL RK_SOLVER ( K_start, K_end, NX, NY, NZ, TH_BC_Lower, TH_BC_Upper, &
U_BC_Lower,U_BC_Upper, V_BC_Lower,V_BC_Upper, W_BC_Lower, W_BC_Upper, &
kx, kz, gamma, zeta, alpha, delta_t, Lx, Lz, n1, n2, n3, plan_bkd, plan_fwd, DY, DYF, Pr, Re, Ri, &
U, V, W, P, TH, THB,mu, mu_dbl_breve, T_ref,Delta_T_Dim, U_wall_Lower,U_wall_Upper, &
 V_wall_Lower,V_wall_Upper, W_wall_Lower,W_wall_Upper,U_store,  V_store, W_store, TH_store )

 iter=iter+1

! Replace this with HDF writing
! writing binary files to store the total U, V, W and TH for the last time step of forward RK solver
write(filename_3,'(A,I3.3,A,I7.7,A)') 'U_total_',iter_big,'_',iter,'.bin'
open(unit=11,file=filepath//filename_3,form='unformatted',status='replace')
write(11)((((U_Store(i1,i2,i3,i4), &
     i1=1,NX),i2=1,NZ),i3=0,NY+1),i4=1,3)
close(11)
write(filename_3,'(A,I3.3,A,I7.7,A)') 'V_total_',iter_big,'_',iter,'.bin'
open(unit=11,file=filepath//filename_3,form='unformatted',status='replace')
write(11)((((V_Store(i1,i2,i3,i4), &
     i1=1,NX),i2=1,NZ),i3=0,NY+1),i4=1,3)
close(11)
write(filename_3,'(A,I3.3,A,I7.7,A)') 'W_total_',iter_big,'_',iter,'.bin'
open(unit=11,file=filepath//filename_3,form='unformatted',status='replace')
write(11)((((W_Store(i1,i2,i3,i4), &
     i1=1,NX),i2=1,NZ),i3=0,NY+1),i4=1,3)
close(11)
write(filename_4,'(A,I3.3,A,I7.7,A)') 'TH_total_',iter_big,'_',iter,'.bin'
open(unit=11,file=filepath//filename_3,form='unformatted',status='replace')
write(11)((((TH_Store(i1,i2,i3,i4), &
     i1=1,NX),i2=1,NZ),i3=0,NY+1),i4=1,3)
close(11)

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

 CALL Dissipation_Calculation(NX, NY, NZ,Lx, Ly, Lz, Ri, Pr, T_ref, &
 plan_bkd, plan_fwd, kx, kz, DY, DYF, u_fluc, v_fluc, w_fluc, TH_fluc, mu, Dissipation)

  ALLOCATE(temp_diss(size(diss_ary)+1 ))
  temp_diss(1:size(diss_ary))=diss_ary
  temp_diss(size(diss_ary)+1)=Dissipation
  DEALLOCATE(diss_ary)
  ALLOCATE(diss_ary(size(temp_diss)))
  diss_ary=temp_diss
  DEALLOCATE(temp_diss)
 ! --------------------------- Initial Energy Start ---------------------------
 CALL Vector_Volume_Integral(U,V,W,U,V,W,Nx,Ny,Nz,DY, DYF,Lx,Ly,Lz,vel_norm_sq)
TH_sq=TH*TH
 CALL Integrate_Volume(TH_sq,Nx,Ny,Ny,Nz,DY,Lx,Ly,Lz,TH_norm_sq)
ETau=(1.0_DP/2.0_DP)*(vel_norm_sq+(Ri/(T_ref**2))*TH_norm_sq)
 ! --------------------------- Initial Energy End ---------------------------

J_old=JJ
 ! ------------------- Obtain The Objective Function: Start ---------------------
 JJ=A1*ETau/Eo
 ! Trapezium rule to integrate A2 cost function
DO K = 2, size(time_ary)
    JJ=JJ+A2*(diss_ary(K-1)+diss_ary(K))/2.0_DP*(time_ary(K)-time_ary(K-1))
 END DO
! ------------------- Obtain The Objective Function: End -----------------------

 ! ----------------------------- Direct Solver End -----------------------------

IF (JJ>J_old) THEN
 ! ------------------- Set Initial Adjoint Fields: Start -----------------------
! Obtain the `initial' condition for the adjoint variables
v1=(A1/Eo* time_final)*u_fluc
v2=(A1/Eo* time_final)*v_fluc
v3=(A1/Eo* time_final)*w_fluc
stau=(A1*Ri* time_final)/(Eo*T_ref**2) * TH_fluc


CALL Poisson_Adjoint_RHS (v1,v2,v3,stau, U, V, W,TH,THB, mu, U_bar, V_bar, W_bar, &
kx, kz, NX, NY, NZ, DY, DYF, K_start, K_end, Re, A2, plan_bkd, plan_fwd, F_Adj_Vel_func)

CALL Initial_Conditions_Pressure(NX, NY, NZ, kx, kz, DY, DYF, plan_bkd, &
                                        F_Adj_Vel_func, Q)
! --------------------- Set Initial Adjoint Fields: End ------------------------

! -------------------------- Adjoint Solver Start ------------------------------
delta_t=time_final-time_ary(size(time_ary)-1)
iter_adj=0

! ------------- Reading the real variables for the adjoint solver -------------
write(filename_3,'(A,I3.3,A,I7.7,A)') 'U_total_',iter_big,'_',iter,'.bin'
open(unit=11,file=filepath//filename_3,form='unformatted',status='old')
read(11) ((((U_total_4(i1,i2,i3,4-i4), &
     i1=1,NX),i2=1,NZ),i3=0,NY+1),i4=1,3)
close(11)

write(filename_3,'(A,I3.3,A,I7.7,A)') 'V_total_',iter_big,'_',iter,'.bin'
open(unit=11,file=filepath//filename_3,form='unformatted',status='old')
read(11) ((((V_total_4(i1,i2,i3,4-i4), &
     i1=1,NX),i2=1,NZ),i3=0,NY+1),i4=1,3)
close(11)

write(filename_3,'(A,I3.3,A,I7.7,A)') 'W_total_',iter_big,'_',iter,'.bin'
open(unit=11,file=filepath//filename_3,form='unformatted',status='old')
read(11) ((((W_total_4(i1,i2,i3,4-i4), &
     i1=1,NX),i2=1,NZ),i3=0,NY+1),i4=1,3)
close(11)

write(filename_4,'(A,I3.3,A,I7.7,A)') 'TH_total_',iter_big,'_',iter,'.bin'
open(unit=11,file=filepath//filename_4,form='unformatted',status='old')
read(11) ((((TH_total_4(i1,i2,i3,4-i4), &
     i1=1,NX),i2=1,NZ),i3=0,NY+1),i4=1,3)
close(11)


DO I=1,size(time_ary)-1
! Adjoint Solution
iter_adj=iter_adj+1

IF (iter_adj .gt. 1) THEN
U_total_4(:,:,:,1)=U_total_4(:,:,:,4)
U_total_4(:,:,:,2)=U_total_4_next(:,:,:,2)
U_total_4(:,:,:,3)=U_total_4_next(:,:,:,1)


V_total_4(:,:,:,1)=V_total_4(:,:,:,4)
V_total_4(:,:,:,2)=V_total_4_next(:,:,:,2)
V_total_4(:,:,:,3)=V_total_4_next(:,:,:,1)


W_total_4(:,:,:,1)=W_total_4(:,:,:,4)
W_total_4(:,:,:,2)=W_total_4_next(:,:,:,2)
W_total_4(:,:,:,3)=W_total_4_next(:,:,:,1)


TH_total_4(:,:,:,1)=TH_total_4(:,:,:,4)
TH_total_4(:,:,:,2)=TH_total_4_next(:,:,:,2)
TH_total_4(:,:,:,3)=TH_total_4_next(:,:,:,1)
END IF

IF (iter_adj .lt. iter) THEN
iter_n=iter-1
write(filename_3,'(A,I3.3,A,I7.7,A)') 'U_total_',iter_big,'_',iter_n,'.bin'
open(unit=11,file=filepath//filename_3,form='unformatted',status='old')
read(11) ((((U_total_4_next(i1,i2,i3,i4), &
     i1=1,NX),i2=1,NZ),i3=0,NY+1),i4=1,2)
read(11) (((U_total_4(i1,i2,i3,4), &
     i1=1,NX),i2=1,NZ),i3=0,NY+1)
close(11)


write(filename_3,'(A,I3.3,A,I7.7,A)') 'V_total_',iter_big,'_',iter_n,'.bin'
open(unit=11,file=filepath//filename_3,form='unformatted',status='old')
read(11) ((((V_total_4_next(i1,i2,i3,i4), &
     i1=1,NX),i2=1,NZ),i3=0,NY+1),i4=1,2)
read(11) (((V_total_4(i1,i2,i3,4), &
     i1=1,NX),i2=1,NZ),i3=0,NY+1)
close(11)


write(filename_3,'(A,I3.3,A,I7.7,A)') 'W_total_',iter_big,'_',iter_n,'.bin'
open(unit=11,file=filepath//filename_3,form='unformatted',status='old')
read(11) ((((W_total_4_next(i1,i2,i3,i4), &
     i1=1,NX),i2=1,NZ),i3=0,NY+1),i4=1,2)
read(11) (((W_total_4(i1,i2,i3,4), &
     i1=1,NX),i2=1,NZ),i3=0,NY+1)
close(11)


write(filename_4,'(A,I3.3,A,I7.7,A)') 'TH_total_',iter_big,'_',iter_n,'.bin'
open(unit=11,file=filepath//filename_4,form='unformatted',status='old')
read(11) ((((TH_total_4_next(i1,i2,i3,i4), &
     i1=1,NX),i2=1,NZ),i3=0,NY+1),i4=1,2)
read(11) (((TH_total_4(i1,i2,i3,4), &
     i1=1,NX),i2=1,NZ),i3=0,NY+1)
close(11)

ELSE

write(filename_1,'(A,I3.3,A)') 'U_total_0_',iter_big,'.bin'
open(unit=11,file=filepath//filename_1,form='unformatted',status='old')
read(11) (((U_total_4(i1,i2,i3,4), &
     i1=1,NX),i2=1,NZ),i3=0,NY+1)
close(11)

write(filename_1,'(A,I3.3,A)') 'V_total_0_',iter_big,'.bin'
open(unit=11,file=filepath//filename_1,form='unformatted',status='old')
read(11) (((V_total_4(i1,i2,i3,4), &
     i1=1,NX),i2=1,NZ),i3=0,NY+1)
close(11)


write(filename_1,'(A,I3.3,A)') 'W_total_0_',iter_big,'.bin'
open(unit=11,file=filepath//filename_1,form='unformatted',status='old')
read(11) (((W_total_4(i1,i2,i3,4), &
     i1=1,NX),i2=1,NZ),i3=0,NY+1)
close(11)


write(filename_2,'(A,I3.3,A)') 'TH_total_0_',iter_big,'.bin'
open(unit=11,file=filepath//filename_2,form='unformatted',status='old')
read(11) (((TH_total_4(i1,i2,i3,4), &
     i1=1,NX),i2=1,NZ),i3=0,NY+1)
close(11)

END IF

DO K=1,4
 U_bar_4(:,:,:,K)=U_bar
 V_bar_4(:,:,:,K)=V_bar
 W_bar_4(:,:,:,K)=W_bar
TH_bar_4(:,:,:,K)=TH_bar
END DO

CALL RK_Solver_Back( K_start, K_end, NX, NY, NZ, v1_BC_Lower, v1_BC_Upper, &
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
CALL update_initial_real_fields(Nx, Ny, Nz, Lx, Ly, Lz, Ri, T_ref, eps, &
                              DY, DYF,  v1, v2, v3, stau, U, V, W, TH)
! ------------------- Update Initial Real Field: End -----------------------

END IF
IF (iter .ge. max_iter) THEN
STOP
END IF
END DO
CALL FFT_destroy(plan_fwd, plan_bkd)
END PROGRAM Channel_Program
