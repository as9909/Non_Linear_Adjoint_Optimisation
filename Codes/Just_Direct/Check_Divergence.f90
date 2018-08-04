! Run using: gfortran -I/usr/local/include fft.f90 Grid_Definition.f90  Convergence_Check.f90 Channel_Solvers_BC.f90 Channel_IC.f90 Check_Divergence.f90 -lfftw3 -lm -o channel
PROGRAM Check_Divergence
USE Fourier_Spectral
USE Grid_Definition
USE Channel_IC
USE Channel_Solvers_BC
USE Convergence_Check
USE, INTRINSIC :: iso_c_binding
IMPLICIT NONE
include 'fftw3.f03'
include 'fftw3l.f03'
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14), NX=16, NY=32, NZ=32
INTEGER, PARAMETER :: U_BC_Lower=1, U_BC_Upper=1, V_BC_Lower=1, V_BC_Upper=1, &
                      W_BC_Lower=1, W_BC_Upper=1, THB_BC_TYPE_X=1, &
                      THB_BC_TYPE_Y=1, THB_BC_TYPE_Z=1, THB_BC_Lower=1, &
                      THB_BC_Upper=1, Hydro_Background=1,&
                      TH_IC_TYPE_Y=2, TH_BC_Lower=1, TH_BC_Upper=1, &
                      max_iter=40, v1_BC_Lower=1, v1_BC_Upper=1, &
                     v2_BC_Lower=1, v2_BC_Upper=1, v3_BC_Lower=1, v3_BC_Upper=1, &
                     stau_BC_Lower=1, stau_BC_Upper=1

! pi=4.0_DP*ATAN(1.0_DP),
REAL(KIND=DP), PARAMETER ::  pi=4.0_DP*ATAN(1.0_DP), &
 U_bulk=1.0_DP, &
  Kick_ini_vel=0.01_DP, Lx=4.0_DP*pi, Ly=2.0_DP, Lz=2.0_DP*pi, Stretch_y=0.0000005_DP,&
  n1=0.0_DP, n2=1.0_DP, n3=0.0_DP, Kick_ini_temp_fluct=0.00_DP, &
    Dist_amp=0.0_DP,     U_wall_lower=0.0_DP,    U_wall_upper=0.0_DP, &
    V_wall_lower=0.0_DP, V_wall_upper=0.0_DP,    W_wall_lower=0.0_DP, &
    W_wall_upper=0.0_DP, THB_wall_lower=0.0_DP, THB_wall_upper=0.0_DP, &
    CFL=0.5_DP, Ri=0.0_DP, Pr=1.0_DP, Re=100.0_DP, Kick_Dist_amp_P=0.0_DP


REAL(KIND=DP) :: delta_t, time, time_final=100.0_DP, T_ref=290.0_DP,&
                 Delta_T_Dim=THB_wall_upper-THB_wall_lower, Drive_x, Drive_y, Drive_z

REAL(KIND=DP), DIMENSION(1:3) :: gamma, zeta, alpha

REAL(KIND=DP), DIMENSION(1:NX) :: GX
REAL(KIND=DP), DIMENSION(0:NY+1) :: GY, GYF
REAL(KIND=DP), DIMENSION(1:NZ) :: GZ, kz
REAL(KIND=DP), DIMENSION(0:NY) :: DY, DYF
REAL(KIND=DP), DIMENSION(1:NX,1:NZ,0:NY+1) :: U, V, W, P, THB, TH, mu, mu_dbl_breve, &
                                              Uo, Vo, Wo, THo, TH_sq, &
                                              U_bar, V_bar, W_bar, TH_bar, chk_Poss_RHS
COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(1:NX/2+1,1:NZ,0:NY+1) ::    F_Vel_func, F_Adj_Vel_func
REAL(KIND=DP), DIMENSION(1:NX,1:NZ,0:NY+1,2) ::    U_total_4_next, V_total_4_next, W_total_4_next, TH_total_4_next
REAL(KIND=DP), DIMENSION(1:NX,1:NZ,0:NY+1,3) ::    U_store,  V_store, W_store, TH_store, A_st
REAL(KIND=DP), DIMENSION(1:NX,1:NZ,0:NY+1,4) :: U_total_4, V_total_4, W_total_4,&
                            U_bar_4, V_bar_4, W_bar_4, TH_bar_4, TH_total_4
REAL(KIND=DP), DIMENSION(1:Nx/2+1) :: kx
type(C_PTR) :: plan_fwd, plan_bkd
INTEGER :: I,J,K, K_Start, K_End, iter, iter_big, iter_adj, iter_n, i1, i2, i3, i4, getcwd, status
REAL(KIND=DP), DIMENSION(:), ALLOCATABLE :: time_ary, temp_time, Eo_check, Eo_check_temp, &
                                        diss_ary, temp_diss

CHARACTER::filename_1*17, filename_2*18, filename_3*23, filename_4*24
CHARACTER(LEN=:), allocatable :: filepath
!status = getcwd(filepath)
!write(*,*) filepath
!filepath = trim(filepath) // 'IO_Files_Store'
!write(*,*) filepath
!filepath='/home/arjun.sharma/Channel_Flow/IO_Files_Store/'
!filepath='/home/ritabrata.thakur/Desktop/Just_Direct/IO_Files_Store/'
filepath='/Users/arjunsharma/Documents/NextCloud/Non_Linear_Adjoint_Optimisation/Codes/Just_Direct/IO_Files_Store/'

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

Drive_z=(12.0_DP/Re)*(U_bulk/Ly**2) !2.0_DP/Re!
Drive_y=0.0_DP
Drive_X=0.0_DP
! --------------- Initialise fft plans and generate wavenumbers ----------------
CALL FFT_Initialise(NX, NZ, plan_fwd, plan_bkd)
CALL create_wavenumbers_1D( NX, NZ, Lx, Lz, kx, kz)

! ------------------------ Generate Grid and Initialise ------------------------
CALL xyz_grid(Lx, Ly, Lz, Stretch_y, NX, NY, NZ, GX, GY, GYF, GZ, DY, DYF)

CALL Initial_Conditions_velocity (U_BC_Lower, U_BC_Upper, V_BC_Lower, &
V_BC_Upper, W_BC_Lower, W_BC_Upper,U_wall_lower, V_wall_lower, W_wall_lower, &
U_wall_upper, V_wall_upper, W_wall_upper, NX, NY, NZ, Lx, Ly, Lz, kx, kz, DY, &
DYF, plan_fwd, plan_bkd, U_bulk, Kick_ini_vel, GYF, U, V, W)

CALL FFT_destroy(plan_fwd, plan_bkd)
END PROGRAM Check_Divergence
