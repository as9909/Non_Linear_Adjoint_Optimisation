MODULE Channel_IC
! --------- Module for the channel flow initial and boundary conditions --------
! ---------------- It consists of the following subroutines: -------------------
! 1. Initial_Conditions_velocity: Initial velocity generation
! 2. Initial_Conditions_Background_Temperature: Initial background temperature
! 3. Initial_Conditions_Temperature: Initial fluctuating temperature
CONTAINS
! ------------------ 1. Subroutine for initialising velocity -------------------
SUBROUTINE Initial_Conditions_velocity( U_BC_Lower, U_BC_Upper, V_BC_Lower, &
V_BC_Upper, W_BC_Lower, W_BC_Upper, U_wall_lower, V_wall_upper, W_wall_lower, &
U_wall_upper, V_wall_lower, W_wall_upper, NX, NY, NZ, Lx, Ly, Lz, kx, kz, DY, &
DYF, plan_fwd, plan_bkd, U_bulk, Kick, GYF, U, V, W)
USE, INTRINSIC :: iso_c_binding
USE Fourier_Spectral
USE Channel_Solvers_BC
IMPLICIT NONE
include 'fftw3.f03'
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: U_BC_Lower, U_BC_Upper, V_BC_Lower, V_BC_Upper, &
                       W_BC_Lower, W_BC_Upper
REAL(KIND=DP), INTENT(IN) :: U_wall_lower, V_wall_upper, W_wall_lower, &
                             U_wall_upper, V_wall_lower, W_wall_upper
INTEGER, INTENT(IN) :: NX, NY, NZ
REAL(KIND=DP), INTENT(IN) :: U_bulk, Kick, Lx, Ly, Lz
INTEGER :: I, J, K
INTEGER, DIMENSION(:), ALLOCATABLE :: seed
REAL(KIND=DP) :: Rand_num
REAL(KIND=DP), PARAMETER :: alfa_t=1
REAL(KIND=DP), DIMENSION(1:NZ), INTENT(IN) :: kz
REAL(KIND=DP), DIMENSION(1:NX/2+1), INTENT(IN) :: kx
REAL(KIND=DP), DIMENSION(0:NY), INTENT(IN) :: DY, DYF
type(C_PTR), INTENT(IN) :: plan_fwd, plan_bkd
REAL(KIND=DP), DIMENSION(0:NY+1), INTENT(IN)  :: GYF
REAL(KIND=DP), DIMENSION(1:NX,1:NZ,0:NY+1), INTENT(OUT) :: U, V, W
REAL(KIND=DP), DIMENSION(1:NX,1:NZ,0:NY+1) :: P

P = 0.0_DP
DO I  = 1,NX
  DO J = 1, NZ
    DO K = 0, NY+1
U(I,J,K)=3.0_DP/2.0_DP*U_bulk*(1.0_DP-(2*GYF(K)/Ly)**2)
END DO
END DO
END DO
V=0.0_DP
W=0.0_DP
CALL RANDOM_SEED
DO I  = 1,NX
  DO J = 1, NZ
    DO K = 0, NY+1
      CALL RANDOM_NUMBER(Rand_num)
      U(I,J,K) = U(I,J,K) + Kick*(Rand_num-0.5_DP)
      CALL RANDOM_NUMBER(Rand_num)
      V(I,J,K) = V(I,J,K) + Kick*(Rand_num-0.5_DP)
      CALL RANDOM_NUMBER(Rand_num)
      W(I,J,K) = W(I,J,K) + Kick*(Rand_num-0.5_DP)
    END DO
  END DO
END DO
U(:,:,0)=0.0_DP
U(:,:,NY+1)=0.0_DP
V(:,:,0)=0.0_DP
V(:,:,NY+1)=0.0_DP
W(:,:,0)=0.0_DP
W(:,:,NY+1)=0.0_DP
CALL Velocity_IC_Boundary_Conditions(U_BC_Lower, U_BC_Upper, V_BC_Lower, &
  V_BC_Upper, W_BC_Lower, W_BC_Upper,U_wall_Lower, V_wall_Lower, W_wall_Lower, &
   U_wall_Upper, V_wall_Upper, W_wall_Upper, NX, NY, NZ, DY, DYF, U, V, W)
CALL Remove_Divergence(NX, NY, NZ, Lx, Lz, alfa_t, kx, kz, DY, DYF, plan_fwd, &
                              plan_bkd, U, V, W, P)
CALL Velocity_IC_Boundary_Conditions(U_BC_Lower, U_BC_Upper, V_BC_Lower, &
V_BC_Upper, W_BC_Lower, W_BC_Upper,U_wall_lower, V_wall_Lower, W_wall_lower, &
U_wall_upper, V_wall_Upper, W_wall_upper, NX, NY, NZ, DY, DYF, U, V, W)
END SUBROUTINE Initial_Conditions_velocity

! ------------ 2. Subroutine for initialising background temperature -----------
SUBROUTINE Initial_Conditions_Background_Temperature( NX, NY, NZ, THB_BC_TYPE_X,&
 THB_BC_TYPE_Y, THB_BC_TYPE_Z, THB_BC_Lower, THB_BC_Upper, Hydro_Background, &
 Lx, Ly, Lz, THB_wall_lower, THB_wall_upper, GX, GYF, GZ, THB)
IMPLICIT NONE
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
REAL(KIND=DP), PARAMETER :: pi=4.0_DP*ATAN(1.0_DP)
INTEGER, INTENT(IN) :: NX, NY, NZ, THB_BC_TYPE_X, THB_BC_TYPE_Y, THB_BC_TYPE_Z,&
                       THB_BC_Lower, THB_BC_Upper, Hydro_Background
REAL(KIND=DP), INTENT(IN) :: Lx, Ly, Lz, THB_wall_lower, THB_wall_upper
REAL(KIND=DP), DIMENSION(1:NX, 1:NZ, 0:NY+1), INTENT(OUT) :: THB
REAL(KIND=DP), DIMENSION(1:NX/2+1), INTENT(IN)  :: GX
REAL(KIND=DP), DIMENSION(0:NY+1), INTENT(IN)  :: GYF
REAL(KIND=DP), DIMENSION(1:NZ), INTENT(IN)  :: GZ
INTEGER :: I, J, K
IF (Hydro_Background .EQ. 0) THEN
  THB(1:NX,0:NY+1,1:NZ)=0.0_DP
ELSE
! Boundary conditions in y
DO I=1,NX
  DO J=1,NZ
    DO K=0,NY+1
  IF (THB_BC_TYPE_Y .EQ. 1) THEN ! Linear Profile
    IF ((THB_BC_Lower .EQ. 1) .AND. (THB_BC_Upper .EQ. 1)) THEN
    ! Dirichlet on both walls
      THB(I,J,K) = (THB_wall_upper-THB_wall_lower)/LY*(GYF(K)-GYF(1))&
                    +THB_wall_lower
    ELSEIF ((THB_BC_Lower .EQ. 2) .AND. (THB_BC_Upper .EQ. 1)) THEN
    ! Dirichlet on lower wall and Neumann on upper
      THB(I,J,K) = THB_wall_upper*(GYF(K)-GYF(1))+THB_wall_lower
    ELSEIF ((THB_BC_Lower .EQ. 1) .AND. (THB_BC_Upper .EQ. 2)) THEN
    ! Dirichlet on upper wall and Neumann on bottom
      THB(I,J,K) = THB_wall_lower*(GYF(K)-GYF(NY))+THB_wall_upper
    ELSEIF ((THB_BC_Lower .EQ. 2) .AND. (THB_BC_Upper .EQ. 2)) THEN
    ! Neumann on both walls
      STOP "Background temperature boundary conditions entered in y are not &
       appropriate, can't have Neumann on both walls"
    ENDIF
  ELSEIF (THB_BC_TYPE_Y .EQ. 2) THEN ! Tanh profile
    THB(I,J,K) = 0.5_DP*(THB_wall_upper-THB_wall_lower)*TANH(GYF(K))
  ELSE
  STOP 'Background temperature boundary conditions entered in y are not &
   appropriate, enter 1 for linear or 2 for tanh in y'
  ENDIF
END DO
END DO
END DO
FORALL (I=1:NX, J=1:NZ, K=0:NY+1)
  THB(I,J,K) = THB(I,J,K) - THB(I,J,1) &
                              +  THB_wall_lower
END FORALL
! Boundary conditions in x
DO I=1,NX
  DO J=1,NZ
    DO K=0,NY+1
IF (THB_BC_TYPE_X .EQ. 1) THEN ! Constant Profile
  THB(I,J,K) = THB(I,J,K) * 1.0_DP
ELSEIF (THB_BC_TYPE_X .EQ. 2) THEN ! Sin x
  THB(I,J,K) = THB(I,J,K) * SIN(2.0_DP*pi*GX(I)/Lx)
ELSEIF (THB_BC_TYPE_X .EQ. 3) THEN ! Cos x
  THB(I,J,K) = THB(I,J,K) * COS(2.0_DP*pi*GX(I)/Lx)
ENDIF
END DO
END DO
END DO
! Boundary conditions in z
DO I=1,NX
  DO J=1,NZ
    DO K=0,NY+1
IF (THB_BC_TYPE_Z .EQ. 1) THEN ! Constant Profile
  THB(I,J,K) = THB(I,J,K) * 1.0_DP
ELSEIF (THB_BC_TYPE_Z .EQ. 2) THEN ! Sin x
  THB(I,J,K) = THB(I,J,K) * SIN(2.0_DP*pi*GZ(J)/Lz)
ELSEIF (THB_BC_TYPE_Z .EQ. 3) THEN ! Cos x
  THB(I,J,K) = THB(I,J,K) * COS(2.0_DP*pi*GZ(J)/Lz)
ENDIF
END DO
END DO
END DO
ENDIF
END SUBROUTINE Initial_Conditions_Background_Temperature

! ----------- 3. Subroutine for initialising fluctuating temperature -----------
SUBROUTINE Initial_Conditions_Temperature(NX, NY, NZ, TH_IC_TYPE_Y,&
    TH_BC_Lower, TH_BC_Upper,GYF,Kick, Dist_amp, K_Start, K_End, Ly, TH)
IMPLICIT NONE
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER :: I, J, K
INTEGER, DIMENSION(:), ALLOCATABLE :: seed
REAL(KIND=DP), PARAMETER :: pi=4.0_DP*ATAN(1.0_DP)
REAL(KIND=DP) :: Rand_num
REAL(KIND=DP), INTENT(IN) :: Kick, Dist_amp, Ly
INTEGER, INTENT(OUT) :: K_Start, K_End
INTEGER, INTENT(IN) :: NX, NY, NZ, TH_IC_TYPE_Y,&
                       TH_BC_Lower, TH_BC_Upper
REAL(KIND=DP), DIMENSION(0:NY+1), INTENT(IN)  :: GYF
REAL(KIND=DP), DIMENSION(1:NX,1:NZ,0:NY+1), INTENT(OUT) :: TH
CALL RANDOM_SEED
DO I  = 1,NX
  DO J = 1, NZ
    DO K = 0, NY+1
    IF (TH_IC_TYPE_Y .EQ. 1) THEN
      CALL RANDOM_NUMBER(Rand_num)
      TH(I,J,K) = Kick*(Rand_num-0.5_DP)
    ELSEIF (TH_IC_TYPE_Y .EQ. 2) THEN ! Sinusoidal in y
      CALL RANDOM_NUMBER(Rand_num)
      TH(I,J,K) = Kick*(Rand_num-0.5_DP)+Dist_amp*SIN(2.0_DP*pi*GYF(K)/Ly)
    ELSEIF (TH_IC_TYPE_Y .EQ. 3) THEN ! Sinusoidal (cos) in y
      CALL RANDOM_NUMBER(Rand_num)
      TH(I,J,K) = Kick*(Rand_num-0.5_DP)+Dist_amp*COS(2.0_DP*pi*GYF(K)/Ly)
    ENDIF
      ! Later we can add distrubances with x and z variations
    END DO
  END DO
END DO
! Boundary conditions in y
DO I=1,NX
DO J=1,NZ
DO K=0,NY+1
  IF ((TH_BC_Lower .EQ. 1) .AND. (TH_BC_Upper .EQ. 1)) THEN
    ! Dirichlet on both walls
    K_Start=2
    K_End=NY-1
    TH(I,J,1) = 0
    TH(I,J,0) = 0
    TH(I,J,NY) = 0
    TH(I,J,NY+1) = 0
  ELSEIF ((TH_BC_Lower .EQ. 2) .AND. (TH_BC_Upper .EQ. 1)) THEN
    ! Dirichlet on lower wall and Neumann on upper
    K_Start=2
    K_End=NY
    TH(I,J,0) = 0
    TH(I,J,1) = 0
    TH(I,J,NY) = TH(I,J,NY-1)
    TH(I,J,NY+1) = TH(I,J,NY)
  ELSEIF ((TH_BC_Lower .EQ. 1) .AND. (TH_BC_Upper .EQ. 2)) THEN
    ! Dirichlet on upper wall and Neumann on bottom
    K_Start=1
    K_End=NY-1
    TH(I,J,1) = TH(I,J,2)
    TH(I,J,0) = TH(I,J,1)
    TH(I,J,NY) = 0
    TH(I,J,NY+1) = 0
  ELSEIF ((TH_BC_Lower .EQ. 2) .AND. (TH_BC_Upper .EQ. 2)) THEN
    K_Start=1
    K_End=NY
    TH(I,J,1) = TH(I,J,2)
    TH(I,J,0) = TH(I,J,1)
    TH(I,J,NY) = TH(I,J,NY-1)
    TH(I,J,NY+1) = TH(I,J,NY)
  ENDIF
END DO
END DO
END DO
END SUBROUTINE Initial_Conditions_Temperature

! ----------- 3. Subroutine for initialising fluctuating temperature -----------
SUBROUTINE Initial_Conditions_Pressure(NX, NY, NZ, kx, kz, DY, DYF, plan_bkd, &
                                        F_vel_funct, P)
! Use Pressure Poisson to obtain initial pressure
USE, INTRINSIC :: iso_c_binding
USE Fourier_Spectral
USE Channel_Solvers_BC
IMPLICIT NONE
include 'fftw3.f03'
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: NX, NY, NZ
REAL(KIND=DP), DIMENSION(1:NX,1:NZ,0:NY+1), INTENT(OUT) :: P
COMPLEX(KIND=DP), DIMENSION(1:NX/2+1,1:NZ,0:NY+1), INTENT(INOUT) :: F_vel_funct
REAL(KIND=DP), DIMENSION(1:NX/2+1), INTENT(IN) :: kx
REAL(KIND=DP), DIMENSION(1:NZ), INTENT(IN) :: kz
REAL(KIND=DP), DIMENSION(0:NY), INTENT(IN) :: DY, DYF
type(C_PTR), INTENT(IN) :: plan_bkd
INTEGER :: I, J, K
REAL(KIND=DP), DIMENSION(1:NX/2+1,1:NZ,0:NY+1) :: A, B, C
! Initialise the matrix and vector
FORALL (I=1:NX/2+1, J=1:NZ, K=0:NY+1)
A(I,J,K) = 0.0_DP
B(I,J,K) = 1.0_DP
C(I,J,K) = 0.0_DP
END FORALL
DO J = 1, NZ
  DO I = 1, NX/2+1
    DO K = 1, NY
      A(I,J,K) =   1.0_DP/(DYF(K)*DY(K-1))
      C(I,J,K) =   1.0_DP/(DYF(K)*DY(K))
      B(I,J,K) = - A(I,J,K) - C(I,J,K) - kx(I)**2 - kz(J)**2
    END DO
  END DO
  ! Boundary Conditions at the walls
  ! 1. p=0: kx=kz=0 i.e. the mean values are set to zero at the bottom wall
  ! 2. dp/dz=0: derivative at other wave-numbers and at top wall is zero
DO I=1,NX/2+1
  IF ((I.EQ.1) .AND. (J.EQ.1)) THEN
    A(I,J,1)=0.0_DP
    B(I,J,1)=1.0_DP
    C(I,J,1)=0.0_DP
    F_vel_funct(I,J,1)=(0.0_DP,0.0_DP)
    A(I,J,NY)=0.0_DP
    B(I,J,NY)=-1.0_DP
    C(I,J,NY)=1.0_DP
    F_vel_funct(I,J,NY)=(0.0_DP,0.0_DP)
 ELSE
  A(I,J,1)=0.0_DP
  B(I,J,1)=1.0_DP
  C(I,J,1)=-1.0_DP
  F_vel_funct(I,J,1)=(0.0_DP,0.0_DP)
  A(I,J,NY)=1.0_DP
  B(I,J,NY)=-1.0_DP
  C(I,J,NY)=0.0_DP
  F_vel_funct(I,J,NY)=(0.0_DP,0.0_DP)
 END IF
END DO
END DO
  CALL Thomas_Matrix_Algorithm_fft(A,B,C,F_vel_funct,NX,NY,NZ)
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_vel_funct, P)
END SUBROUTINE Initial_Conditions_Pressure


END MODULE Channel_IC



SUBROUTINE Poisson_RHS(Re, Ri, n1, n2, n3,NX, NY, NZ,  DY, DYF, U, V, W, TH, mu, kx, kz, &
K_start, K_end, plan_bkd, plan_fwd, F_Vel_func)
  USE, INTRINSIC :: iso_c_binding
  USE Fourier_Spectral
  IMPLICIT NONE
  include 'fftw3.f03'
  INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
  INTEGER, INTENT(IN) :: K_Start, K_End
  INTEGER, INTENT(IN) :: NX, NY, NZ
  REAL(KIND=DP), INTENT(IN) :: Re, Ri, n1, n2, n3
  REAL(KIND=DP), DIMENSION(1:Nx/2+1), INTENT(IN) :: kx
  REAL(KIND=DP), DIMENSION(1:Nz), INTENT(IN) :: kz
  REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1), INTENT(INOUT) :: U, V, W, TH,mu
  REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1) :: Vel_func
  COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(1:NX/2+1,1:NZ,0:NY+1), INTENT(OUT) :: F_Vel_func
REAL(KIND=DP), DIMENSION(0:NY), INTENT(IN) :: DY, DYF
  REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1) :: U_dbl_breve,W_dbl_breve, &
   TH_dbl_breve, mu_dbl_breve, V_bar, Ux, Uz, V_barx, V_barz, Wx, Wz, THx, &
   THz, mux, muz, Uxx, Uzz, V_barxx, V_barzz, Wxx, Wzz, muxx, muzz, muxz, &
   mux_dbl_breve, muz_dbl_breve, V_bary, Uy, Wy, THy, muy, muxy, muzy, Uyy, &
   V_baryy, Wyy, muyy

  COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(1:NX/2+1,1:NZ,0:NY+1) :: F_U, F_V_bar, F_W, F_TH, F_mu,&
  F_Ux, F_Uz, &
  F_V_barx, F_V_barz, F_Wx, F_Wz, F_THx, F_THz, F_mux, F_muz,&
  F_Uxx,F_Uzz,F_V_barxx,F_V_barzz,F_Wxx,F_Wzz,F_muxx, F_muzz, F_muxz

  COMPLEX(C_DOUBLE_COMPLEX), PARAMETER :: ii=(0.0_DP, 1.0_DP)
  type(C_PTR), INTENT(IN) :: plan_bkd, plan_fwd
  INTEGER :: I, J, K

  CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, U, F_U )
  CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, W, F_W )
  CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, mu, F_mu )
  CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, TH, F_TH )
  V_bar=0.0_DP
  FORALL (I=1:NX/2+1,J=1:NZ,K=0:NY)
    V_bar(I,J,K)=(V(I,J,K+1)+V(I,J,K))/2.0_DP
  END FORALL
  CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, V_bar, F_V_bar )
  FORALL  (I=1:NX, J=1:NZ, K=K_start:K_end)
    U_dbl_breve(I,J,K)=(U(I,J,K)*DYF(K-1)+U(I,J,K-1)*DYF(K))/(2.0_DP*DY(K-1))
    W_dbl_breve(I,J,K)=(W(I,J,K)*DYF(K-1)+W(I,J,K-1)*DYF(K))/(2.0_DP*DY(K-1))
    TH_dbl_breve(I,J,K)=(TH(I,J,K)*DYF(K-1)+TH(I,J,K-1)*DYF(K))/(2.0_DP*DY(K-1))
    mu_dbl_breve(I,J,K)=(mu(I,J,K)*DYF(K-1)+mu(I,J,K-1)*DYF(K))/(2.0_DP*DY(K-1))
  END FORALL

  FORALL (I=1:NX/2+1,J=1:NZ,K=0:NY+1)
    F_Ux(I,J,K)=ii*kx(I)*F_U(I,J,K)
    F_Uz(I,J,K)=ii*kz(J)*F_U(I,J,K)
    F_V_barx(I,J,K)=ii*kx(I)*F_V_bar(I,J,K)
    F_V_barz(I,J,K)=ii*kz(J)*F_V_bar(I,J,K)
    F_Wx(I,J,K)=ii*kx(I)*F_W(I,J,K)
    F_Wz(I,J,K)=ii*kz(J)*F_W(I,J,K)
    F_THx(I,J,K)=ii*kx(I)*F_TH(I,J,K)
    F_THz(I,J,K)=ii*kz(J)*F_TH(I,J,K)
    F_mux(I,J,K)=ii*kx(I)*F_mu(I,J,K)
    F_muz(I,J,K)=ii*kz(J)*F_mu(I,J,K)


    F_Uxx(I,J,K)=(-1.0_DP)*(kx(I)**2)*F_U(I,J,K)
    F_Uzz(I,J,K)=(-1.0_DP)*(kz(J)**2)*F_U(I,J,K)
    F_V_barxx(I,J,K)=(-1.0_DP)*(kx(I)**2)*F_V_bar(I,J,K)
    F_V_barzz(I,J,K)=(-1.0_DP)*(kz(J)**2)*F_V_bar(I,J,K)
    F_Wxx(I,J,K)=(-1.0_DP)*(kx(I)**2)*F_W(I,J,K)
    F_Wzz(I,J,K)=(-1.0_DP)*(kz(J)**2)*F_W(I,J,K)
    F_muxx(I,J,K)=(-1.0_DP)*(kx(I)**2)*F_mu(I,J,K)
    F_muzz(I,J,K)=(-1.0_DP)*(kz(J)**2)*F_mu(I,J,K)
    F_muxz(I,J,K)=(-1.0_DP)*(kx(I)*kz(J))*F_mu(I,J,K)
  END FORALL

  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Ux, Ux )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Uz, Uz )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_V_barx, V_barx )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_V_barz, V_barz )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Wx, Wx )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Wz, Wz )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_THx, THx )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_THz, THz )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_mux, mux )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_muz, muz )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Uxx, Uxx )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Uzz, Uzz )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_V_barxx, V_barxx )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_V_barzz, V_barzz )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Wxx, Wxx )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Wzz, Wzz )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_muxx, muxx )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_muzz, muzz )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_muxz, muxz )

  FORALL  (I=1:NX, J=1:NZ, K=K_start:K_end)
    mux_dbl_breve(I,J,K)=(mux(I,J,K)*DYF(K-1)+mux(I,J,K-1)*DYF(K))/(2.0_DP*DY(K-1))
    muz_dbl_breve(I,J,K)=(muz(I,J,K)*DYF(K-1)+muz(I,J,K-1)*DYF(K))/(2.0_DP*DY(K-1))
  END FORALL

  FORALL  (I=1:NX, J=1:NZ, K=K_start:K_end)
    V_bary(I,J,K)=(V(I,J,K+1)-V(I,J,K))/DYF(K)
    Uy(I,J,K)=(U_dbl_breve(I,J,K+1)-U_dbl_breve(I,J,K))/DYF(K)
    Wy(I,J,K)=(W_dbl_breve(I,J,K+1)-W_dbl_breve(I,J,K))/DYF(K)
    THy(I,J,K)=(TH_dbl_breve(I,J,K+1)-TH_dbl_breve(I,J,K))/DYF(K)
    muy(I,J,K)=(mu_dbl_breve(I,J,K+1)-mu_dbl_breve(I,J,K))/DYF(K)
    muxy(I,J,K)=(mux_dbl_breve(I,J,K+1)-mux_dbl_breve(I,J,K))/DYF(K)
    muzy(I,J,K)=(muz_dbl_breve(I,J,K+1)-muz_dbl_breve(I,J,K))/DYF(K)

    Uyy(I,J,K)= ((U(I,J,K+1)-U(I,J,K))/DY(K)-(U(I,J,K)-U(I,J,K-1))/DY(K-1))/DYF(K)
    V_baryy(I,J,K)= ((V_bar(I,J,K+1)-V_bar(I,J,K))/DY(K)-(V_bar(I,J,K)-V_bar(I,J,K-1))/DY(K-1))/DYF(K)
    Wyy(I,J,K)= ((W(I,J,K+1)-W(I,J,K))/DY(K)-(W(I,J,K)-W(I,J,K-1))/DY(K-1))/DYF(K)
    muyy(I,J,K)= ((mu(I,J,K+1)-mu(I,J,K))/DY(K)-(mu(I,J,K)-mu(I,J,K-1))/DY(K-1))/DYF(K)
  END FORALL

Vel_func=(-1.0_DP)*(Ux**2+V_bary**2+Wz**2+2*(Uy*V_barx+Uz*Wx+V_barz*Wy))

Vel_func=Vel_func+Ri*(THx*n1+THy*n2+THz*n3)

Vel_func=Vel_func+ (2.0_DP/Re)*(&
muxx*Ux+muyy*V_bary+muzz*Wz+muxy*(Uy+V_barx)+muxz*(Uz+Wx)+muzy*(Wy+V_barz)+&
mux*(Uxx+Uyy+Uzz)+muy*(V_barxx+V_baryy+V_barzz)+muz*(Wxx+Wyy+Wzz))

CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, Vel_func, F_Vel_func )

END SUBROUTINE Poisson_RHS
