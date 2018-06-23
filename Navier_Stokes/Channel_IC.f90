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
V_BC_Upper, W_BC_Lower, W_BC_Upper,U_wall_lower, V_wall_upper, W_wall_lower, &
U_wall_upper, V_wall_lower, W_wall_upper, NX, NY, NZ, DY, DYF, U, V, W)
CALL Remove_Divergence(NX, NY, NZ, Lx, Lz, alfa_t, kx, kz, DY, DYF, plan_fwd, &
                              plan_bkd, U, V, W, P)
CALL Velocity_IC_Boundary_Conditions(U_BC_Lower, U_BC_Upper, V_BC_Lower, &
V_BC_Upper, W_BC_Lower, W_BC_Upper,U_wall_lower, V_wall_upper, W_wall_lower, &
U_wall_upper, V_wall_lower, W_wall_upper, NX, NY, NZ, DY, DYF, U, V, W)
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


END MODULE Channel_IC
