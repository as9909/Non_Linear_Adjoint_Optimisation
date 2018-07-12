MODULE Channel_Solvers_BC
! -------- Module comprising various solvers and boundary conditions -----------
! ---------------- It consists of the following subroutines: -------------------
! 1. Remove_Divergence: Removes divergence from initial noisy velocity (it may
! be used to remove divergence fractionally when momentum equations are added
! later)
! 2. Thomas_Matrix_Algorithm_fft: Implements Thomas algorithm for a complex RHS
! 3. Thomas_Matrix_Algorithm_real: Implements Thomas algorithm for a real RHS
! 4. Velocity_Boundary_Conditions: Boundary conditions for velocity
! 5. Scalar_Boundary_Conditions: Boundary conditions for fluctuating temperature
! 6. courant: gets the time step value
! 7. RK_SOLVER: performs the time stepping
CONTAINS
SUBROUTINE Remove_Divergence(NX, NY, NZ, Lx, Lz, alfa_t, kx, kz, DY, DYF, &
                             plan_fwd, plan_bkd, U, V, W, P)
! -------------- This is used to remove divergence from velocity ---------------
USE, INTRINSIC :: iso_c_binding
USE Fourier_Spectral
IMPLICIT NONE
include 'fftw3.f03'
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
COMPLEX(KIND=DP),PARAMETER :: ii=(0.d0, 1.d0)
INTEGER, INTENT(IN) :: NX, NY, NZ
INTEGER :: I, J, K
COMPLEX(KIND=DP), DIMENSION(1:NX/2+1,1:NZ,0:NY+1) :: F_U, F_V, F_W, F_P, Q
REAL(KIND=DP), DIMENSION(1:NX/2,1:NZ,0:NY+1), INTENT(INOUT) :: U, V, W, P
REAL(KIND=DP), DIMENSION(1:NX/2+1,1:NZ,0:NY+1) :: A, B, C
!COMPLEX(KIND=DP), DIMENSION(1:NX/2+1,0:NY+1):: Vel_Div
REAL(KIND=DP), INTENT(IN) :: Lx, Lz, alfa_t
REAL(KIND=DP), DIMENSION(1:NX/2+1), INTENT(IN) :: kx
REAL(KIND=DP), DIMENSION(1:NZ), INTENT(IN) :: kz
REAL(KIND=DP), DIMENSION(0:NY), INTENT(IN) :: DY, DYF
type(C_PTR), INTENT(IN) :: plan_fwd, plan_bkd
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, U, F_U)
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, V, F_V)
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, W, F_W)
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, P, F_P)
! Initialise the matrix and vector
FORALL (I=1:NX/2+1, J=1:NZ, K=0:NY+1)
A(I,J,K) = 0.0_DP
B(I,J,K) = 1.0_DP
C(I,J,K) = 0.0_DP
Q(I,J,K) = (0.0_DP,0.0_DP)
END FORALL
DO J = 1, NZ
  DO I = 1, NX/2+1
    DO K = 1, NY
      A(I,J,K) = 1/(DYF(K)*DY(K-1))
      C(I,J,K) = 1/(DYF(K)*DY(K))
      B(I,J,K) = - A(I,J,K) - C(I,J,K) - kx(I)**2 - kz(J)**2
      Q(I,J,K) = (ii*kx(I)*F_U(I,J,K)+(F_V(I,J,K+1)-F_V(I,J,K))/DYF(K) &
                     +ii*kz(J)*F_W(I,J,K))/alfa_t
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
    Q(I,J,1)=(0.0_DP,0.0_DP)
    A(I,J,NY)=0.0_DP
    B(I,J,NY)=-1.0_DP
    C(I,J,NY)=1.0_DP
    Q(I,J,NY)=(0.0_DP,0.0_DP)
 ELSE
  A(I,J,1)=0.0_DP
  B(I,J,1)=1.0_DP
  C(I,J,1)=-1.0_DP
  Q(I,J,1)=(0.0_DP,0.0_DP)
  A(I,J,NY)=1.0_DP
  B(I,J,NY)=-1.0_DP
  C(I,J,NY)=0.0_DP
  Q(I,J,NY)=(0.0_DP,0.0_DP)
 END IF
END DO
END DO
  CALL Thomas_Matrix_Algorithm_fft(A,B,C,Q,NX,NY,NZ)
! Remove Divergence from the original velocity field
DO J = 1, NZ
  DO I = 1, NX/2+1
    DO K = 1, NY
      F_U(I,J,K) = F_U(I,J,K) - ii * kx(I) * Q(I,J,K) * alfa_t
      IF (K .ge. 2) THEN ! Ignore the ghost cells
        F_V(I,J,K) = F_V(I,J,K) - (Q(I,J,K)-Q(I,J,K-1))/DY(K-1) * alfa_t
      END IF
      F_W(I,J,K) = F_W(I,J,K) - ii*kz(J)*Q(I,J,K) * alfa_t
      F_P(I,J,K) = F_P(I,J,K) + Q(I,J,K)
    END DO
  END DO
END DO
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_U, U)
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_V, V)
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_W, W)
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_P, P)
END SUBROUTINE Remove_Divergence

SUBROUTINE Thomas_Matrix_Algorithm_fft(A,B,C,D,NX,NY,NZ)
! -------- This is used to implement thomas algorithm for complex RHS ----------
IMPLICIT NONE
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: NX, NY, NZ
REAL(KIND=8), INTENT(INOUT), DIMENSION(1:NX/2+1,1:NZ,0:NY+1) :: A, B, C
COMPLEX(KIND=8), INTENT(INOUT), DIMENSION(1:NX/2+1,1:NZ,0:NY+1) :: D
INTEGER :: I, J, K
DO I = 1,NX/2+1
DO J = 1,NZ
  C(I,J,0) = C(I,J,0)/B(I,J,0)
  D(I,J,0) = D(I,J,0)/B(I,J,0)
  DO K = 1,NY+1 ! Don't Ignore the ghost cells
    C(I,J,K) =  C(I,J,K)/(B(I,J,K)-A(I,J,K)*C(I,J,K-1))
    D(I,J,K) =  (D(I,J,K)-A(I,J,K)*D(I,J,K-1))/(B(I,J,K)-A(I,J,K)*C(I,J,K-1))
  END DO
  DO K = NY,0,-1 ! Don't Ignore the ghost cells
    D(I,J,K) = D(I,J,K)-C(I,J,K)*D(I,J,K+1)
  END DO
  END DO
  END DO
END SUBROUTINE Thomas_Matrix_Algorithm_fft
SUBROUTINE Thomas_Matrix_Algorithm_real(A,B,C,D,NX,NY,NZ)
! -------- This is used to implement thomas algorithm for real RHS ----------
IMPLICIT NONE
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: NX, NY, NZ
REAL(KIND=8), INTENT(INOUT), DIMENSION(1:NX,1:NZ,0:NY+1) :: A, B, C, D
INTEGER :: I, J, K
DO I = 1,NX
DO J = 1,NZ
  C(I,J,0) = C(I,J,0)/B(I,J,0)
  D(I,J,0) = D(I,J,0)/B(I,J,0)
  DO K = 1,NY+1 ! Don't Ignore the ghost cells
    C(I,J,K) =  C(I,J,K)/(B(I,J,K)-A(I,J,K)*C(I,J,K-1))
    D(I,J,K) =  (D(I,J,K)-A(I,J,K)*D(I,J,K-1))/(B(I,J,K)-A(I,J,K)*C(I,J,K-1))
  END DO
  DO K = NY,0,-1 ! Don't Ignore the ghost cells
    D(I,J,K) = D(I,J,K)-C(I,J,K)*D(I,J,K+1)
  END DO
END DO
END DO
END SUBROUTINE Thomas_Matrix_Algorithm_real

! --------------- 4. Subroutine for velocity boundary conditions ---------------
SUBROUTINE Velocity_IC_Boundary_Conditions(U_BC_Lower, U_BC_Upper, V_BC_Lower, &
  V_BC_Upper, W_BC_Lower, W_BC_Upper,U_wall_Lower, V_wall_Lower, W_wall_Lower, &
   U_wall_Upper, V_wall_Upper, W_wall_Upper, NX, NY, NZ, DY, DYF, U, V, W)
IMPLICIT NONE
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: U_BC_Lower, U_BC_Upper, V_BC_Lower, V_BC_Upper, &
                       W_BC_Lower, W_BC_Upper
REAL(KIND=DP), INTENT(IN) :: U_wall_Lower, V_wall_Upper, W_wall_Lower, &
                             U_wall_Upper, V_wall_Lower, W_wall_Upper
REAL(KIND=DP), DIMENSION(1:NX,1:NZ,0:NY+1), INTENT(INOUT) :: U, V, W
INTEGER, INTENT(IN) :: NX, NY, NZ
INTEGER :: I, J
REAL(KIND=DP), DIMENSION(0:NY), INTENT(IN) :: DY, DYF
! --------------------- Set Lower Wall Boundary Conditions ---------------------
IF (U_BC_Lower .EQ. 1) THEN ! Dirichlet Conditions
  DO I = 1,NX
    DO J = 1, NZ
      U(I,J,0)  = U_wall_Lower
      U(I,J,1)  = U_wall_Lower
      END DO
    END DO
ELSE IF (U_BC_Lower .EQ. 2) THEN ! Neumann Conditions
  DO I = 1,NX
    DO J = 1, NZ
      U(I,J,1)  = U(I,J,2)-DY(1)*U_wall_Lower
      U(I,J,0)  = U(I,J,1)-DY(0)*U_wall_Lower
    END DO
  END DO
ELSE
  STOP 'Boundary Conditions entered are not approriate'
END  IF
IF (W_BC_Lower .EQ. 1) THEN ! Dirichlet Conditions
  DO I = 1,NX
    DO J = 1, NZ
      W(I,J,0)  = W_wall_Lower
      W(I,J,1)  = W_wall_Lower
      END DO
    END DO
ELSE IF (W_BC_Lower .EQ. 2) THEN ! Neumann Conditions
  DO I = 1,NX
    DO J = 1, NZ
      W(I,J,1)  = W(I,J,2)-DY(1)*W_wall_Lower
      W(I,J,0)  = W(I,J,1)-DY(0)*W_wall_Lower
    END DO
  END DO
ELSE
  STOP 'Boundary Conditions entered are not approriate'
END  IF
IF (V_BC_Lower .EQ. 1) THEN ! Dirichlet Conditions
  DO I = 1,NX
    DO J = 1, NZ
      V(I,J,1)  = 2.0_DP*V_wall_Lower-V(I,J,2)
      V(I,J,0)  = V(I,J,1)
    END DO
  END DO
ELSE IF (V_BC_Lower .EQ. 2) THEN ! Neumann Conditions
  DO I = 1,NX
    DO J = 1, NZ
      V(I,J,1)  = V(I,J,2)-DYF(1)*V_wall_Lower
      V(I,J,0)  = V(I,J,1)-DYF(0)*V_wall_Lower
    END DO
  END DO
ELSE
  STOP 'Boundary Conditions entered are not approriate'
END IF
! --------------------- Set Upper Wall Boundary Conditions ---------------------
IF (U_BC_Upper .EQ. 1) THEN ! Dirichlet Conditions
  DO I = 1,NX
    DO J = 1, NZ
      U(I,J,NY)   = U_wall_Upper
      U(I,J,NY+1) = U_wall_Upper
      END DO
    END DO
ELSE IF (U_BC_Upper .EQ. 2) THEN ! Neumann Conditions
  DO I = 1,NX
    DO J = 1, NZ
      U(I,J,NY)   = U(I,J,NY-1)+DY(NY-1)*U_wall_Upper
      U(I,J,NY+1) = U(I,J,NY)+DY(NY)*U_wall_Upper
    END DO
  END DO
ELSE
  STOP 'Boundary Conditions entered are not approriate'
END  IF
IF (W_BC_Upper .EQ. 1) THEN ! Dirichlet Conditions
  DO I = 1,NX
    DO J = 1, NZ
      W(I,J,NY)   =  W_wall_Upper
      W(I,J,NY+1) =  W_wall_Upper
      END DO
    END DO
ELSE IF (W_BC_Upper .EQ. 2) THEN ! Neumann Conditions
  DO I = 1,NX
    DO J = 1, NZ
      W(I,J,NY)   = W(I,J,NY-1)+DY(NY-1)*W_wall_Upper
      W(I,J,NY+1) = W(I,J,NY)+DY(NY)*W_wall_Upper
    END DO
  END DO
ELSE
  STOP 'Boundary Conditions entered are not approriate'
END  IF
IF (V_BC_Upper .EQ. 1) THEN ! Dirichlet Conditions
  DO I = 1,NX
    DO J = 1, NZ
      V(I,J,NY+1) = 2.0_DP*V_wall_Upper-V(I,J,NY)
    END DO
  END DO
ELSE IF (V_BC_Upper .EQ. 2) THEN ! Neumann Conditions
  DO I = 1,NX
    DO J = 1, NZ
      V(I,J,NY+1) = V(I,J,NY)+DYF(NY)*V_wall_Upper
    END DO
  END DO
ELSE
  STOP 'Boundary Conditions entered are not approriate'
END IF
END SUBROUTINE Velocity_IC_Boundary_Conditions

! ---------------- 5. Subroutine for scalar boundary conditions ----------------
SUBROUTINE Scalar_Boundary_Conditions(A,B,C,D,NX,NY,NZ,TH_BC_Lower,TH_BC_Upper)
IMPLICIT NONE
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: NX, NY, NZ, TH_BC_Lower,TH_BC_Upper
REAL(KIND=DP), INTENT(INOUT), DIMENSION(1:NX,1:NZ,0:NY+1) :: A, B, C, D
INTEGER :: I,J
DO J = 1, NZ
  DO I = 1, NX
IF (TH_BC_Lower .EQ. 1) THEN
! Dirichlet
A(I,J,0)=0
A(I,J,1)=0
B(I,J,0)=1
B(I,J,1)=1
C(I,J,0)=0
C(I,J,1)=0
D(I,J,0)=0 ! Drichlet will follow from the initial condition of background
D(I,J,1)=0 ! Drichlet will follow from the initial condition of background
ELSE IF (TH_BC_Lower .EQ. 2) THEN
! Neumann
A(I,J,0)=0
A(I,J,1)=0
B(I,J,0)=1
B(I,J,1)=1  !Y(1)
C(I,J,0)=0
C(I,J,1)=-1 ! Y(2)
D(I,J,0)=0 ! Neumann will follow from the initial condition of background
D(I,J,1)=0 ! (Y(2)-Y(1)=0)
END IF
IF (TH_BC_Upper .EQ. 1) THEN
! Dirichlet
A(I,J,NY+1)=0
A(I,J,NY)=0
B(I,J,NY+1)=1
B(I,J,NY)=1
C(I,J,NY+1)=0
C(I,J,NY)=0
D(I,J,NY+1)=0 ! Drichlet will follow from the initial condition of background
D(I,J,NY)=0 ! Drichlet will follow from the initial condition of background
ELSE IF (TH_BC_Upper .EQ. 2) THEN
! Neumann
A(I,J,NY+1)=0
A(I,J,NY)=-1 ! Y(NY-1)
B(I,J,NY+1)=1
B(I,J,NY)=1 ! Y(NY)
C(I,J,NY+1)=0
C(I,J,NY)=0
D(I,J,NY+1)=0 ! Neumann will follow from the initial condition of background
D(I,J,NY)=0 ! (Y(NY)-Y(NY-1)=0)
END IF
END DO
END DO
END SUBROUTINE Scalar_Boundary_Conditions


! ---------------- 5. Subroutine for velocity boundary conditions ----------------
SUBROUTINE V_Boundary_Conditions(A,B,C,D,NX,NY,NZ,DYF,V_BC_Lower,V_BC_Upper,V_Lower,V_Upper)
IMPLICIT NONE
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: NX, NY, NZ, V_BC_Lower,V_BC_Upper
REAL(KIND=DP), INTENT(IN) :: V_Lower,V_Upper
REAL(KIND=DP), INTENT(IN), DIMENSION(0:NY) :: DYF
REAL(KIND=DP), INTENT(INOUT), DIMENSION(1:NX,1:NZ,0:NY+1) :: A, B, C, D
INTEGER :: I,J
DO J = 1, NZ
  DO I = 1, NX
IF (V_BC_Lower .EQ. 1) THEN
! Dirichlet
A(I,J,0)=0
A(I,J,1)=0
B(I,J,0)=1
B(I,J,1)=1
C(I,J,0)=-1
C(I,J,1)=1
D(I,J,0)=0
D(I,J,1)=2*V_Lower
ELSE IF (V_BC_Lower .EQ. 2) THEN
! Neumann
A(I,J,0)=0
A(I,J,1)=0
B(I,J,0)=-1
B(I,J,1)=-1  !Y(1)
C(I,J,0)=1
C(I,J,1)=1 ! Y(2)
D(I,J,0)=DYF(0)*V_Lower ! (Y(1)-Y(0)=V_Lower*DYF(0))
D(I,J,1)=DYF(1)*V_Lower ! (Y(2)-Y(1)=V_Lower*DYF(1))
END IF
IF (V_BC_Upper .EQ. 1) THEN
! Dirichlet
A(I,J,NY+1)=1
B(I,J,NY+1)=1
C(I,J,NY+1)=0
D(I,J,NY+1)=2*V_Upper
ELSE IF (V_BC_Upper .EQ. 2) THEN
! Neumann
A(I,J,NY+1)=-1
B(I,J,NY+1)=1
C(I,J,NY+1)=0
D(I,J,NY+1)=DYF(NY)*V_Upper
END IF
END DO
END DO
END SUBROUTINE V_Boundary_Conditions

SUBROUTINE UW_Boundary_Conditions(A,B,C,D,NX,NY,NZ,DY,UW_vel_BC_Lower,UW_vel_BC_Upper,UW_vel_Lower,UW_vel_Upper)

IMPLICIT NONE
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: NX, NY, NZ, UW_vel_BC_Lower,UW_vel_BC_Upper
REAL(KIND=DP), INTENT(IN) :: UW_vel_Lower,UW_vel_Upper
REAL(KIND=DP), INTENT(IN), DIMENSION(0:NY) :: DY
REAL(KIND=DP), INTENT(INOUT), DIMENSION(1:NX,1:NZ,0:NY+1) :: A, B, C, D
INTEGER :: I,J
DO J = 1, NZ
  DO I = 1, NX
IF (UW_vel_BC_Lower .EQ. 1) THEN
! Dirichlet
A(I,J,0)=0
A(I,J,1)=0
B(I,J,0)=1
B(I,J,1)=1
C(I,J,0)=0
C(I,J,1)=0
D(I,J,0)=0
D(I,J,1)=UW_vel_Lower
ELSE IF (UW_vel_BC_Lower .EQ. 2) THEN
! Neumann
A(I,J,0)=0
A(I,J,1)=0
B(I,J,0)=-1
B(I,J,1)=-1
C(I,J,0)=1
C(I,J,1)=1 !
D(I,J,0)=DY(0)*UW_vel_Lower ! (Y(1)-Y(0)=UW_vel_Lower*DYF(0))
D(I,J,1)=DY(1)*UW_vel_Lower ! (Y(2)-Y(1)=UW_vel_Lower*DYF(1))
END IF
IF (UW_vel_BC_Upper .EQ. 1) THEN
! Dirichlet
A(I,J,NY)=0
A(I,J,NY+1)=0
B(I,J,NY)=1
B(I,J,NY+1)=1
C(I,J,NY)=0
C(I,J,NY+1)=0
D(I,J,NY)=UW_vel_Upper
D(I,J,NY+1)=UW_vel_Upper
ELSE IF (UW_vel_BC_Upper .EQ. 2) THEN
! Neumann
A(I,J,NY)=-1
A(I,J,NY+1)=-1
B(I,J,NY)=1
B(I,J,NY+1)=1
C(I,J,NY)=0
C(I,J,NY+1)=0
D(I,J,NY)=DY(NY-1)*UW_vel_Upper
D(I,J,NY+1)=DY(NY)*UW_vel_Upper
END IF
END DO
END DO
END SUBROUTINE UW_Boundary_Conditions

SUBROUTINE Viscosity_Temperature(NX, NY, NZ, TH, THB, T_ref,Delta_T_Dim, DY, DYF,K_start, K_end, mu, mu_dbl_breve)
! Arrhenius- Type viscosity model
IMPLICIT NONE
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: NX, NY, NZ, K_start, K_end
REAL(KIND=DP), INTENT(IN) :: T_ref,Delta_T_Dim
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1), INTENT(IN) :: TH, THB
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1) :: T_total
REAL (KIND=DP), DIMENSION(0:NY), INTENT(IN)  :: DY, DYF
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1), INTENT(OUT) :: mu, mu_dbl_breve
REAL(KIND=DP) :: a=-2.10_DP, b=-4.45_DP, c=6.55_DP ! Model Paramaters
INTEGER :: I, J, K
! make temperature dimensional and then divide by reference temperature (input)
FORALL (I=1:NX, J=1:NZ, K=0:NY+1)
T_Total(I,J,K)=(TH(I,J,K)+THB(I,J,K))*Delta_T_Dim/T_ref
mu(I,J,K) = EXP(a + b/T_Total(I,J,K) + c/T_Total(I,J,K)**2)
END FORALL
FORALL  (I=1:NX, J=1:NZ, K=K_start:K_end)
! viscosity interpolated at the base grid (second order)
mu_dbl_breve(I,J,K)=(mu(I,J,K)*DYF(K-1)+mu(I,J,K-1)*DYF(K))/(2*DY(K-1))
END FORALL
END SUBROUTINE Viscosity_Temperature

SUBROUTINE courant( NX, NY, NZ, Lx, Ly, Lz, DYF, CFL, Ri, Pr, Re, U, V, W, &
                THB_wall_lower, delta_t )
! --------------------- This is used to obtain time step -----------------------
IMPLICIT NONE
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: NX, NY, NZ
REAL (KIND=DP), INTENT(OUT) :: delta_t
REAL (KIND=DP), INTENT(IN)  :: CFL, Ri, Pr, Re, Lx, Ly, Lz, THB_wall_lower
REAL (KIND=DP), DIMENSION(0:NY), INTENT(IN)  :: DYF
REAL(KIND=DP), PARAMETER :: pi=4.0_DP*ATAN(1.0_DP)
REAL (KIND=DP)  :: DX, DZ, Nmax, delta_t_x, delta_t_y, delta_t_z
REAL(KIND=DP), DIMENSION(1:NX,1:NZ,0:NY+1), INTENT(IN) :: U, V, W
INTEGER :: I, J, K
DX = Lx/NX
DZ = Lz/NZ
! Set the initial delta_t to some arbitrary large number
delta_t=999.d0
! Viscosity and diffusivity
delta_t=MIN( delta_t, 0.5_DP*MIN( DX,DZ )*Re )
delta_t=MIN( delta_t, 0.5_DP*MINVAL( DYF )*Re )
delta_t=MIN( delta_t, delta_t*Pr )
! Buoyancy period (for stratified flows)
IF (RI .NE. 0.0_DP) then
    Nmax    = sqrt(abs( Ri * THB_wall_lower ))
    delta_t = MIN(delta_t,0.1_DP*2.0_DP*pi/Nmax)
END IF
DO I=1,NX
DO J=1,NZ
DO K=1,NY
  delta_t_x=CFL*DX    /abs(U(I,J,K))
  delta_t_y=CFL*DYF(J)/abs(V(I,J,K))
  delta_t_z=CFL*DZ    /abs(W(I,J,K))
  delta_t=MIN( delta_t, delta_t_x, delta_t_y, delta_t_z )
END DO
END DO
END DO
IF (delta_t .le. 0.0_DP) then
  delta_t=0.0001_DP
END IF
END SUBROUTINE courant

SUBROUTINE RK_SOLVER ( K_start, K_end, NX, NY, NZ, TH_BC_Lower, TH_BC_Upper, &
U_BC_Lower,U_BC_Upper, V_BC_Lower,V_BC_Upper, W_BC_Lower, W_BC_Upper, &
kx, kz, gamma, zeta, alpha, delta_t, Lx, Lz, n1, n2, n3, plan_bkd, plan_fwd,&
 DY, DYF, Pr, Re, Ri, &
U, V, W,P, TH,  THB,mu, mu_dbl_breve, T_ref,Delta_T_Dim, U_wall_Lower,U_wall_Upper, &
 V_wall_Lower,V_wall_Upper, W_wall_Lower,W_wall_Upper,U_store,  V_store, W_store, TH_store )
! --------------- This is used to time step relevant equations -----------------
USE, INTRINSIC :: iso_c_binding
USE Fourier_Spectral
IMPLICIT NONE
include 'fftw3.f03'
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: K_start, K_end, NX, NY, NZ, TH_BC_Lower,TH_BC_Upper,&
 U_BC_Lower,U_BC_Upper, V_BC_Lower,V_BC_Upper, W_BC_Lower, W_BC_Upper
COMPLEX(C_DOUBLE_COMPLEX), PARAMETER :: ii=(0.0_DP, 1.0_DP)
COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(1:NX/2+1,1:NZ,0:NY+1) :: F_UTH, F_WTH, &
              F_TH, F_Exp_TH, F_Exp_U, F_Exp_V, F_Exp_W, F_V, F_Vx, F_Vz, &
              F_MU_Uy_p_Vx, F_MU_Wy_p_Vz, F_UV, F_VW, F_UU, F_UW, F_WW, F_U, F_W, &
              F_Ux, F_Uz, F_Wx, F_Wz, F_MU_Ux, F_MU_Wz, F_MU_Uz_p_Wx, F_P, F_Px, F_Pz
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1) :: Exp_TH_m1, Exp_TH, U_THpTHB, &
                Cranck_Exp_Th, THpTHB, THpTHB_Int, R_TH, A, B, C, D, W_THpTHB, &
                Exp_U, Exp_U_m1, Exp_V, Exp_V_m1, Exp_W, Exp_W_m1, U_breve, W_breve, &
                Uy, Wy, Vx, Vz, MU_Uy_p_Vx, MU_Wy_p_Vz, UV, VW, Cranck_Exp_U, &
                Cranck_Exp_V, Cranck_Exp_W, UU, UW, WW, Ux, Uz, Wx, Wz, MU_Ux, &
                MU_Wz, MU_Uz_p_Wx, Px, Pz
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1), INTENT(INOUT) :: U, V, W, TH, P,&
                          mu, mu_dbl_breve
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1,3), INTENT(OUT) ::U_store,  V_store, W_store, TH_store
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1), INTENT(IN)    :: THB
REAL(KIND=DP), DIMENSION (0:NY+1), INTENT(IN) :: DY, DYF
REAL(KIND=DP), DIMENSION (1:3), INTENT(IN) :: gamma, zeta, alpha
REAL(KIND=DP), INTENT(IN) :: delta_t, Lx, Lz, n1, n2, n3,Pr, Re, Ri, T_ref,Delta_T_Dim, &
  U_wall_Lower,U_wall_Upper, V_wall_Lower,V_wall_Upper, W_wall_Lower,W_wall_Upper
REAL(KIND=DP) :: alfa_t
type(C_PTR), INTENT(IN) :: plan_bkd, plan_fwd
REAL(KIND=DP), DIMENSION(1:Nx/2+1), INTENT(IN) :: kx
REAL(KIND=DP), DIMENSION(1:Nz), INTENT(IN) :: kz
INTEGER :: I, J, K, RK_step

Exp_TH_m1=0.0_DP
Exp_U_m1 = 0.0_DP
Exp_V_m1 = 0.0_DP
Exp_W_m1 = 0.0_DP

DO RK_step = 1,3
! 1. 1st solve the temperature equation
F_Exp_TH(1:NX/2+1,1:NZ,0:NY+1)=(0.0_DP, 0.0_DP)
THpTHB=TH+THB
THpTHB_Int=0.0_DP ! Must Initialise all variables
U_THpTHB=0.0_DP
W_THpTHB=0.0_DP
FORALL (I=1:NX, J=1:NZ, K=K_start:K_end)
  THpTHB_Int(I,J,K)=(THpTHB(I,J,K-1)*DYF(K)+THpTHB(I,J,K)*DYF(K-1))/(2*DY(K-1)) ! \breve\breve(T+To)
  U_THpTHB(I,J,K)=U(I,J,K)*THpTHB(I,J,K)
  W_THpTHB(I,J,K)=W(I,J,K)*THpTHB(I,J,K)
END FORALL
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, U_THpTHB, F_UTH )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, W_THpTHB, F_WTH )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, TH,       F_TH  )
FORALL (I=1:NX/2+1,J=1:NZ,K=0:NY+1)
F_Exp_TH(I,J,K) = -1*ii*kx(I)*F_UTH(I,J,K)-1*ii*kz(J)*F_WTH(I,J,K) + &
                  -1/(Pr*Re)*(kx(I)**2+kz(J)**2)*F_TH(I,J,K)
END FORALL
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Exp_TH, Exp_Th )
Exp_Th=0.0_DP
Cranck_Exp_Th=0.0_DP
FORALL (I=1:NX, J=1:NZ, K=K_start:K_end)
Exp_Th(I,J,K) = -1*(THpTHB_Int(I,J,K+1)*V(I,J,K+1) &
                -THpTHB_Int(I,J,K)*V(I,J,K))/DYF(K)+Exp_Th(I,J,K)
Cranck_Exp_Th(I,J,K) = 1/(Pr*Re*DYF(K))*((TH(I,J,K+1)-TH(I,J,K))/DY(K) &
                -(TH(I,J,K)-TH(I,J,K-1))/DY(K-1))
END FORALL
U_breve=0.0_DP
W_breve=0.0_DP
Uy=0.0_DP
Wy=0.0_DP
! 2. Solve the y- velocity
FORALL (I=1:NX, J=1:NZ, K=1:NY+1)
  ! U and W interpolated at the base grid (quasi second order)
  U_breve(I,J,K)= (U(I,J,K)*DYF(K)+U(I,J,K-1)*DYF(K-1))/(2*DY(K-1))
  W_breve(I,J,K)= (W(I,J,K)*DYF(K)+W(I,J,K-1)*DYF(K-1))/(2*DY(K-1))
  ! U and W wall normal gradients at the base grid
  Uy(I,J,K)=(U(I,J,K)-U(I,J,K-1))/DY(K-1)
  Wy(I,J,K)=(W(I,J,K)-W(I,J,K-1))/DY(K-1)
END FORALL

CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, V, F_V )
FORALL (I=1:NX/2+1,J=1:NZ,K=0:NY+1)
  F_Vx(I,J,K)=ii*kx(I)*F_V(I,J,K)
  F_Vz(I,J,K)=ii*kz(J)*F_V(I,J,K)
END FORALL
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Vx, Vx )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Vz, Vz )
MU_Uy_p_Vx=(Uy+Vx)*mu_dbl_breve
MU_Wy_p_Vz=(Wy+Vz)*mu_dbl_breve
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, MU_Uy_p_Vx, F_MU_Uy_p_Vx )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, MU_Wy_p_Vz, F_MU_Wy_p_Vz )
UV=U_breve*V
VW=V*W_breve
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, UV, F_UV )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, VW, F_VW )

FORALL (I=1:NX/2+1,J=1:NZ,K=0:NY+1)
  F_Exp_V(I,J,K)=ii*((kx(I)*F_MU_Uy_p_Vx(I,J,K)+kz(J)*F_MU_Wy_p_Vz(I,J,K))/Re &
    -(kx(I)*F_UV(I,J,K)+kz(J)*F_VW(I,J,K)))
END FORALL
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Exp_V, Exp_V )
Cranck_Exp_V=0.0_DP
FORALL (I=1:NX, J=1:NZ, K=K_start:K_end)
! Adding the temperature, y laplacian and the explicit pressure gradient
  Cranck_Exp_V(I,J,K) = Ri*n2*(TH(I,J,K-1)*DYF(K)+TH(I,J,K)*DYF(K-1))/(2.0_DP*DY(K-1)) + &
  2.0_DP/(Re*DY(K-1))*(mu(I,J,K)*(V(I,J,K+1)-V(I,J,K))/DYF(K)    &
  - mu(I,J,K-1)*(V(I,J,K)-V(I,J,K-1))/DYF(K-1)  )   &
  -2.0_DP*(P(I,J,K)-P(I,J,K-1))/DY(K-1)
END FORALL
! 3. Solve u and w - equation
UU = U*U
UW = U*W
WW = W*W
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, UU, F_UU )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, WW, F_WW )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, UW, F_UW )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, U, F_U )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, W, F_W )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, P, F_P )

FORALL (I=1:NX/2+1,J=1:NZ,K=0:NY+1)
  F_Ux(I,J,K)=ii*kx(I)*F_U(I,J,K)
  F_Uz(I,J,K)=ii*kz(J)*F_U(I,J,K)
  F_Wz(I,J,K)=ii*kz(J)*F_W(I,J,K)
  F_Wx(I,J,K)=ii*kx(I)*F_W(I,J,K)
  F_Px(I,J,K)=ii*kx(I)*F_P(I,J,K)
  F_Pz(I,J,K)=ii*kz(J)*F_P(I,J,K)
END FORALL

CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Ux, Ux )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Uz, Uz )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Px, Px )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Pz, Pz )

MU_Ux = 2*mu*Ux
MU_Wz = 2*mu*Wz
MU_Uz_p_Wx = mu*(Uz+Wx)
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, MU_Ux, F_MU_Ux )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, MU_Wz, F_MU_Wz )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, MU_Uz_p_Wx, F_MU_Uz_p_Wx )

! explicit terms for both u and w equations
FORALL (I=1:NX/2+1,J=1:NZ,K=0:NY+1)
F_Exp_U(I,J,K) = ii*(-kz(J)*F_UW(I,J,K) - kx(I)*F_UU(I,J,K) &
                + (kx(I)*F_MU_Ux(I,J,K) + kz(J)*F_MU_Uz_p_Wx(I,J,K))/Re)
F_Exp_W(I,J,K) = ii*(-kx(I)*F_UW(I,J,K) - kz(J)*F_WW(I,J,K) &
              + (kz(J)*F_MU_Wz(I,J,K) + kx(I)*F_MU_Uz_p_Wx(I,J,K))/Re)
END FORALL
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Exp_U, Exp_U )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Exp_W, Exp_W )
Cranck_Exp_U=0.0_DP
Cranck_Exp_W=0.0_DP
FORALL (I=1:NX, J=1:NZ, K=K_start:K_end)
Cranck_Exp_U(I,J,K) = (1/Re)*(mu_dbl_breve(I,J,K+1)*Vx(I,J,K+1)  - &
      mu_dbl_breve(I,J,K)*Vx(I,J,K))/DYF(K) - ((U(I,J,K+1)+U(I,J,K))*V(I,J,K+1) - &
      (U(I,J,K) + U(I,J,K-1))*V(I,J,K))/(2.0_DP*DYF(K)) + Ri*n1*TH(I,J,K) + &
      (mu_dbl_breve(I,J,K+1)*(U(I,J,K+1)-U(I,J,K))/DY(K) &
      - mu_dbl_breve(I,J,K+1)*(U(I,J,K)-U(I,J,K-1))/DY(K-1))/(Re*DYF(K)) &
      -2.0_DP*Px(I,J,K)
Cranck_Exp_W(I,J,K) = (1/Re)*(mu_dbl_breve(I,J,K+1)*Vz(I,J,K+1)  - &
      mu_dbl_breve(I,J,K)*Vz(I,J,K))/DYF(K) - ((W(I,J,K+1)+W(I,J,K))*V(I,J,K+1) - &
      (W(I,J,K) + W(I,J,K-1))*V(I,J,K))/(2.0_DP*DYF(K)) + Ri*n3*TH(I,J,K) + &
      (mu_dbl_breve(I,J,K+1)*(W(I,J,K+1)-W(I,J,K))/DY(K) &
      - mu_dbl_breve(I,J,K+1)*(W(I,J,K)-W(I,J,K-1))/DY(K-1))/(Re*DYF(K))&
      -2.0_DP*Pz(I,J,K)
END FORALL
!! ------------------- Finishing off the temperature equation ------------------
! IMPLICIT PART IS NOW ADDED
A=0.0_DP
B=0.0_DP
C=0.0_DP
DO J = 1, NZ
  DO I = 1, NX
    DO K = K_Start, K_End
      A(I,J,K) = -1/(Pr*Re)*1/(DYF(K)*DY(K-1))*(alpha(RK_step)/2.0_DP)*delta_t
      C(I,J,K) = -1/(Pr*Re)*1/(DYF(K)*DY(K))  *(alpha(RK_step)/2.0_DP)*delta_t
    END DO
  END DO
END DO
B = - A - C + 1.0_DP
TH = TH+delta_t*(gamma(RK_step)*Exp_Th + zeta(RK_step)*Exp_TH_m1 + &
                    (alpha(RK_step)/2.0_DP)*Cranck_Exp_Th)
Exp_TH_m1=Exp_TH
CALL Scalar_Boundary_Conditions(A,B,C,TH,NX,NY,NZ,TH_BC_Lower,TH_BC_Upper)
CALL Thomas_Matrix_Algorithm_real(A,B,C,TH,NX,NY,NZ)
!! -----------------------------------------------------------------------------
! Obtain the Viscosity from the temperature
CALL Viscosity_Temperature(NX, NY, NZ, TH, THB, T_ref,Delta_T_Dim, DY, DYF,K_start, K_end, mu, mu_dbl_breve)

! We already know TH at new time step, therefore we can keep it on the right
FORALL  (I=1:NX, J=1:NZ, K=K_start:K_end)
Cranck_Exp_V(I,J,K)=Cranck_Exp_V(I,J,K)+Ri*n2*(TH(I,J,K-1)*DYF(K)+TH(I,J,K)*DYF(K-1))/(2.0_DP*DY(K-1))
END FORALL
Cranck_Exp_U=Cranck_Exp_U+Ri*n1*TH
Cranck_Exp_W=Cranck_Exp_W+Ri*n3*TH
!! ------------------------ Finishing off the v equation -----------------------
! IMPLICIT PART IS NOW ADDED
A=0.0_DP
B=0.0_DP
C=0.0_DP
DO J = 1, NZ
  DO I = 1, NX
    DO K = 0, Ny+1
    A(I,J,K) = -2.0_DP*((V(I,J,K-1)/2+1/Re*mu(I,J,K-1)/DYF(K-1))/DY(K-1))*(alpha(RK_step)/2.0_DP)*delta_t
    C(I,J,K) = 2.0_DP*((V(I,J,K+1)/2-1/Re*mu(I,J,K)/DYF(K))/DY(K-1))*(alpha(RK_step)/2.0_DP)*delta_t
    B(I,J,K) = 1.0_DP+2.0_DP/Re*(mu(I,J,K)/DYF(K)+mu(I,J,K-1)/DYF(K-1))/DY(K-1)*(alpha(RK_step)/2.0_DP)*delta_t
      END DO
  END DO
END DO
V = V + delta_t * (gamma(RK_step)*Exp_V + zeta(RK_step)*Exp_V_m1 + (alpha(RK_step)/2.0_DP)*Cranck_Exp_V)
Exp_V_m1=Exp_V
CALL V_Boundary_Conditions(A,B,C,V,NX,NY,NZ,DYF,V_BC_Lower,V_BC_Upper,V_wall_Lower,V_wall_Upper)
CALL Thomas_Matrix_Algorithm_real(A,B,C,V,NX,NY,NZ)
!! -----------------------------------------------------------------------------
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, V, F_V )
FORALL (I=1:NX/2+1,J=1:NZ,K=0:NY+1)
  F_Vx(I,J,K)=ii*kx(I)*F_V(I,J,K)
  F_Vz(I,J,K)=ii*kz(J)*F_V(I,J,K)
END FORALL
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Vx, Vx )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Vz, Vz )

FORALL  (I=1:NX, J=1:NZ, K=K_start:K_end)
Cranck_Exp_U(I,J,K)=Cranck_Exp_U(I,J,K)+(1/Re)*(mu_dbl_breve(I,J,K+1)*Vx(I,J,K+1)&
      - mu_dbl_breve(I,J,K)*Vx(I,J,K))/DYF(K)
Cranck_Exp_W(I,J,K)=Cranck_Exp_W(I,J,K)+(1/Re)*(mu_dbl_breve(I,J,K+1)*Vz(I,J,K+1)&
      - mu_dbl_breve(I,J,K)*Vz(I,J,K))/DYF(K)
END FORALL
!! ------------------------ Finishing off the u equation -----------------------
! IMPLICIT PART IS NOW ADDED
A=0.0_DP
B=0.0_DP
C=0.0_DP
DO J = 1, NZ
  DO I = 1, NX
    DO K = K_Start, K_End
    A(I,J,K) = ( V(I,J,K)/2   - mu_dbl_breve(I,J,K)/(Re*DY(K-1)) ) &
    /(DYF(k))*(alpha(RK_step)/2.0_DP)*delta_t
    C(I,J,K) = (-V(I,J,K+1)/2 - mu_dbl_breve(I,J,K+1)/(Re*DY(K)) ) &
    /(DYF(k))*(alpha(RK_step)/2.0_DP)*delta_t
    B(I,J,K) = 1.0_DP+(((-V(I,J,K+1)+V(I,J,K))/2+ &
    mu_dbl_breve(I,J,K+1)/(Re*DY(K)) &
    + mu_dbl_breve(I,J,K)/(Re*DY(K-1)))/DYF(k))*(alpha(RK_step)/2.0_DP)*delta_t
      END DO
  END DO
END DO
U = U + delta_t * (gamma(RK_step)*Exp_U + zeta(RK_step)*Exp_U_m1 +(alpha(RK_step)/2.0_DP)*Cranck_Exp_U)
                   Exp_U_m1=Exp_U
CALL UW_Boundary_Conditions(A,B,C,U,NX,NY,NZ,DY,U_BC_Lower,U_BC_Upper,U_wall_Lower,U_wall_Upper)
CALL Thomas_Matrix_Algorithm_real(A,B,C,U,NX,NY,NZ)
!! ------------------------ Finishing off the w equation -----------------------
! IMPLICIT PART IS NOW ADDED
A=0.0_DP
B=0.0_DP
C=0.0_DP
DO J = 1, NZ
  DO I = 1, NX
    DO K = K_Start, K_End
    A(I,J,K) = ( V(I,J,K)/2   - mu_dbl_breve(I,J,K)/(Re*DY(K-1)) ) &
    /(DYF(k))*(alpha(RK_step)/2.0_DP)*delta_t
    C(I,J,K) = (-V(I,J,K+1)/2 - mu_dbl_breve(I,J,K+1)/(Re*DY(K)) ) &
    /(DYF(k))*(alpha(RK_step)/2.0_DP)*delta_t
    B(I,J,K) = 1.0_DP+(((-V(I,J,K+1)+V(I,J,K))/2+ mu_dbl_breve(I,J,K+1)/(Re*DY(K)) &
    + mu_dbl_breve(I,J,K)/(Re*DY(K-1)))/DYF(k))*(alpha(RK_step)/2.0_DP)*delta_t
      END DO
  END DO
END DO
W = W + delta_t * (gamma(RK_step)*Exp_W + zeta(RK_step)*Exp_W_m1 + (alpha(RK_step)/2.0_DP)*Cranck_Exp_W)
                    Exp_W_m1=Exp_W
CALL UW_Boundary_Conditions(A,B,C,W,NX,NY,NZ,DY,W_BC_Lower,W_BC_Upper,W_wall_Lower,W_wall_Upper)
CALL Thomas_Matrix_Algorithm_real(A,B,C,W,NX,NY,NZ)
!! -----------------------------------------------------------------------------
alfa_t=alpha(RK_step)*delta_t
CALL Remove_Divergence(NX, NY, NZ, Lx, Lz, alfa_t, kx, kz, DY, DYF, plan_fwd, &
                       plan_bkd, U, V, W, P)

U_store(:,:,:,RK_step)=U
V_store(:,:,:,RK_step)=V
W_store(:,:,:,RK_step)=W
TH_store(:,:,:,RK_step)=TH

END DO
END SUBROUTINE RK_SOLVER

SUBROUTINE Dissipation_Calculation( plan_bkd, plan_fwd, kx, kz, DY, DYF, u, v, w, TH, mutot,Dissipation)
USE Convergence_Check
USE Integrate_Volume
USE Fourier_Spectral
IMPLICIT NONE
include 'fftw3.f03'


type(C_PTR), INTENT(IN) :: plan_bkd, plan_fwd
REAL(KIND=DP), DIMENSION(1:Nx/2+1), INTENT(IN) :: kx
REAL(KIND=DP), DIMENSION(1:Nz), INTENT(IN) :: kz
REAL(KIND=DP), DIMENSION (0:NY), INTENT(IN) :: DY, DYF
COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(1:NX/2+1,1:NZ,0:NY+1) :: F_U, F_V, F_W, F_TH, &
F_Ux, F_Uz, F_Vx, F_Vz, F_Wx, F_Wz, F_THx, F_THz
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1), INTENT(IN) :: u, v, w, TH, mutot
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1):: ux, uy, uz, vx, vy, vz, THx, &
THy, THz, U_dbl_breve, W_dbl_breve, TH_dbl_breve, mutot_dbl_breve, TH_terms

REAL(KIND=DP) :: Diss_1, Diss_2, Diss_3, TH_Dis
REAL(KIND=DP), INTENT(OUT) :: Dissipation
COMPLEX(C_DOUBLE_COMPLEX), PARAMETER :: ii=(0.0_DP, 1.0_DP)

CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, u, F_U )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, v, F_V )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, w, F_W )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, TH, F_TH )
FORALL (I=1:NX/2+1,J=1:NZ,K=0:NY+1)
F_Ux(I,J,K)=ii*kx(I)*F_U(I,J,K)
F_Uz(I,J,K)=ii*kz(J)*F_U(I,J,K)
F_Vx(I,J,K)=ii*kx(I)*F_V(I,J,K)
F_Vz(I,J,K)=ii*kz(J)*F_V(I,J,K)
F_Wx(I,J,K)=ii*kx(I)*F_W(I,J,K)
F_Wz(I,J,K)=ii*kz(J)*F_W(I,J,K)

F_THx(I,J,K)=ii*kx(I)*F_TH(I,J,K)
F_THz(I,J,K)=ii*kz(J)*F_TH(I,J,K)
END FORALL
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Ux, Ux )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Uz, Uz )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Vx, Vx )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Vz, Vz )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Wx, Wx )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Wz, Wz )

CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_THx, THx )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_THz, THz )

U_dbl_breve=0.0_DP
W_dbl_breve=0.0_DP
Uy=0.0_DP
Wy=0.0_DP

FORALL (I=1:NX, J=1:NZ, K=1:NY+1)
  ! U and W interpolated at the base grid (second order)
  U_dbl_breve(I,J,K)=(U(I,J,K)*DYF(K-1)+U(I,J,K-1)*DYF(K))/(2*DY(K-1))
  W_dbl_breve(I,J,K)=(W(I,J,K)*DYF(K-1)+W(I,J,K-1)*DYF(K))/(2*DY(K-1))
  TH_dbl_breve(I,J,K)=(TH(I,J,K)*DYF(K-1)+TH(I,J,K-1)*DYF(K))/(2*DY(K-1))
  mutot_dbl_breve(I,J,K)=(mutot(I,J,K)*DYF(K-1)+mutot(I,J,K-1)*DYF(K))/(2*DY(K-1))
  ! U and W wall normal gradients at the base grid
  Uy(I,J,K)=(U_dbl_breve(I,J,K)-U_dbl_breve(I,J,K-1))/DYF(K-1)
  Wy(I,J,K)=(W_dbl_breve(I,J,K)-W_dbl_breve(I,J,K-1))/DYF(K-1)
  THy(I,J,K)=(TH_dbl_breve(I,J,K)-TH_dbl_breve(I,J,K-1))/DYF(K-1)
END FORALL

FORALL (I=1:NX, J=1:NZ, K=1:NY+1)
Vy(I,J,K)=(V(I,J,K+1)-V(I,J,K-1))/(2.0_DP*DYF(K-1))
END FORALL
FORALL (I=1:NX, J=1:NZ)
Vy(I,J,NY+1)=(V(I,J,NY+1)-V(I,J,NY))/DYF(NY)
END FORALL

FORALL (I=1:NX, J=1:NZ, K=1:NY+1)
Ux(I,J,K)=Ux(I,J,K)*sqrt(mutot(I,J,K))
Uy(I,J,K)=Uy(I,J,K)*sqrt(mutot(I,J,K))
Uz(I,J,K)=Uz(I,J,K)*sqrt(mutot(I,J,K))
Wx(I,J,K)=Wx(I,J,K)*sqrt(mutot(I,J,K))
Wy(I,J,K)=Wy(I,J,K)*sqrt(mutot(I,J,K))
Wz(I,J,K)=Wz(I,J,K)*sqrt(mutot(I,J,K))
Vx(I,J,K)=Vx(I,J,K)*sqrt(mutot_dbl_breve(I,J,K))
Vy(I,J,K)=Vy(I,J,K)*sqrt(mutot_dbl_breve(I,J,K))
Vz(I,J,K)=Vz(I,J,K)*sqrt(mutot_dbl_breve(I,J,K))
END FORALL
CALL Vector_Volume_Integral(Ux, Vx, Wx,Ux, Vx, Wx,Nx,Ny,Nz,DY, DYF,Lx,Ly,Lz,Diss_1)
CALL Vector_Volume_Integral(Uy, Vy, Wy,Uy, Vy, Wy,Nx,Ny,Nz,DY, DYF,Lx,Ly,Lz,Diss_2)
CALL Vector_Volume_Integral(Uz, Vz, Wz,Uz, Vz, Wz,Nx,Ny,Nz,DY, DYF,Lx,Ly,Lz,Diss_3)
TH_terms=(THx*THx+THy*THy+THz*THz)*Ri/(Pr*T_ref**2)
CALL Integrate_Volume(TH_terms,Nx,Ny,Ny,Nz,DYF_mod,Lx,Ly,Lz,TH_Dis)
Dissipation=Diss_1+Diss_2+Diss_3+TH_Dis
END SUBROUTINE Dissipation_Calculation

END MODULE Channel_Solvers_BC
