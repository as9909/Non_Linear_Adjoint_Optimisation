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
SUBROUTINE Remove_Divergence(NX, NY, NZ, Lx, Lz, kx, kz, DY, DYF, plan_fwd, &
                              plan_bkd, U, V, W)
! -------------- This is used to remove divergence from velocity ---------------
USE, INTRINSIC :: iso_c_binding
USE Fourier_Spectral
IMPLICIT NONE
include 'fftw3.f03'
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
COMPLEX(KIND=DP),PARAMETER :: ii=(0.d0, 1.d0)
INTEGER, INTENT(IN) :: NX, NY, NZ
INTEGER :: I, J, K
COMPLEX(KIND=DP), DIMENSION(1:NX/2+1,1:NZ,0:NY+1) :: F_U, F_V, F_W, Q
REAL(KIND=DP), DIMENSION(1:NX/2,1:NZ,0:NY+1), INTENT(INOUT) :: U, V, W
REAL(KIND=DP), DIMENSION(1:NX/2+1,1:NZ,0:NY+1) :: A, B, C
!COMPLEX(KIND=DP), DIMENSION(1:NX/2+1,0:NY+1):: Vel_Div
REAL(KIND=DP), INTENT(IN) :: Lx, Lz
REAL(KIND=DP), DIMENSION(1:NX/2+1), INTENT(IN) :: kx
REAL(KIND=DP), DIMENSION(1:NZ), INTENT(IN) :: kz
REAL(KIND=DP), DIMENSION(0:NY), INTENT(IN) :: DY, DYF
type(C_PTR), INTENT(IN) :: plan_fwd, plan_bkd
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, U, F_U)
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, V, F_V)
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, W, F_W)
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
      Q(I,J,K) = ii*kx(I)*F_U(I,J,K)+(F_V(I,J,K+1)-F_V(I,J,K))/DYF(K) &
                     +ii*kz(J)*F_W(I,J,K)
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
  !FORALL (I = 1: NX/2+1, K = 1: NY) ! Ignore the ghost cells
  !    Q(I,J,K)=Vel_Div(I,K)
  !END FORALL
! Remove Divergence from the original velocity field
DO J = 1, NZ
  DO I = 1, NX/2+1
    DO K = 1, NY
      F_U(I,J,K) = F_U(I,J,K) - ii * kx(I) * Q(I,J,K)
      IF (K .ge. 2) THEN ! Ignore the ghost cells
        F_V(I,J,K) = F_V(I,J,K) - (Q(I,J,K)-Q(I,J,K-1))/DY(K-1)
      END IF
      F_W(I,J,K) = F_W(I,J,K) - ii*kz(J)*Q(I,J,K)
    END DO
  END DO
END DO
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_U, U)
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_V, V)
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_W, W)
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
  C(I,J,1) = C(I,J,1)/B(I,J,1)
  D(I,J,1) = D(I,J,1)/B(I,J,1)
  DO K = 2,NY ! Ignore the ghost cells
    C(I,J,K) =  C(I,J,K)/(B(I,J,K)-A(I,J,K)*C(I,J,K-1))
    D(I,J,K) =  (D(I,J,K)-A(I,J,K)*D(I,J,K-1))/(B(I,J,K)-A(I,J,K)*C(I,J,K-1))
  END DO
  DO K = NY-1,1,-1 ! Ignore the ghost cells
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
  C(I,J,1) = C(I,J,1)/B(I,J,1)
  D(I,J,1) = D(I,J,1)/B(I,J,1)
  DO K = 2,NY ! Ignore the ghost cells
    C(I,J,K) =  C(I,J,K)/(B(I,J,K)-A(I,J,K)*C(I,J,K-1))
    D(I,J,K) =  (D(I,J,K)-A(I,J,K)*D(I,J,K-1))/(B(I,J,K)-A(I,J,K)*C(I,J,K-1))
  END DO
  DO K = NY-1,1,-1 ! Ignore the ghost cells
    D(I,J,K) = D(I,J,K)-C(I,J,K)*D(I,J,K+1)
  END DO
END DO
END DO
END SUBROUTINE Thomas_Matrix_Algorithm_real

! --------------- 4. Subroutine for velocity boundary conditions ---------------
SUBROUTINE Velocity_Boundary_Conditions(U_BC_Lower, U_BC_Upper, V_BC_Lower, &
  V_BC_Upper, W_BC_Lower, W_BC_Upper,U_wall_lower, V_wall_upper, W_wall_lower, &
   U_wall_upper, V_wall_lower, W_wall_upper, NX, NY, NZ, DY, DYF, U, V, W)
IMPLICIT NONE
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: U_BC_Lower, U_BC_Upper, V_BC_Lower, V_BC_Upper, &
                       W_BC_Lower, W_BC_Upper
REAL(KIND=DP), INTENT(IN) :: U_wall_lower, V_wall_upper, W_wall_lower, &
                             U_wall_upper, V_wall_lower, W_wall_upper
REAL(KIND=DP), DIMENSION(1:NX,1:NZ,0:NY+1), INTENT(INOUT) :: U, V, W
INTEGER, INTENT(IN) :: NX, NY, NZ
INTEGER :: I, J
REAL(KIND=DP), DIMENSION(0:NY), INTENT(IN) :: DY, DYF
! --------------------- Set Lower Wall Boundary Conditions ---------------------
IF (U_BC_Lower .EQ. 1) THEN ! Dirichlet Conditions
  DO I = 1,NX
    DO J = 1, NZ
      U(I,J,0)  = U_wall_lower
      U(I,J,1)  = U_wall_lower
      END DO
    END DO
ELSE IF (U_BC_Lower .EQ. 2) THEN ! Neumann Conditions
  DO I = 1,NX
    DO J = 1, NZ
      U(I,J,1)  = U(I,J,2)-DY(1)*U_wall_lower
      U(I,J,0)  = U(I,J,1)-DY(0)*U_wall_lower
    END DO
  END DO
ELSE
  STOP 'Boundary Conditions entered are not approriate'
END  IF
IF (W_BC_Lower .EQ. 1) THEN ! Dirichlet Conditions
  DO I = 1,NX
    DO J = 1, NZ
      W(I,J,0)  = W_wall_lower
      W(I,J,1)  = W_wall_lower
      END DO
    END DO
ELSE IF (U_BC_Lower .EQ. 2) THEN ! Neumann Conditions
  DO I = 1,NX
    DO J = 1, NZ
      W(I,J,1)  = W(I,J,2)-DY(1)*W_wall_lower
      W(I,J,0)  = W(I,J,1)-DY(0)*W_wall_lower
    END DO
  END DO
ELSE
  STOP 'Boundary Conditions entered are not approriate'
END  IF
IF (V_BC_Lower .EQ. 1) THEN ! Dirichlet Conditions
  DO I = 1,NX
    DO J = 1, NZ
      V(I,J,1)  = 2.0_DP*V_wall_lower-V(I,J,2)
      V(I,J,0)  = V(I,J,1)
    END DO
  END DO
ELSE IF (V_BC_Lower .EQ. 2) THEN ! Neumann Conditions
  DO I = 1,NX
    DO J = 1, NZ
      V(I,J,1)  = V(I,J,2)-DYF(1)*V_wall_lower
      V(I,J,1)  = V(I,J,1)-DYF(0)*V_wall_lower
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
ELSE IF (U_BC_Lower .EQ. 2) THEN ! Neumann Conditions
  DO I = 1,NX
    DO J = 1, NZ
      U(I,J,NY)   = U(I,J,NY-1)+DY(NY-1)*U_wall_lower
      U(I,J,NY+1) = U(I,J,NY)+DY(NY)*U_wall_lower
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
ELSE IF (W_BC_Lower .EQ. 2) THEN ! Neumann Conditions
  DO I = 1,NX
    DO J = 1, NZ
      W(I,J,NY)   = W(I,J,NY-1)+DY(NY-1)*W_wall_lower
      W(I,J,NY+1) = W(I,J,NY)+DY(NY)*W_wall_lower
    END DO
  END DO
ELSE
  STOP 'Boundary Conditions entered are not approriate'
END  IF
IF (V_BC_Upper .EQ. 1) THEN ! Dirichlet Conditions
  DO I = 1,NX
    DO J = 1, NZ
      V(I,J,NY+1) = 2.0_DP*V_wall_upper-V(I,J,NY)
    END DO
  END DO
ELSE IF (V_BC_Upper .EQ. 2) THEN ! Neumann Conditions
  DO I = 1,NX
    DO J = 1, NZ
      V(I,J,NY+1) = V(I,J,NY)+DYF(NY)*V_wall_upper
    END DO
  END DO
ELSE
  STOP 'Boundary Conditions entered are not approriate'
END IF
END SUBROUTINE Velocity_Boundary_Conditions

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
kx, kz, gamma, zeta, alpha, delta_t, plan_bkd, plan_fwd, DY, DYF, Pr, Re, Ri, &
U, V, W, TH, THB )
! --------------- This is used to time step relevant equations -----------------
USE, INTRINSIC :: iso_c_binding
USE Fourier_Spectral
IMPLICIT NONE
include 'fftw3.f03'
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: K_start, K_end, NX, NY, NZ, TH_BC_Lower,TH_BC_Upper
COMPLEX(C_DOUBLE_COMPLEX), PARAMETER :: ii=(0.0_DP, 1.0_DP)
COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(1:NX/2+1,1:NZ,0:NY+1) :: F_UTH, F_WTH, &
              F_TH, F_Exp_TH
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1) :: Exp_TH_m1, Exp_TH, U_THpTHB, &
                Cranck_Exp_Th, THpTHB, THpTHB_Int, R_TH, A, B, C, D, W_THpTHB
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1), INTENT(INOUT) :: U, V, W, TH
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1), INTENT(IN)    :: THB
REAL(KIND=DP), DIMENSION (0:NY+1), INTENT(IN) :: DY, DYF
REAL(KIND=DP), DIMENSION (1:3), INTENT(IN) :: gamma, zeta, alpha
REAL(KIND=DP), INTENT(IN) :: delta_t
REAL(KIND=DP), INTENT(IN) :: Pr, Re, Ri
type(C_PTR), INTENT(IN) :: plan_bkd, plan_fwd
REAL(KIND=DP), DIMENSION(1:Nx/2+1), INTENT(IN) :: kx
REAL(KIND=DP), DIMENSION(1:Nz), INTENT(IN) :: kz
INTEGER :: I, J, K, RK_step
Exp_TH_m1=0.0_DP
DO RK_step = 1,3
F_Exp_TH(1:NX/2+1,1:NZ,0:NY+1)=(0.0_DP, 0.0_DP)
THpTHB=TH+THB
FORALL (I=1:NX, J=1:NZ, K=K_start:K_end)
THpTHB_Int(I,J,K)=(THpTHB(I,J,K+1)*DYF(K)+THpTHB(I,J,K)*DYF(K+1))/(2*DY(K))
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

FORALL (I=1:NX, J=1:NZ, K=K_start:K_end)
  Exp_Th(I,J,K) = -1*(THpTHB_Int(I,J,K+1)*V(I,J,K+1) &
                  -THpTHB_Int(I,J,K)*V(I,J,K))/DYF(K)+Exp_Th(I,J,K)
  Cranck_Exp_Th(I,J,K) = 1/(Pr*Re*DYF(K))*((TH(I,J,K+1)-TH(I,J,K))/DY(K)-(TH(I,J,K)-TH(I,J,K-1))/DY(K-1))
END FORALL
TH = TH+delta_t*(gamma(RK_step)*Exp_Th + zeta(RK_step)*Exp_TH_m1 + &
                    (alpha(RK_step)/2.0_DP)*Cranck_Exp_Th)
Exp_TH_m1=Exp_TH
DO J = 1, NZ
  DO I = 1, NX
    DO K = K_Start, K_End
      A(I,J,K) = -1/(Pr*Re)*1/(DYF(K)*DY(K-1))*(alpha(RK_step)/2.0_DP)*delta_t
      C(I,J,K) = -1/(Pr*Re)*1/(DYF(K)*DY(K))  *(alpha(RK_step)/2.0_DP)*delta_t
    END DO
  END DO
END DO
B = - A - C + 1.0_DP
CALL Scalar_Boundary_Conditions(A,B,C,TH,NX,NY,NZ,TH_BC_Lower,TH_BC_Upper)
CALL Thomas_Matrix_Algorithm_real(A,B,C,TH,NX,NY,NZ)
END DO
END SUBROUTINE RK_SOLVER

END MODULE Channel_Solvers_BC
