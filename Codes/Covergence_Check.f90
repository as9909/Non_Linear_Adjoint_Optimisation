MODULE Covergence_Check

SUBROUTINE Base_to_Fraction(A, Nx, Ny, Nz, DYF, DYF_mod)
IMPLICIT NONE
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: Nx, Ny, Nz
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1), INTENT(INOUT):: A
REAL(KIND=DP), DIMENSION (0:NY), INTENT(IN):: DYF
REAL(KIND=DP), DIMENSION (0:NY), INTENT(OUT):: DYF_mod
REAL(KIND=DP), INTENT(IN):: delta_t

INTEGER :: I, J, K

FORALL  (I=1:NX, J=1:NZ)
A(I,J,1)=(A(I,J,1)+A(I,J,2))/2
A(I,J,NY+1)=(A(I,J,NY+1)+A(I,J,NY))/2
END FORALL
DYF_mod(1)=DYF(1)/2
DYF_mod(NY)=DYF(NY)/2
END SUBROUTINE Base_to_Fraction

SUBROUTINE Integrate_Volume(A,Nx,Ny,N_y,Nz,DY_F,Lx,Ly,Lz,A_int)
IMPLICIT NONE
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: Nx, Ny, Nz, N_y
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1), INTENT(IN):: A
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ):: A_bar
REAL(KIND=DP), DIMENSION (1:NX):: A_bar_bar
REAL(KIND=DP), DIMENSION (1), INTENT(IN):: Lx, Ly, Lz
REAL(KIND=DP), DIMENSION (0:NY), INTENT(IN) :: DY_F
A_bar=0.0_DP
A_bar_bar=0.0_DP
A_int=0.0_DP
FORALL  (I=1:NX, J=1:NZ, K=2:N_y)
A_bar(I,J)=A_bar(I,J)+DY_F(K-1)*(A(I,J,K-1)+A(I,J,K))/2.0_DP
END FORALL
FORALL  (I=1:NX)
A_bar_bar(I)=A_bar_bar(I)+(Lz/Nz)*(A_bar(I,J-1)+A_bar(I,J))/2.0_DP
END FORALL
FORALL  (I=1:NX)
A_int=A_int+(Lx/Nx)*(A_bar_bar(I-1)+A_bar_bar(I))/2.0_DP
END FORALL
A_int=A_int/(Lx*Ly*Lz)
END SUBROUTINE Integrate_Volume

SUBROUTINE Vector_Volume_Integral(a1, a2, a3, b1, b2, b3 ,Nx,Ny,Nz,DY, DYF,Lx,Ly,Lz,Vec_norm)

INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: Nx, Ny, Nz
INTEGER :: Nyp1
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1), INTENT(IN):: a1, a2, a3, b1, b2, b3
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1):: a1b1, a2b2, a3b3
REAL(KIND=DP), DIMENSION (0:NY), INTENT(IN) :: DY, DYF
REAL(KIND=DP), DIMENSION (1) :: Ax_int, Ay_int, Az_int
REAL(KIND=DP), DIMENSION (1), INTENT(OUT) :: Vec_norm
Nyp1=Ny+1
a1b1=a1*b1
a2b2=a2*b2
a3b3=a3*b3
CALL Integrate_Volume(a1b1,Nx,Ny,Ny,Nz,DY,Lx,Ly,Lz,Ax_int)
CALL Integrate_Volume(a3b3,Nx,Ny,Ny,Nz,DY,Lx,Ly,Lz,Az_int)

CALL Base_to_Fraction(a2b2, Nx, Ny, Nz, DYF, DYF_mod)
CALL Integrate_Volume(a2b2,Nx,Ny,Nyp1,Nz,DYF_mod,Lx,Ly,Lz,Ay_int)
Vec_norm=Ax_int+Az_int+Ay_int
END SUBROUTINE Vector_Volume_Integral

SUBROUTINE update_initial_real_fields(Nx, Ny, Nz, Lx, Ly, Lz, Ri, T_ref, eps, &
                              DY, DYF,  v1o, v2o, v3o, stauo, Uo, Vo, Wo, THo)
IMPLICIT NONE
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: Nx, Ny, Nz
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1), INTENT(IN):: v1o, v2o, v3o, stauo
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1), INTENT(INOUT):: Uo, Vo, Wo, THo
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1):: THo_sq, stauo_THo&
check_c1_1, check_c1_2, check_c1_3, check_c1_1, check_c1_2, check_c1_3, &
check_c2_1, check_c2_2, check_c2_3, check_c2_1, check_c2_2, check_c2_3
REAL(KIND=DP), DIMENSION (1), INTENT(IN):: Lx, Ly, Lz, Ri, T_ref
REAL(KIND=DP), DIMENSION (1), INTENT(INOUT):: eps
REAL(KIND=DP), DIMENSION (0:NY), INTENT(IN) :: DY, DYF
REAL(KIND=DP), DIMENSION(1) :: c1, c2, det, c1_norm, c2_norm, c, vel_norm_sq,  &
TH_norm_sq, Eo, stauoTHo_int, vouo_int, adj_vel_norm_sq, adj_TH_norm_sq

CALL Vector_Volume_Integral(Uo, Vo, Wo, Uo, Vo, Wo,Nx,Ny,Nz,DY, DYF,Lx,Ly,Lz,vel_norm_sq)
THo_sq=THo*THo
Integrate_Volume(THo_sq,Nx,Ny,Ny,Nz,DY_F,Lx,Ly,Lz,TH_norm_sq)
Eo=(1.0_DP/2.0_DP)*(vel_norm_sq+(Ri/(T_ref**2))*TH_norm_sq)

CALL Vector_Volume_Integral(v1o,v2o,v3o,v1o,v2o,v3o,Nx,Ny,Nz,DY, DYF,Lx,Ly,Lz,adj_vel_norm_sq)
stauo_sq=stauo*stauo
Integrate_Volume(stauo_sq,Nx,Ny,Ny,Nz,DY_F,Lx,Ly,Lz,adj_TH_norm_sq)
Eo_adj=(1.0_DP/2.0_DP)*(adj_vel_norm_sq+(Ri/(T_ref**2))*adj_TH_norm_sq)

CALL Vector_Volume_Integral(Uo, Vo, Wo, v1o, v2o, v3o,Nx,Ny,Nz,DY, DYF,Lx,Ly,Lz,vouo_int)
stauo_THo=stauo*THo
Integrate_Volume(stauo_THo,Nx,Ny,Ny,Nz,DY_F,Lx,Ly,Lz,stauoTHo_int)

det=-1.0_DP
DO WHILE (det .lt. 0.0_DP)
alpha=2.0_DP*eps*Eo
beta=-2.0_DP*(2.0_DP*Eo+eps*(vouo_int+stauoTHo_int))
gamma=2.0_DP*(vouo_int+stauoTHo_int)+eps*Eo_adj
det=beta**2-4.0_DP*alpha*gamma
eps=eps/4.0_DP
END DO
c1=(-1.0*beta-SQRT(det))/2.0_DP*alpha)
c2=(-1.0*beta+SQRT(det))/2.0_DP*alpha)

check_c1_1=eps*(v1o-2*c1*Uo)
check_c1_2=eps*(v2o-2*c1*Vo)
check_c1_3=eps*(v3o-2*c1*Wo)
Vector_Volume_Integral(check_c1_1, check_c1_2, check_c1_3, check_c1_1, check_c1_2, &
check_c1_3 ,Nx,Ny,Nz,DY, DYF,Lx,Ly,Lz,c1_norm)
check_c2_1=eps*(v1o-2*c2*Uo)
check_c2_2=eps*(v2o-2*c2*Vo)
check_c2_3=eps*(v3o-2*c2*Wo)
Vector_Volume_Integral(check_c2_1, check_c2_2, check_c2_3, check_c2_1, check_c2_2, &
check_c2_3 ,Nx,Ny,Nz,DY, DYF,Lx,Ly,Lz, c2_norm)
c=c1_norm
IF (c2_norm .lt. c1_norm)
c=c2_norm
END IF

Uo=Uo+eps*(v1o-c*Uo)
Vo=Vo+eps*(v2o-c*Vo)
Wo=Wo+eps*(v3o-c*Wo)
THo=THo+((eps*T_ref**2)/Ri)*(stauo-c*THo*Ri/(T_ref**2))

END SUBROUTINE update_initial_real_fields
END MODULE Covergence_Check
