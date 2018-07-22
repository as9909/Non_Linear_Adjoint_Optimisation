MODULE Convergence_Check
CONTAINS

SUBROUTINE Base_to_Fraction(A, Nx, Ny, Nz, DYF, DYF_mod)
IMPLICIT NONE
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: Nx, Ny, Nz
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1), INTENT(INOUT):: A
REAL(KIND=DP), DIMENSION (0:NY), INTENT(IN):: DYF
REAL(KIND=DP), DIMENSION (0:NY), INTENT(OUT):: DYF_mod
INTEGER :: I, J, K
DYF_mod=0.0_DP
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
INTEGER :: I, J, K
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1), INTENT(IN):: A
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ):: A_bar
REAL(KIND=DP), DIMENSION (1:NX):: A_bar_bar
REAL(KIND=DP)  , INTENT(IN):: Lx, Ly, Lz
REAL(KIND=DP), DIMENSION (0:NY), INTENT(IN) :: DY_F
REAL(KIND=DP), INTENT(OUT) :: A_int
A_bar=0.0_DP
A_bar_bar=0.0_DP
A_int=0.0_DP
DO I=1,NX
DO J=1,NZ
DO K=2,N_y
A_bar(I,J)=A_bar(I,J)+DY_F(K-1)*(A(I,J,K-1)+A(I,J,K))/2.0_DP
END DO
END DO
END DO
DO I=1,NX
DO J=1,NZ
A_bar_bar(I)=A_bar_bar(I)+(Lz/Nz)*(A_bar(I,J-1)+A_bar(I,J))/2.0_DP
END DO
END DO
DO I=1,NX
A_int=A_int+(Lx/Nx)*(A_bar_bar(I-1)+A_bar_bar(I))/2.0_DP
END DO
A_int=A_int/(Lx*Ly*Lz)
END SUBROUTINE Integrate_Volume

SUBROUTINE Vector_Volume_Integral(a1, a2, a3, b1, b2, b3 ,Nx,Ny,Nz,DY, DYF,Lx,Ly,Lz,Vec_norm)

INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: Nx, Ny, Nz
INTEGER :: Nyp1
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1), INTENT(IN):: a1, a2, a3, b1, b2, b3
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1):: a1b1, a2b2, a3b3
REAL(KIND=DP), DIMENSION (0:NY), INTENT(IN) :: DY, DYF
REAL(KIND=DP), DIMENSION (0:NY) :: DYF_mod
REAL(KIND=DP)  , INTENT(IN):: Lx, Ly, Lz
REAL(KIND=DP)   :: Ax_int, Ay_int, Az_int
REAL(KIND=DP)  , INTENT(OUT) :: Vec_norm
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
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1):: THo_sq, stauo_THo, stauo_sq, &
check_c1_1, check_c1_2, check_c1_3, &
check_c2_1, check_c2_2, check_c2_3
REAL(KIND=DP)  , INTENT(IN):: Lx, Ly, Lz, Ri, T_ref
REAL(KIND=DP)  , INTENT(INOUT):: eps
REAL(KIND=DP), DIMENSION (0:NY), INTENT(IN) :: DY, DYF
REAL(KIND=DP) :: kc1, kc2, det, c1_norm, c2_norm, kc, vel_norm_sq,  &
TH_norm_sq, Eo, stauoTHo_int, vouo_int, adj_vel_norm_sq, adj_TH_norm_sq, &
alpha, beta, gamma, Eo_adj

CALL Vector_Volume_Integral(Uo, Vo, Wo, Uo, Vo, Wo,Nx,Ny,Nz,DY,DYF,Lx,Ly,Lz,vel_norm_sq)
THo_sq=THo*THo
CALL Integrate_Volume(THo_sq,Nx,Ny,Ny,Nz,DY,Lx,Ly,Lz,TH_norm_sq)
Eo=(1.0_DP/2.0_DP)*(vel_norm_sq+(Ri/(T_ref**2))*TH_norm_sq)

CALL Vector_Volume_Integral(v1o,v2o,v3o,v1o,v2o,v3o,Nx,Ny,Nz,DY, DYF,Lx,Ly,Lz,adj_vel_norm_sq)
stauo_sq=stauo*stauo
CALL Integrate_Volume(stauo_sq,Nx,Ny,Ny,Nz,DY,Lx,Ly,Lz,adj_TH_norm_sq)
Eo_adj=(1.0_DP/2.0_DP)*(adj_vel_norm_sq+(Ri/(T_ref**2))*adj_TH_norm_sq)

CALL Vector_Volume_Integral(Uo, Vo, Wo, v1o, v2o, v3o,Nx,Ny,Nz,DY, DYF,Lx,Ly,Lz,vouo_int)
stauo_THo=stauo*THo
CALL Integrate_Volume(stauo_THo,Nx,Ny,Ny,Nz,DY,Lx,Ly,Lz,stauoTHo_int)

det=-1.0_DP
DO WHILE (det .lt. 0.0_DP)
alpha=2.0_DP*eps*Eo
beta=-2.0_DP*(2.0_DP*Eo+eps*(vouo_int+stauoTHo_int))
gamma=2.0_DP*(vouo_int+stauoTHo_int)+eps*Eo_adj
det=beta**2-4.0_DP*alpha*gamma
eps=eps/4.0_DP
END DO
 kc1=(-1.0*beta-SQRT(det))/(2.0_DP*alpha)
 kc2=(-1.0*beta+SQRT(det))/(2.0_DP*alpha)

check_c1_1=eps*(v1o-2*kc1*Uo)
check_c1_2=eps*(v2o-2*kc1*Vo)
check_c1_3=eps*(v3o-2*kc1*Wo)
CALL Vector_Volume_Integral(check_c1_1, check_c1_2, check_c1_3, check_c1_1, check_c1_2, &
check_c1_3 ,Nx,Ny,Nz,DY, DYF,Lx,Ly,Lz,c1_norm)
check_c2_1=eps*(v1o-2*kc2*Uo)
check_c2_2=eps*(v2o-2*kc2*Vo)
check_c2_3=eps*(v3o-2*kc2*Wo)
CALL Vector_Volume_Integral(check_c2_1, check_c2_2, check_c2_3, check_c2_1, check_c2_2, &
check_c2_3 ,Nx,Ny,Nz,DY, DYF,Lx,Ly,Lz, c2_norm)
 kc=c1_norm
IF (c2_norm .lt. c1_norm) THEN
 kc=c2_norm
END IF

Uo=Uo+eps*(v1o-kc*Uo)
Vo=Vo+eps*(v2o-kc*Vo)
Wo=Wo+eps*(v3o-kc*Wo)
THo=THo+((eps*T_ref**2)/Ri)*(stauo-kc*THo*Ri/(T_ref**2))

END SUBROUTINE update_initial_real_fields

SUBROUTINE Poisson_Adjoint_RHS (v1,v2,v3,stau,u_tot,v_tot,w_tot,TH_tot,THB, mu_tot, U_bar, V_bar, W_bar, &
kx, kz, NX, NY, NZ, DY, DYF, K_start, K_end, Re, A2, plan_bkd, plan_fwd, F_Adj_Vel_func)
  USE, INTRINSIC :: iso_c_binding
  USE Fourier_Spectral
  IMPLICIT NONE
  include 'fftw3.f03'
  INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
  INTEGER, INTENT(IN) :: NX, NY, NZ
  INTEGER, INTENT(IN) :: K_Start, K_End
  REAL(KIND=DP), INTENT(IN) :: Re, A2
  REAL(KIND=DP), DIMENSION(1:Nx/2+1), INTENT(IN) :: kx
  REAL(KIND=DP), DIMENSION(1:Nz), INTENT(IN) :: kz
  REAL(KIND=DP), DIMENSION(0:NY), INTENT(IN) :: DY, DYF
  REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1), INTENT(INOUT) :: u_tot, v_tot, w_tot, U_bar, V_bar, W_bar, &
 	mu_tot, v1, v2, v3, stau, TH_tot, THB
 COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(1:NX/2+1,1:NZ,0:NY+1), INTENT(OUT) :: F_Adj_Vel_func
 type(C_PTR), INTENT(IN) :: plan_bkd, plan_fwd

 REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1) :: u, v, w, TH_tot_p_THB, v_frac, V_bar_frac, v2_frac, &
  ux, uz, U_barx, U_barz, v_fracx, v_fracz, V_bar_fracx, V_bar_fracz, wx, wz, &
W_barx, W_barz, TH_tot_p_THBx, TH_tot_p_THBz, mu_totx, mu_totz, uxx, uzz, v_fracxx, v_fraczz, &
wxx, wzz, TH_tot_p_THBxx, TH_tot_p_THBzz, mu_totxx, mu_totzz, mu_totxz, v1x, v1z, v2_fracx, &
v2_fracz, v3x, v3z, staux, stauz, mu_totx_dbl_breve, mu_totz_dbl_breve, v_fracy, uy, wy, &
 TH_tot_p_THBy, mu_toty, mu_totxy, mu_totzy, V_bar_fracy, U_bary, W_bary, v2_fracy, v1y, v3y, &
stauy, uyy, wyy, v_fracyy, TH_tot_p_THByy, mu_totyy, &
u_dbl_breve, w_dbl_breve, TH_tot_p_THB_dbl_breve, &
mu_tot_dbl_breve, U_bar_dbl_breve, W_bar_dbl_breve, v1_dbl_breve, v3_dbl_breve, &
stau_dbl_breve,Adj_Vel_func
COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(1:NX/2+1,1:NZ,0:NY+1) :: F_u, F_w, F_U_bar, &
	F_W_bar, F_mu_tot, F_TH_tot_p_THB, F_v1, F_v3, F_stau, F_v_frac, F_V_bar_frac, F_v2_frac, &
F_ux, F_uz, F_U_barx, F_U_barz, F_v_fracx, F_v_fracz, F_V_bar_fracx, F_V_bar_fracz, F_wx, F_wz, &
 F_W_barx, F_W_barz, F_TH_tot_p_THBx, F_TH_tot_p_THBz, F_mu_totx, F_mu_totz, F_uxx, F_uzz, &
 F_v_fracxx, F_v_fraczz, F_wxx, F_wzz, F_TH_tot_p_THBxx, F_TH_tot_p_THBzz, F_mu_totxx, F_mu_totzz, &
 F_mu_totxz, F_v1x, F_v1z, F_v2_fracx, F_v2_fracz, F_v3x, F_v3z, F_staux, F_stauz
 COMPLEX(C_DOUBLE_COMPLEX), PARAMETER :: ii=(0.0_DP, 1.0_DP)
 INTEGER :: I, J, K
TH_tot_p_THB=TH_tot+THB
u=u_tot-U_bar
v=v_tot-V_bar
w=w_tot-W_bar
  CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, u, F_u )
  CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, w, F_w )
  CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, U_bar, F_U_bar )
  CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, W_bar, F_W_bar )
  CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, mu_tot, F_mu_tot )
  CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, TH_tot_p_THB, F_TH_tot_p_THB )

  CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, v1, F_v1 )
  CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, v3, F_v3 )
  CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, stau, F_stau )

  v_frac=0.0_DP
  V_bar_frac=0.0_DP
  v2_frac=0.0_DP
  FORALL (I=1:NX/2+1,J=1:NZ,K=0:NY)
    v_frac(I,J,K)=(v(I,J,K+1)+v(I,J,K))/2.0_DP
    V_bar_frac(I,J,K)=(V_bar(I,J,K+1)+V_bar(I,J,K))/2.0_DP

    v2_frac(I,J,K)=(v2(I,J,K+1)+v2(I,J,K))/2.0_DP
  END FORALL
  CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, v_frac, F_v_frac )
  CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, V_bar_frac, F_V_bar_frac )
  CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, v2_frac, F_v2_frac )

  FORALL  (I=1:NX, J=1:NZ, K=K_start:K_end)
    u_dbl_breve(I,J,K)=(u(I,J,K)*DYF(K-1)+u(I,J,K-1)*DYF(K))/(2.0_DP*DY(K-1))
    w_dbl_breve(I,J,K)=(w(I,J,K)*DYF(K-1)+w(I,J,K-1)*DYF(K))/(2.0_DP*DY(K-1))
    TH_tot_p_THB_dbl_breve(I,J,K)=(TH_tot_p_THB(I,J,K)*DYF(K-1)+TH_tot_p_THB(I,J,K-1)*DYF(K))/(2.0_DP*DY(K-1))
    mu_tot_dbl_breve(I,J,K)=(mu_tot(I,J,K)*DYF(K-1)+mu_tot(I,J,K-1)*DYF(K))/(2.0_DP*DY(K-1))

    U_bar_dbl_breve(I,J,K)=(U_bar(I,J,K)*DYF(K-1)+U_bar(I,J,K-1)*DYF(K))/(2.0_DP*DY(K-1))
    W_bar_dbl_breve(I,J,K)=(W_bar(I,J,K)*DYF(K-1)+W_bar(I,J,K-1)*DYF(K))/(2.0_DP*DY(K-1))

    v1_dbl_breve(I,J,K)=(v1(I,J,K)*DYF(K-1)+v1(I,J,K-1)*DYF(K))/(2.0_DP*DY(K-1))
    v3_dbl_breve(I,J,K)=(v3(I,J,K)*DYF(K-1)+v3(I,J,K-1)*DYF(K))/(2.0_DP*DY(K-1))
    stau_dbl_breve(I,J,K)=(stau(I,J,K)*DYF(K-1)+stau(I,J,K-1)*DYF(K))/(2.0_DP*DY(K-1))
  END FORALL

  FORALL (I=1:NX/2+1,J=1:NZ,K=0:NY+1)
    F_ux(I,J,K)=ii*kx(I)*F_u(I,J,K)
    F_uz(I,J,K)=ii*kz(J)*F_u(I,J,K)
    F_U_barx(I,J,K)=ii*kx(I)*F_U_bar(I,J,K)
    F_U_barz(I,J,K)=ii*kz(J)*F_U_bar(I,J,K)

    F_v_fracx(I,J,K)=ii*kx(I)*F_v_frac(I,J,K)
    F_v_fracz(I,J,K)=ii*kz(J)*F_v_frac(I,J,K)
    F_V_bar_fracx(I,J,K)=ii*kx(I)*V_bar_frac(I,J,K)
    F_V_bar_fracz(I,J,K)=ii*kz(J)*V_bar_frac(I,J,K)

    F_wx(I,J,K)=ii*kx(I)*F_w(I,J,K)
    F_wz(I,J,K)=ii*kz(J)*F_w(I,J,K)
    F_W_barx(I,J,K)=ii*kx(I)*F_W_bar(I,J,K)
    F_W_barz(I,J,K)=ii*kz(J)*F_W_bar(I,J,K)

    F_TH_tot_p_THBx(I,J,K)=ii*kx(I)*F_TH_tot_p_THB(I,J,K)
    F_TH_tot_p_THBz(I,J,K)=ii*kz(J)*F_TH_tot_p_THB(I,J,K)

    F_mu_totx(I,J,K)=ii*kx(I)*F_mu_tot(I,J,K)
    F_mu_totz(I,J,K)=ii*kz(J)*F_mu_tot(I,J,K)


    F_uxx(I,J,K)=(-1.0_DP)*(kx(I)**2)*F_u(I,J,K)
    F_uzz(I,J,K)=(-1.0_DP)*(kz(J)**2)*F_u(I,J,K)

    F_v_fracxx(I,J,K)=(-1.0_DP)*(kx(I)**2)*F_v_frac(I,J,K)
    F_v_fraczz(I,J,K)=(-1.0_DP)*(kz(J)**2)*F_v_frac(I,J,K)

    F_wxx(I,J,K)=(-1.0_DP)*(kx(I)**2)*F_w(I,J,K)
    F_wzz(I,J,K)=(-1.0_DP)*(kz(J)**2)*F_w(I,J,K)

    F_TH_tot_p_THBxx(I,J,K)=(-1.0_DP)*(kx(I)**2)*F_TH_tot_p_THB(I,J,K)
    F_TH_tot_p_THBzz(I,J,K)=(-1.0_DP)*(kz(J)**2)*F_TH_tot_p_THB(I,J,K)

    F_mu_totxx(I,J,K)=(-1.0_DP)*(kx(I)**2)*F_mu_tot(I,J,K)
    F_mu_totzz(I,J,K)=(-1.0_DP)*(kz(J)**2)*F_mu_tot(I,J,K)
    F_mu_totxz(I,J,K)=(-1.0_DP)*(kx(I)*kz(J))*F_mu_tot(I,J,K)

    F_v1x(I,J,K)=ii*kx(I)*F_v1(I,J,K)
    F_v1z(I,J,K)=ii*kz(J)*F_v1(I,J,K)
    F_v2_fracx(I,J,K)=ii*kx(I)*F_v2_frac(I,J,K)
    F_v2_fracz(I,J,K)=ii*kz(J)*F_v2_frac(I,J,K)
    F_v3x(I,J,K)=ii*kx(I)*F_v3(I,J,K)
    F_v3z(I,J,K)=ii*kz(J)*F_v3(I,J,K)

    F_staux(I,J,K)=ii*kx(I)*F_stau(I,J,K)
    F_stauz(I,J,K)=ii*kz(J)*F_stau(I,J,K)
  END FORALL

  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_ux, ux )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_uz, uz )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_U_barx, U_barx )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_U_barz, U_barz )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_v_fracx, v_fracx )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_v_fracz, v_fracz)
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_V_bar_fracx, V_bar_fracx )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_V_bar_fracz, V_bar_fracz )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_wx, wx )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_wz, wz)
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_W_barx, W_barx )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_W_barz, W_barz )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_TH_tot_p_THBx, TH_tot_p_THBx )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_TH_tot_p_THBz, TH_tot_p_THBz)
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_mu_totx, mu_totx )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_mu_totz, mu_totz )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_uxx, uxx )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_uzz, uzz)
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_v_fracxx, v_fracxx )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_v_fraczz, v_fraczz )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_wxx, wxx )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_wzz, wzz)
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_TH_tot_p_THBxx, TH_tot_p_THBxx )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_TH_tot_p_THBzz, TH_tot_p_THBzz )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_mu_totxx, mu_totxx )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_mu_totzz, mu_totzz )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_mu_totxz, mu_totxz )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_v1x, v1x)
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_v1z, v1z )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_v2_fracx, v2_fracx )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_v2_fracz, v2_fracz )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_v3x, v3x)
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_v3z, v3z )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_staux, staux )
  CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_stauz, stauz )

  FORALL  (I=1:NX, J=1:NZ, K=K_start:K_end)
    mu_totx_dbl_breve(I,J,K)=(mu_totx(I,J,K)*DYF(K-1)+mu_totx(I,J,K-1)*DYF(K))/(2.0_DP*DY(K-1))
    mu_totz_dbl_breve(I,J,K)=(mu_totz(I,J,K)*DYF(K-1)+mu_totz(I,J,K-1)*DYF(K))/(2.0_DP*DY(K-1))
  END FORALL

  FORALL  (I=1:NX, J=1:NZ, K=K_start:K_end)
    v_fracy(I,J,K)=(v(I,J,K+1)-v(I,J,K))/DYF(K)
    uy(I,J,K)=(u_dbl_breve(I,J,K+1)-u_dbl_breve(I,J,K))/DYF(K)
    wy(I,J,K)=(w_dbl_breve(I,J,K+1)-w_dbl_breve(I,J,K))/DYF(K)
    TH_tot_p_THBy(I,J,K)=(TH_tot_p_THB_dbl_breve(I,J,K+1)-TH_tot_p_THB_dbl_breve(I,J,K))/DYF(K)
    mu_toty(I,J,K)=(mu_tot_dbl_breve(I,J,K+1)-mu_tot_dbl_breve(I,J,K))/DYF(K)
    mu_totxy(I,J,K)=(mu_totx_dbl_breve(I,J,K+1)-mu_totx_dbl_breve(I,J,K))/DYF(K)
    mu_totzy(I,J,K)=(mu_totz_dbl_breve(I,J,K+1)-mu_totz_dbl_breve(I,J,K))/DYF(K)

    V_bar_fracy(I,J,K)=(V_bar(I,J,K+1)-V_bar(I,J,K))/DYF(K)
    U_bary(I,J,K)=(U_bar_dbl_breve(I,J,K+1)-U_bar_dbl_breve(I,J,K))/DYF(K)
    W_bary(I,J,K)=(W_bar_dbl_breve(I,J,K+1)-W_bar_dbl_breve(I,J,K))/DYF(K)

    v2_fracy(I,J,K)=(v2(I,J,K+1)-v2(I,J,K))/DYF(K)
    v1y(I,J,K)=(v1_dbl_breve(I,J,K+1)-v1_dbl_breve(I,J,K))/DYF(K)
    v3y(I,J,K)=(v3_dbl_breve(I,J,K+1)-v3_dbl_breve(I,J,K))/DYF(K)
    stauy(I,J,K)=(stau_dbl_breve(I,J,K+1)-stau_dbl_breve(I,J,K))/DYF(K)

    uyy(I,J,K)= ((u(I,J,K+1)-u(I,J,K))/DY(K)-(u(I,J,K)-u(I,J,K-1))/DY(K-1))/DYF(K)
    wyy(I,J,K)= ((w(I,J,K+1)-w(I,J,K))/DY(K)-(w(I,J,K)-w(I,J,K-1))/DY(K-1))/DYF(K)
    v_fracyy(I,J,K)= ((v_frac(I,J,K+1)-v_frac(I,J,K))/DY(K)-(v_frac(I,J,K)-v_frac(I,J,K-1))/DY(K-1))/DYF(K)
    TH_tot_p_THByy(I,J,K)= ((TH_tot_p_THB(I,J,K+1)-TH_tot_p_THB(I,J,K))/DY(K)&
    -(TH_tot_p_THB(I,J,K)-TH_tot_p_THB(I,J,K-1))/DY(K-1))/DYF(K)
    mu_totyy(I,J,K)= ((mu_tot(I,J,K+1)-mu_tot(I,J,K))/DY(K)-(mu_tot(I,J,K)-mu_tot(I,J,K-1))/DY(K-1))/DYF(K)
  END FORALL

Adj_Vel_func=(-1.0_DP)*(v1*uxx+ v1x*ux+ v2_frac*v_fracyy+ v2_fracy*v_fracy+ v3*wzz+ v3z*wz+ &
			v2_frac*v_fracxx+ v2_fracx*v_fracx+ v1*uyy+ v1y*uy+ v3*wxx+ v3x*wx+ &
			v1*uzz+ v1z*uz+ v3*wyy+ v3y*wy+ v2_frac*v_fraczz+ v2_fracz*v_fracz)

Adj_Vel_func=Adj_Vel_func + (-1.0_DP)*(&
		(U_barx+ux)*v1x+ (V_bar_fracy+v_fracy)*v2_fracy+ (W_barz+wz)*v3z+ &
		(U_bary+uy)*v2_fracx+ (V_bar_fracx+v_fracx)*v1y+ (W_barx+wx)*v1z+ &
		(U_barz+uz)*v3x+ (V_bar_fracz+v_fracz)*v3y+ (W_bary+wy)*v2_fracz)

Adj_Vel_func=Adj_Vel_func + &
		stau*(TH_tot_p_THBxx+TH_tot_p_THByy+TH_tot_p_THBzz)+&
		staux*TH_tot_p_THBx+stauy*TH_tot_p_THBy+stauz*TH_tot_p_THBz

Adj_Vel_func=Adj_Vel_func + A2/Re*(&
		mu_totxx*ux+ mu_totxy*v_fracx+ mu_totxz*wx+&
		mu_totxy*uy+ mu_totyy*v_fracy+ mu_totzy*wy+&
		mu_totxz*uz+ mu_totzy*v_fracz+ mu_totzz*wz+&
		mu_totx*(uxx+uyy+uzz)+&
		mu_toty*(v_fracxx+v_fracyy+v_fraczz)+&
		mu_totz*(wxx+wyy+wzz))

Adj_Vel_func = Adj_Vel_func - (2.0_DP/Re)*( &
mu_totxx*ux+mu_totyy*v_fracy+mu_totzz*wz+mu_totxy*(uy+v_fracx)+mu_totxz*(uz+wx)+mu_totzy*(wy+v_fracz)+&
mu_totx*(uxx+uyy+uzz)+mu_toty*(v_fracxx+v_fracyy+v_fraczz)+mu_totz*(wxx+wyy+wzz))

CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, Adj_Vel_func, F_Adj_Vel_func )

END SUBROUTINE Poisson_Adjoint_RHS

END MODULE Convergence_Check
