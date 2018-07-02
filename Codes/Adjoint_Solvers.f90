MODULE Adjoint_Solvers
! Make current step rk and the previous rkm1 instead of how it is presently: where current step is rkp1 and previous is rk
SUBROUTINE cubic_spline(F, alpha, delta_t, Nx, Ny, Nz)
IMPLICIT NONE
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: Nx, Ny, Nz
REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1, 4), INTENT(INOUT): F
REAL(KIND=DP), DIMENSION (3), INTENT(IN): alpha
REAL(KIND=DP), INTENT(IN): delta_t
REAL(KIND=DP): a, b, c, l1, l2, l3, n1, n2, beta
n1=alpha(3)*delta_t
n2=(alpha(2)+alpha(3))*delta_t
beta=delta_t*n1*n2*(delta_t-n1)*(delta_t-n2)*(n1-n2)
FORALL (I=1:NX, J=1:NZ, K=0:NZ+1)
l1=F(I,J,K,1)-F(I,J,K,2)
l2=F(I,J,K,1)-F(I,J,K,3)
l3=F(I,J,K,1)-F(I,J,K,4)
a=(l1*n2*delta_t*(delta_t-n2)-l2*n1*delta_t*(delta_t-n1)+l3*n1*n2*(n2-n1))/beta
b=(l2*n1*delta_t*(delta_t**2-n1**2)-l1*n2*delta_t*(delta_t**2-n2**2)&
  +l3*n1*n2*(n1**2-n2**2))/beta
c=(l1*(n2**2)*(delta_t**2)*(delta_t-n2)-l2*(n1**2)*(delta_t**2)*(delta_t-n1)&
  +l3*(n1**2)*(n2**2)*(n2-n1))/beta
F(I,J,K,2)= a*(alpha(1)*delta_t)**3 + b*(alpha(1)*delta_t)**2&
           +c*alpha(1)*delta_t + F(I,J,K,1)
F(I,J,K,3)= a*((alpha(1)+alpha(2))*delta_t)**3 &
           +b*((alpha(1)+alpha(2))*delta_t)**2 &
           +c*(alpha(1)+alpha(2))*delta_t + F(I,J,K,1)
END FORALL
END SUBROUTINE cubic_spline

SUBROUTINE Adjoint_uw_Boundary_Conditions(A, B, C, v, NX, NY, NZ, v_BC_Lower,&
                                        v_BC_Upper, v_wall_Lower, v_wall_Upper)
IMPLICIT NONE
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: Nx, Ny, Nz

END SUBROUTINE Adjoint_v_Boundary_Conditions

SUBROUTINE Adjoint_v_Boundary_Conditions(A, B, C, v2, NX, NY, NZ, v2_BC_Lower,&
                                      v2_BC_Upper, v2_wall_Lower, v2_wall_Upper)
IMPLICIT NONE
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: Nx, Ny, Nz

END SUBROUTINE Adjoint_v_Boundary_Conditions


SUBROUTINE Adjoint_scalar_Boundary_Conditions(A, B, C, stau, NX, NY, NZ, &
                stau_BC_Lower, stau_BC_Upper, stau_wall_Lower, stau_wall_Upper)
IMPLICIT NONE
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: Nx, Ny, Nz

END SUBROUTINE Adjoint_scalar_Boundary_Conditions


SUBROUTINE RK_Solver_Back( K_start, K_end, NX, NY, NZ, v1_BC_Lower, v1_BC_Upper, &
v2_BC_Lower, v2_BC_Upper, v3_BC_Lower, v3_BC_Upper, stau_BC_Lower, stau_BC_Upper, &
U_total_4, V_total_4, W_total_4,  U_bar_4, V_bar_4, W_bar_4, TH_bar_4, TH_total_4, &
THB, DY, DYF, gamma, zeta, alpha, delta_t, Lx, Lz, n1, n2, n3, Pr, Re, Ri, T_ref, &
A2, v1_wall_Lower, v1_wall_Upper, v2_wall_Lower, v2_wall_Upper, v3_wall_Lower, &
v3_wall_Upper, stau_wall_Lower, stau_wall_Upper, plan_bkd, plan_fwd, kx, kz, &
v1, v2, v3, Q, stau )
USE Remove_Divergence
USE Thomas_Matrix_Algorithm_real
USE, INTRINSIC :: iso_c_binding
USE Fourier_Spectral
IMPLICIT NONE
include 'fftw3.f03'
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER :: I, J, K, RK_step

INTEGER, INTENT(IN) :: K_start, K_end, NX, NY, NZ, v1_BC_Lower, v1_BC_Upper, &
v2_BC_Lower, v2_BC_Upper, v3_BC_Lower, v3_BC_Upper, stau_BC_Lower, stau_BC_Upper

REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1, 4), INTENT(IN):: U_total_4, &
V_total_4, W_total_4,  U_bar_4, V_bar_4, W_bar_4, TH_bar_4, TH_total_4

REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1), INTENT(IN):: THB

REAL(KIND=DP), DIMENSION (0:NY+1), INTENT(IN) :: DY, DYF

REAL(KIND=DP), DIMENSION (1:3), INTENT(IN) :: gamma, zeta, alpha

REAL(KIND=DP), INTENT (IN) :: delta_t, Lx, Lz, n1, n2, n3, Pr, Re, Ri, T_ref, &
A2, v1_wall_Lower, v1_wall_Upper, v2_wall_Lower, v2_wall_Upper, v3_wall_Lower, &
v3_wall_Upper, stau_wall_Lower, stau_wall_Upper

type(C_PTR), INTENT(IN) :: plan_bkd, plan_fwd

REAL(KIND=DP), DIMENSION(1:Nx/2+1), INTENT(IN) :: kx

REAL(KIND=DP), DIMENSION(1:Nz), INTENT(IN) :: kz

REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1), INTENT(OUT) :: v1, v2, v3, Q, stau

REAL(KIND=DP), DIMENSION (1:NX, 1:NZ, 0:NY+1) ::  U_bar, V_bar, W_bar, TH_bar, &
U_bar_rkp1, V_bar_rkp1, W_bar_rkp1, TH_bar_rkp1, &
u, v, w, TH, THpTHB, u_rkp1, v_rkp1, w_rkp1, TH_rkp1, THpTHB_rkp1, &
u_dbl_breve, w_dbl_breve, u_dbl_breve_rkp1, w_dbl_breve_rkp1, &
U_bar_dbl_breve, W_bar_dbl_breve, U_bar_dbl_breve_rkp1, W_bar_dbl_breve_rkp1, &
ux, vx, wx, uz, vz, wz, ux_rkp1, vx_rkp1, wx_rkp1, uz_rkp1, vz_rkp1, wz_rkp1, &
THBpTHBx, THBpTHBz, THBpTHBx_rkp1, THBpTHBz_rkp1, &
U_barx, V_barx, W_barx, U_barz, V_barz, W_barz, &
U_barx_rkp1, V_barx_rkp1, W_barx_rkp1, U_barz_rkp1, V_barz_rkp1, W_barz_rkp1,&
v1x, v2x, v3x, Qx, v1z, v2z, v3z, Qz, v1_dbl_breve, v3_dbl_breve, stau_dblbreve, &
mu, mu_bar, mutot_dbl_breve, mu_rkp1, mu_bar_rkp1, mutot_dbl_breve_rkp1, &
mu_T, mu_bar_T, mu_T_rkp1, mu_bar_T_rkp1, &
v1Utot_p_mutot_v1x, v1Wtot_p_mutot_v1z_v3x, v3Utot_p_mutot_v3x, v3Wtot_p_mutot_v3z, &
v2Utot_p_mutot_v2x, v2Wtot_p_mutot_v2z, &
tau_Utot, tau_Wtot, Exp_v1, Exp_v2, Exp_v3, Exp_stau, &
exp_spec_terms_w1_x, exp_spec_terms_w1_z, exp_spec_terms_w2_x, exp_spec_terms_w2_z, &
 exp_spec_terms_w3_x, exp_spec_terms_w3_z, &
Cranck_Exp_v1, Cranck_Exp_v2, Cranck_Exp_v3, Cranck_Exp_stau, Cranck_Exp_v2_temp, &
Exp_v1_m1, Exp_v2_m1, Exp_v3_m1, Exp_stau_m1, A, B, C


COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(1:NX/2+1,1:NZ,0:NY+1) :: F_u_rkp1, F_v_rkp1, &
F_w_rkp1, F_TH_rkp1,F_TH, F_THpTHB_rkp1, F_U_bar_rkp1, F_V_bar_rkp1, F_W_bar_rkp1, &
F_ux_rkp1, F_vx_rkp1, F_wx_rkp1, F_uz_rkp1, F_vz_rkp1, F_wz_rkp1, &
F_THBpTHBx_rkp1, F_THBpTHBz_rkp1, F_U_barx_rkp1, F_U_barz_rkp1, F_V_barx_rkp1, &
F_V_barz_rkp1, F_W_barx_rkp1, F_W_barz_rkp1, &
F_v1, F_v2, F_v3, F_Q, F_stau, F_v1x, F_v2x, F_v3x, F_Qx, F_v1z, F_v2z, F_v3z, F_Qz, &
F_v1Utot_p_mutot_v1x, F_v1Wtot_p_mutot_v1z_v3x, F_v3Utot_p_mutot_v3x, F_v3Wtot_p_mutot_v3z, &
F_v2Utot_p_mutot_v2x, F_v2Wtot_p_mutot_v2z, F_tau_Utot, F_tau_Wtot, F_Exp_v1, F_Exp_v3, F_Exp_v2, &
F_Exp_stau, F_exp_spec_terms_w1_x, F_exp_spec_terms_w1_z, F_exp_spec_terms_w2_x, F_exp_spec_terms_w2_z&
, F_exp_spec_terms_w3_x, F_exp_spec_terms_w3_z, &
F_Cranck_Exp_v1, F_Cranck_Exp_v2, F_Cranck_Exp_v3, F_Cranck_Exp_stau,


CALL cubic_spline(U_total_4, alpha, delta_t, Nx, Ny, Nz)
CALL cubic_spline(V_total_4, alpha, delta_t, Nx, Ny, Nz)
CALL cubic_spline(W_total_4, alpha, delta_t, Nx, Ny, Nz)
CALL cubic_spline(U_bar_4, alpha, delta_t, Nx, Ny, Nz)
CALL cubic_spline(V_bar_4, alpha, delta_t, Nx, Ny, Nz)
CALL cubic_spline(W_bar_4, alpha, delta_t, Nx, Ny, Nz)
CALL cubic_spline(TH_bar_4, alpha, delta_t, Nx, Ny, Nz)
CALL cubic_spline(TH_total_4, alpha, delta_t, Nx, Ny, Nz)
CALL cubic_spline(TH_bar_4, alpha, delta_t, Nx, Ny, Nz)
CALL cubic_spline(TH_total_4, alpha, delta_t, Nx, Ny, Nz)

Exp_v1_m1 =0.0_DP
Exp_v2_m1 =0.0_DP
Exp_v3_m1 =0.0_DP
Exp_stau_m1 =0.0_DP

DO RK_step = 1,3

U_bar=U_bar_4(:,:,:,RK_step)
U_bar_rkp1=U_bar_4(:,:,:,RK_step+1)
V_bar=V_bar_4(:,:,:,RK_step)
V_bar_rkp1=V_bar_4(:,:,:,RK_step+1)
W_bar=W_bar_4(:,:,:,RK_step)
W_bar_rkp1=W_bar_4(:,:,:,RK_step+1)
TH_bar=TH_bar_4(:,:,:,RK_step)
TH_bar_rkp1=TH_bar_4(:,:,:,RK_step+1)

u=U_total_4(:,:,:,RK_step)-U_bar
u_rkp1=U_total_4(:,:,:,RK_step+1)-U_bar_rkp1
v=V_total_4(:,:,:,RK_step)-V_bar
v_rkp1=V_total_4(:,:,:,RK_step+1)-V_bar_rkp1
w=W_total_4(:,:,:,RK_step)-W_bar
w_rkp1=W_total_4(:,:,:,RK_step+1)-W_bar_rkp1
TH=TH_total_4(:,:,:,RK_step)-TH_bar
TH_rkp1=TH_total_4(:,:,:,RK_step+1)-TH_bar_rkp1
THpTHB=TH_total_4(:,:,:,RK_step)+THB
THpTHB_rkp1=TH_total_4(:,:,:,RK_step+1)+THB

!!!! call a subroutine to find mu, mubar, mu_rkp1 and mubar_rkp1, mu_bar_T, and mu_bar_T_rkp1 as a function of TH, THbar, THrkp1 and THbarrkp1

stau_dblbreve=0.0_DP
FORALL  (I=1:NX, J=1:NZ, K=K_start:K_end)
stau_dblbreve(I,J,K)=(stau(I,J,K)*DYF(K-1)+stau(I,J,K-1)*DYF(K))/(2*DY(K-1))
mutot_dbl_breve(I,J,K)=((mu(I,J,K)+mu_bar(I,J,K))*DYF(K-1)+(mu(I,J,K-1)+mu_bar(I,J,K-1))*DYF(K))/(2*DY(K-1))
u_dbl_breve(I,J,K)=(u(I,J,K)*DYF(K-1)+u(I,J,K-1)*DYF(K))/(2*DY(K-1))
U_bar_dbl_breve(I,J,K)=(U_bar(I,J,K)*DYF(K-1)+U_bar(I,J,K-1)*DYF(K))/(2*DY(K-1))
w_dbl_breve(I,J,K)=(w(I,J,K)*DYF(K-1)+w(I,J,K-1)*DYF(K))/(2*DY(K-1))
W_bar_dbl_breve(I,J,K)=(W_bar(I,J,K)*DYF(K-1)+W_bar(I,J,K-1)*DYF(K))/(2*DY(K-1))
u_dbl_breve_rkp1(I,J,K)=(u_rkp1(I,J,K)*DYF(K-1)+u_rkp1(I,J,K-1)*DYF(K))/(2*DY(K-1))
U_bar_dbl_breve_rkp1(I,J,K)=(U_bar_rkp1(I,J,K)*DYF(K-1)+U_bar_rkp1(I,J,K-1)*DYF(K))/(2*DY(K-1))
w_dbl_breve_rkp1(I,J,K)=(w_rkp1(I,J,K)*DYF(K-1)+w_rkp1(I,J,K-1)*DYF(K))/(2*DY(K-1))
W_bar_dbl_breve_rkp1(I,J,K)=(W_bar_rkp1(I,J,K)*DYF(K-1)+W_bar_rkp1(I,J,K-1)*DYF(K))/(2*DY(K-1))
mutot_dbl_breve_rkp1(I,J,K)=((mu_rkp1(I,J,K)+mu_bar_rkp1(I,J,K))*DYF(K-1)+(mu_rkp1(I,J,K-1)+mu_bar_rkp1(I,J,K-1))*DYF(K))/(2*DY(K-1))
v1_dbl_breve(I,J,K)=(v1(I,J,K)*DYF(K-1)+v1(I,J,K-1)*DYF(K))/(2*DY(K-1))
v3_dbl_breve(I,J,K)=(v3(I,J,K)*DYF(K-1)+v3(I,J,K-1)*DYF(K))/(2*DY(K-1))


END FORALL

CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, u_rkp1, F_u_rkp1 )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, v_rkp1 ,F_v_rkp1 )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, w_rkp1, F_w_rkp1 )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, TH_rkp1,F_TH_rkp1 )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, TH,F_TH )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, THpTHB_rkp1, F_THpTHB_rkp1 )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, U_bar_rkp1, F_U_bar_rkp1 )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, V_bar_rkp1, F_V_bar_rkp1 )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, W_bar_rkp1, F_W_bar_rkp1 )


CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, v1, F_v1 )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, v2, F_v2 )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, v3, F_v3 )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, Q, F_Q )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, stau, F_stau )

FORALL (I=1:NX/2+1,J=1:NZ,K=0:NY+1)
F_ux_rkp1(I,J,K)=ii*kx(I)*F_u_rkp1(I,J,K)
F_vx_rkp1(I,J,K)=ii*kx(I)*F_v_rkp1(I,J,K)
F_wx_rkp1(I,J,K)=ii*kx(I)*F_w_rkp1(I,J,K)
F_uz_rkp1(I,J,K)=ii*kz(J)*F_u_rkp1(I,J,K)
F_vz_rkp1(I,J,K)=ii*kz(J)*F_v_rkp1(I,J,K)
F_wz_rkp1(I,J,K)=ii*kz(J)*F_w_rkp1(I,J,K)
F_THBpTHBx_rkp1(I,J,K)=ii*kx(I)*F_THpTHB_rkp1(I,J,K)
F_THBpTHBz_rkp1(I,J,K)=ii*kz(J)*F_THpTHB_rkp1(I,J,K)
F_U_barx_rkp1(I,J,K)=ii*kx(I)*F_U_bar_rkp1(I,J,K)
F_U_barz_rkp1(I,J,K)=ii*kz(J)*F_U_bar_rkp1(I,J,K)
F_V_barx_rkp1(I,J,K)=ii*kx(I)*F_V_bar_rkp1(I,J,K)
F_V_barz_rkp1(I,J,K)=ii*kz(J)*F_V_bar_rkp1(I,J,K)
F_W_barx_rkp1(I,J,K)=ii*kx(I)*F_W_bar_rkp1(I,J,K)
F_W_barz_rkp1(I,J,K)=ii*kz(J)*F_W_bar_rkp1(I,J,K)

F_v1x(I,J,K)=ii*kx(I)*F_v1(I,J,K)
F_v2x(I,J,K)=ii*kx(I)*F_v2(I,J,K)
F_v3x(I,J,K)=ii*kx(I)*F_v3(I,J,K)
F_Qx(I,J,K) =ii*kx(I)*F_Q(I,J,K)
F_v1z(I,J,K)=ii*kz(J)*F_v1(I,J,K)
F_v2z(I,J,K)=ii*kz(J)*F_v2(I,J,K)
F_v3z(I,J,K)=ii*kz(J)*F_v3(I,J,K)
F_Qz(I,J,K) =ii*kz(J)*F_Q(I,J,K)

END FORALL

ux=ux_rkp1
vx=vx_rkp1
wx=wx_rkp1
uz=uz_rkp1
vz=vz_rkp1
wz=wz_rkp1
THBpTHBx=THBpTHBx_rkp1
THBpTHBz=THBpTHBz_rkp1
U_barx=U_barx_rkp1
V_barx=V_barx_rkp1
W_barx=W_barx_rkp1
U_barz=U_barz_rkp1
V_barz=V_barz_rkp1
W_barz=W_barz_rkp1

CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_ux_rkp1, ux_rkp1 )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_vx_rkp1, vx_rkp1 )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_wx_rkp1, wx_rkp1 )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_uz_rkp1, uz_rkp1 )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_vz_rkp1, vz_rkp1 )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_wz_rkp1, wz_rkp1 )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_THBpTHBx_rkp1, THBpTHBx_rkp1 )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_THBpTHBz_rkp1, THBpTHBz_rkp1 )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_U_barx_rkp1, U_barx_rkp1 )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_U_barz_rkp1, U_barz_rkp1 )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_V_barx_rkp1, V_barx_rkp1 )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_V_barz_rkp1, V_barz_rkp1 )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_W_barx_rkp1, W_barx_rkp1 )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_W_barz_rkp1, W_barz_rkp1 )



CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_v1x, v1x )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_v2x, v2x )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_v3x, v3x )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Qx, Qx )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_v1z, v1z )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_v2z, v2z )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_v3z, v3z )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Qz, Qz )


v1Utot_p_mutot_v1x    =v1*(U_bar+u)+2.0_DP/Re*(mu+mu_bar)*v1x
v1Wtot_p_mutot_v1z_v3x=v1*(W_bar+w)+1.0_DP/Re*(mu+mu_bar)*(v1z+v3x)
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, v1Utot_p_mutot_v1x, F_v1Utot_p_mutot_v1x )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, v1Wtot_p_mutot_v1z_v3x, F_v1Wtot_p_mutot_v1z_v3x )

v3Utot_p_mutot_v3x    =v3*(U_bar+u)+1.0_DP/Re*(mu+mu_bar)*v3x
v3Wtot_p_mutot_v3z    =v3*(W_bar+w)+2.0_DP/Re*(mu+mu_bar)*v3z
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, v3Utot_p_mutot_v3x, F_v3Utot_p_mutot_v3x )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, v3Wtot_p_mutot_v3z, F_v3Wtot_p_mutot_v3z )

v2Utot_p_mutot_v2x    =v2*(U_bar+u)+1.0_DP/Re*(mu+mu_bar)*v2x
v2Wtot_p_mutot_v2z    =v2*(W_bar+w)+1.0_DP/Re*(mu+mu_bar)*v2z
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, v2Utot_p_mutot_v2x, F_v2Utot_p_mutot_v2x )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, v2Wtot_p_mutot_v2z, F_v2Wtot_p_mutot_v2z )

tau_Utot=stau*(u+U_bar)
tau_Wtot=stau*(w+W_bar)
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, tau_Utot, F_tau_Utot )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, tau_Wtot, F_tau_Wtot )

FORALL (I=1:NX/2+1,J=1:NZ,K=0:NY+1)
F_Exp_v1(I,J,K)=ii*(kx(I)*F_v1Utot_p_mutot_v1x(I,J,K)+kz(J)*F_v1Wtot_p_mutot_v1z_v3x(I,J,K))
F_Exp_v3(I,J,K)=ii*(kx(I)*F_v3Utot_p_mutot_v3x(I,J,K)+kz(J)*F_v3Wtot_p_mutot_v3z(I,J,K))
F_Exp_v2(I,J,K)=ii*(kx(I)*F_v2Utot_p_mutot_v2x(I,J,K)+kz(J)*F_v2Wtot_p_mutot_v2z(I,J,K))
F_Exp_stau(I,J,K)=ii*(kx(I)*F_tau_Utot(I,J,K)+kz(J)*F_tau_Wtot(I,J,K))-1.0_DP/(Re*Pr)*(kx(I)**2+kz(J)**2)*F_stau(I,J,K)
END FORALL
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Exp_v1, Exp_v1 )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Exp_v3, Exp_v3 )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Exp_v2, Exp_v2 )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Exp_stau, Exp_stau )

Exp_v1 = Exp_v1 + v1*ux + v2*vx + v3*wx -stau*THBpTHBx
Exp_v3 = Exp_v3 + v1*uz + v2*vz + v3*wz -stau*THBpTHBz

! mutot_dbl_breve is double breve of mu + mu_bar
FORALL (I=1:NX, J=1:NZ, K=1:NY+1)
Exp_v1(I,J,K)=Exp_v1(I,J,K)+1.0_DP/(Re*DYF(K))*(mutot_dbl_breve(I,J,K+1)*v2x(I,J,K+1)-mutot_dbl_breve(I,J,K)*v2x(I,J,K))
Exp_v3(I,J,K)=Exp_v3(I,J,K)+1.0_DP/(Re*DYF(K))*(mutot_dbl_breve(I,J,K+1)*v2z(I,J,K+1)-mutot_dbl_breve(I,J,K)*v2z(I,J,K))
Exp_v2(I,J,K)=Exp_v2(I,J,K)-stau_dblbreve(I,J,K)*(THpTHB(I,J,K)-THpTHB(I,J,K-1))/(DY(K-1))
END FORALL

exp_spec_terms_w1_x=-1.0_DP*A2/Re*((mu+mu_bar)*ux+(mu_rkp1+mu_bar_rkp1)*ux_rkp1)
exp_spec_terms_w1_z=-1.0_DP*A2/Re*((mu+mu_bar)*uz+(mu_rkp1+mu_bar_rkp1)*uz_rkp1)
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, exp_spec_terms_w1_x, F_exp_spec_terms_w1_x )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, exp_spec_terms_w1_z, F_exp_spec_terms_w1_z )

FORALL (I=1:NX/2+1,J=1:NZ,K=0:NY+1)

F_Cranck_Exp_v1(I,J,K)=ii*(kx(I)*F_exp_spec_terms_w1_x(I,J,K)+kz(J)*F_exp_spec_terms_w1_z(I,J,K))

F_Cranck_Exp_stau(I,J,K)=2.0_DP*A2*Ri/(Re*Pr*(T_ref**2))*(kx(I)**2+kz(J)**2)*(F_TH(I,J,K)+F_TH_rkp1(I,J,K))
END FORALL

CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Cranck_Exp_v1, Cranck_Exp_v1 )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Cranck_Exp_stau, Cranck_Exp_stau )


FORALL (I=1:NX, J=1:NZ, K=1:NY+1)
! mutot_dbl_breve is double breve of mu + mu_bar
Cranck_Exp_v1(I,J,K)=Cranck_Exp_v1(I,J,K)+ &
  ( (v1(I,J,K+1)+v1(I,J,K)) * (V_bar(I,J,K+1)+v(I,J,K+1)) &
   -(v1(I,J,K)+v1(I,J,K-1)) * (V_bar(I,J,K)  +v(I,J,K  )) )/(2.0_DP*DYF(K)) &
   +1.0_DP/(Re*DYF(K))*(mutot_dbl_breve(I,J,K+1)*(v1(I,J,K+1)-v1(I,J,K))/DY(K) &
   - mutot_dbl_breve(I,J,K)*(v1(I,J,K)-v1(I,J,K-1))/DY(K-1)&
   -1.0_DP*A2*(mutot_dbl_breve(I,J,K+1)*(u(I,J,K+1)-u(I,J,K))/DY(K) &
   - mutot_dbl_breve(I,J,K)*(u(I,J,K)-u(I,J,K-1))/DY(K-1))&
   -1.0_DP*A2*(mutot_dbl_breve_rkp1(I,J,K+1)*(u_rkp1(I,J,K+1)-u_rkp1(I,J,K))/DY(K) &
   - mutot_dbl_breve_rkp1(I,J,K)*(u_rkp1(I,J,K)-u_rkp1(I,J,K-1))/DY(K-1)))

Cranck_Exp_v2(I,J,K) = (v1(I,J,K)*DYF(K)+v1(I,J,K-1)*DYF(K-1))*(u(I,J,K)-u(I,J,K-1))/(2.0_DP*(DY(K-1)**2)) &
                      +(v3(I,J,K)*DYF(K)+v3(I,J,K-1)*DYF(K-1))*(w(I,J,K)-w(I,J,K-1))/(2.0_DP*(DY(K-1)**2)) &
                      +(v2(I,J,K)*(v(I,J,K+1)-v(I,J,K)))/DY(K-1) &
                      +((v2(I,J,K+1)+v2(I,J,K))*(V_bar(I,J,K+1)+v(I,J,K+1)+V_bar(I,J,K)+v(I,J,K)) &
                      - (v2(I,J,K)+v2(I,J,K-1))*(V_bar(I,J,K)+v(I,J,K)+V_bar(I,J,K-1)+v(I,J,K-1)))/(4.0_DP*DY(K-1)) &
                      +(2.0_DP/Re)*((mu(I,J,K)+mu_bar(I,J,K))*(v2(I,J,K+1)-v2(I,J,K))/DYF(K)&
                                  - (mu(I,J,K-1)+mu_bar(I,J,K-1))*(v2(I,J,K)-v2(I,J,K-1))/DYF(K-1))/DY(K-1) &
                      -(A2/Re)*((mu(I,J,K)+mu_bar(I,J,K))*(v(I,J,K+1)-v(I,J,K))/DYF(K)&
                                 - (mu(I,J,K-1)+mu_bar(I,J,K-1))*(v(I,J,K)-v(I,J,K-1))/DYF(K-1))/DY(K-1) &
                      -(A2/Re)*((mu_rkp1(I,J,K)+mu_bar_rkp1(I,J,K))*(v_rkp1(I,J,K+1)-v_rkp1(I,J,K))/DYF(K)&
                                 - (mu_rkp1(I,J,K-1)+mu_bar_rkp1(I,J,K-1))*(v_rkp1(I,J,K)-v_rkp1(I,J,K-1))/DYF(K-1))/DY(K-1)

exp_spec_terms_w2_x(I,J,K) =1.0_DP/Re*mutot_dbl_breve(I,J,K)*(v1(I,J,K+1)-v1(I,J,K-1))/(2.0_DP*DY(K-1)) &
                            -A2/Re*(mutot_dbl_breve(I,J,K)*vx(I,J,K)+mutot_dbl_breve_rkp1(I,J,K)*vx_rkp1(I,J,K))
exp_spec_terms_w2_z(I,J,K) =1.0_DP/Re*mutot_dbl_breve(I,J,K)*(v3(I,J,K+1)-v3(I,J,K-1))/(2.0_DP*DY(K-1)) &
                            -A2/Re*(mutot_dbl_breve(I,J,K)*vz(I,J,K)+mutot_dbl_breve_rkp1(I,J,K)*vz_rkp1(I,J,K))


Cranck_Exp_stau(I,J,K)= Cranck_Exp_stau(I,J,K)&
&
+(stau_dblbreve(I,J,K+1)*(V_bar(I,J,K+1)+v(I,J,K+1))&
                  -stau_dblbreve(I,J,K)*(V_bar(I,J,K)+v(I,J,K)))/DYF(K) &
&
+ (1.0_DP/(Re*Pr))*((stau(I,J,K+1)-stau(I,J,K))/DY(K)&
                    -(stau(I,J,K)-stau(I,J,K-1))/DY(K-1))/DYF(K) &
&
-2.0_DP*A2*Ri/(Re*Pr*(T_ref**2))*&
((TH(I,J,K+1)-TH(I,J,K))/DY(K)-(TH(I,J,K)-TH(I,J,K-1))/DY(K-1)&
&
+(TH_rkp1(I,J,K+1)-TH_rkp1(I,J,K))/DY(K)-(TH_rkp1(I,J,K)-TH_rkp1(I,J,K-1))/DY(K-1))/DYF(K) &
&
+A2*(mu_T(I,J,K)+mu_bar_T(I,J,K))/Re*(ux(I,J,K)**2+((u_dbl_breve(I,J,K+1)-u_dbl_breve(I,J,K))/DYF(K))**2+uz(I,J,K)**2&
+vx(I,J,K)**2+(v(I,J,K+1)-v(I,J,K))/DYF(K))**2+vz(I,J,K)**2)&
+wx(I,J,K)**2+(w_dbl_breve(I,J,K+1)-w_dbl_breve(I,J,K))/DYF(K))**2+wz(I,J,K)**2)&
&
+A2*(mu_T_rkp1(I,J,K)+mu_bar_T_rkp1(I,J,K))/Re*(ux_rkp1(I,J,K)**2&
+((u_dbl_breve_rkp1(I,J,K+1)-u_dbl_breve_rkp1(I,J,K))/DYF(K))**2+uz_rkp1(I,J,K)**2&
+vx_rkp1(I,J,K)**2+(v_rkp1(I,J,K+1)-v_rkp1(I,J,K))/DYF(K))**2+vz_rkp1(I,J,K)**2)&
+wx_rkp1(I,J,K)**2+((w_dbl_breve_rkp1(I,J,K+1)-w_dbl_breve_rkp1(I,J,K))/DYF(K))**2+wz_rkp1(I,J,K)**2)&
&
-2.0_DP/Re(&
(mu_T(I,J,K)*(ux(I,J,K)+U_barx(I,J,K))+mu_bar_T(I,J,K)*ux(I,J,K))*v1x(I,J,K) + &
&
(mu_T(I,J,K)*((u_dbl_breve(I,J,K+1)+U_bar_dbl_breve(I,J,K+1)&
-u_dbl_breve(I,J,K)-U_bar_dbl_breve(I,J,K))/(2.0_DP*DYF(K))&
+0.25_DP*(vx(I,J,K+1)+vx(I,J,K)+V_barx(I,J,K+1)+V_barx(I,J,K)))+&
mu_bar_T(I,J,K)*((u_dbl_breve(I,J,K+1)-u_dbl_breve(I,J,K))/(2.0_DP*DYF(K))&
+0.25_DP*(vx(I,J,K+1)+vx(I,J,K))))*(0.5_DP*(v2x(I,J,K+1)+v2x(I,J,K))&
+(v1_dbl_breve(I,J,K+1)-v1_dbl_breve(I,J,K))/DYF(K) )  +&
&
((mu_T(I,J,K)*(uz(I,J,K)+U_barz(I,J,K)+wx(I,J,K)+W_barx(I,J,K))/2.0_DP+&
  mu_bar_T(I,J,K)*(uz(I,J,K)+wx(I,J,K))/2.0_DP)*(v3x(I,J,K)+v1z(I,J,K)))+&
&
(mu_T(I,J,K)*(wz(I,J,K)+W_barz(I,J,K))+mu_bar_T(I,J,K)*wz(I,J,K))*v3z(I,J,K) + &
&
(mu_T(I,J,K)*((w_dbl_breve(I,J,K+1)+W_bar_dbl_breve(I,J,K+1)&
-w_dbl_breve(I,J,K)-W_bar_dbl_breve(I,J,K))/(2.0_DP*DYF(K))&
+0.25_DP*(vz(I,J,K+1)+vz(I,J,K)+V_barz(I,J,K+1)+V_barz(I,J,K)))+&
mu_bar_T(I,J,K)*((w_dbl_breve(I,J,K+1)-w_dbl_breve(I,J,K))/(2.0_DP*DYF(K))&
+0.25_DP*(vz(I,J,K+1)+vz(I,J,K))))*(0.5_DP*(v2z(I,J,K+1)+v2z(I,J,K))&
+(v3_dbl_breve(I,J,K+1)-v3_dbl_breve(I,J,K))/DYF(K) )  +&
&
(mu_T(I,J,K)*(v(I,J,K+1)+V_bar(I,J,K+1)-v(I,J,K)-V_bar(I,J,K))/DYF(K)+&
mu_bar_T(I,J,K)*(v(I,J,K+1)-v(I,J,K))/DYF(K))*((v2(I,J,K+1)-v2(I,J,K))/DYF(K))&
)&
&
Ri*(n1*v1(I,J,K)+n2*v2(I,J,K)+n3*v3(I,J,K))

END FORALL
Cranck_Exp_v1=Cranck_Exp_v1+2.0_DP*Qx

exp_spec_terms_w3_x= -1.0_DP/Re*((mu+mu_bar)*(A2*wx-v1z)+A2*(mu_rkp1+mu_bar_rkp1)*wx_rkp1)
exp_spec_terms_w3_z= -1.0_DP*A2/Re*((mu+mu_bar)*wz+(mu_rkp1+mu_bar_rkp1)*wz_rkp1)
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, exp_spec_terms_w3_z, F_exp_spec_terms_w3_z )


!! ------------------- Finishing off the adjoint u equation ------------------
! IMPLICIT PART IS NOW ADDED
A=0.0_DP
B=0.0_DP
C=0.0_DP
DO J = 1, NZ
  DO I = 1, NX
    DO K = K_Start, K_End
    A(I,J,K) =  (V_bar_rkp1(I,J,K)+v_rkp1(I,J,K))/(2.0_DP*DYF(K))&
               - mutot_dbl_breve_rkp1(I,J,K)/(Re*DYF(K)*DY(K-1))
    C(I,J,K) = -(V_bar_rkp1(I,J,K+1)+v_rkp1(I,J,K+1))/(2.0_DP*DYF(K))&
               - mutot_dbl_breve_rkp1(I,J,K+1)/(Re*DYF(K)*DY(K))
    B(I,J,K) =1.0_DP - (V_bar_rkp1(I,J,K+1)+v_rkp1(I,J,K+1))/(2.0_DP*DYF(K))&
              + (V_bar_rkp1(I,J,K)+v_rkp1(I,J,K))/(2.0_DP*DYF(K)) &
              + mutot_dbl_breve_rkp1(I,J,K)/(Re*DYF(K)*DY(K-1)) &
              + mutot_dbl_breve_rkp1(I,J,K+1)/(Re*DYF(K)*DY(K))
END DO
END DO
END DO

v1 = v1+ delta_t*(gamma(RK_step)*Exp_v1 + zeta(RK_step)*Exp_v1_m1 + &
                    (alpha(RK_step)/2.0_DP)*Cranck_Exp_v1)
Exp_v1_m1 = Exp_v1

CALL Adjoint_uw_Boundary_Conditions(A,B,C,v1,NX,NY,NZ,v1_BC_Lower,v1_BC_Upper,v1_wall_Lower,v1_wall_Upper)
CALL Thomas_Matrix_Algorithm_real(A,B,C,v1,NX,NY,NZ)

CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, v1, F_v1 )
FORALL (I=1:NX/2+1,J=1:NZ,K=0:NY+1)
F_v1x(I,J,K)=ii*kx(I)*F_v1(I,J,K)
F_v1z(I,J,K)=ii*kz(J)*F_v1(I,J,K)
END FORALL
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_v1x, v1x )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_v1z, v1z )

exp_spec_terms_w3_x=exp_spec_terms_w3_x+1.0/Re*(mu_rkp1+mu_bar_rkp1)*v1z
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, exp_spec_terms_w3_x, F_exp_spec_terms_w3_x )
FORALL (I=1:NX/2+1,J=1:NZ,K=0:NY+1)
F_Cranck_Exp_v3(I,J,K)=ii*(kx(I)*F_exp_spec_terms_w3_x(I,J,K)+kz(J)*F_exp_spec_terms_w3_z(I,J,K))
END FORALL
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Cranck_Exp_v3, Cranck_Exp_v3 )

FORALL (I=1:NX, J=1:NZ, K=1:NY+1)
! mutot_dbl_breve is double breve of mu + mu_bar
Cranck_Exp_v3(I,J,K)=Cranck_Exp_v3(I,J,K)+ &
( (v3(I,J,K+1)+v3(I,J,K)) * (V_bar(I,J,K+1)+v(I,J,K+1)) &
 -(v3(I,J,K)+v3(I,J,K-1)) * (V_bar(I,J,K)  +v(I,J,K  )) )/(2.0_DP*DYF(K)) &
 +1.0_DP/(Re*DYF(K))*(mutot_dbl_breve(I,J,K+1)*(v3(I,J,K+1)-v3(I,J,K))/DY(K) &
 - mutot_dbl_breve(I,J,K)*(v3(I,J,K)-v3(I,J,K-1))/DY(K-1)&
 -1.0_DP*A2*(mutot_dbl_breve(I,J,K+1)*(w(I,J,K+1)-w(I,J,K))/DY(K) &
 - mutot_dbl_breve(I,J,K)*(w(I,J,K)-w(I,J,K-1))/DY(K-1))&
 -1.0_DP*A2*(mutot_dbl_breve_rkp1(I,J,K+1)*(w_rkp1(I,J,K+1)-w_rkp1(I,J,K))/DY(K) &
 - mutot_dbl_breve_rkp1(I,J,K)*(w_rkp1(I,J,K)-w_rkp1(I,J,K-1))/DY(K-1)))
END FORALL
Cranck_Exp_v3=Cranck_Exp_v3+2.0_DP*Qz

!! ------------------- Finishing off the adjoint w equation ------------------
! IMPLICIT PART IS NOW ADDED
A=0.0_DP
B=0.0_DP
C=0.0_DP
DO J = 1, NZ
  DO I = 1, NX
    DO K = K_Start, K_End
A(I,J,K) =  (V_bar_rkp1(I,J,K)+v_rkp1(I,J,K))/(2.0_DP*DYF(K))&
           - mutot_dbl_breve_rkp1(I,J,K)/(Re*DYF(K)*DY(K-1))
C(I,J,K) = -(V_bar_rkp1(I,J,K+1)+v_rkp1(I,J,K+1))/(2.0_DP*DYF(K))&
           - mutot_dbl_breve_rkp1(I,J,K+1)/(Re*DYF(K)*DY(K))
B(I,J,K) =1.0_DP - (V_bar_rkp1(I,J,K+1)+v_rkp1(I,J,K+1))/(2.0_DP*DYF(K))&
          + (V_bar_rkp1(I,J,K)+v_rkp1(I,J,K))/(2.0_DP*DYF(K)) &
          + mutot_dbl_breve_rkp1(I,J,K)/(Re*DYF(K)*DY(K-1)) &
          + mutot_dbl_breve_rkp1(I,J,K+1)/(Re*DYF(K)*DY(K))
END DO
END DO
END DO

v3 = v3+ delta_t*(gamma(RK_step)*Exp_v3 + zeta(RK_step)*Exp_v3_m1 + &
                    (alpha(RK_step)/2.0_DP)*Cranck_Exp_v3)
Exp_v3_m1 = Exp_v3

CALL Adjoint_uw_Boundary_Conditions(A,B,C,v3,NX,NY,NZ,v3_BC_Lower,v3_BC_Upper,v3_wall_Lower,v3_wall_Upper)
CALL Thomas_Matrix_Algorithm_real(A,B,C,v3,NX,NY,NZ)

! ------


FORALL (I=1:NX, J=1:NZ, K=1:NY+1)

Cranck_Exp_v2(I,J,K) = Cranck_Exp_v2(I,J,K) &
                      +(v1(I,J,K)*DYF(K)+v1(I,J,K-1)*DYF(K-1))*(u(I,J,K)-u(I,J,K-1))/(2.0_DP*(DY(K-1)**2)) &
                      +(v3(I,J,K)*DYF(K)+v3(I,J,K-1)*DYF(K-1))*(w(I,J,K)-w(I,J,K-1))/(2.0_DP*(DY(K-1)**2)) &
                      + 2.0_DP*(Q(I,J,K)-Q(I,J,K-1))/DY(K-1)

exp_spec_terms_w2_x(I,J,K) =exp_spec_terms_w2_x(I,J,K)+1.0_DP/Re*mutot_dbl_breve_rkp1(I,J,K)*(v1(I,J,K+1)-v1(I,J,K-1))/(2.0_DP*DY(K-1))
exp_spec_terms_w2_z(I,J,K) =exp_spec_terms_w2_z(I,J,K)+1.0_DP/Re*mutot_dbl_breve_rkp1(I,J,K)*(v3(I,J,K+1)-v3(I,J,K-1))/(2.0_DP*DY(K-1))

END FORALL
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, exp_spec_terms_w2_x, F_exp_spec_terms_w2_x )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, exp_spec_terms_w2_z, F_exp_spec_terms_w2_z )

FORALL (I=1:NX/2+1,J=1:NZ,K=0:NY+1)
F_Cranck_Exp_v2_temp(I,J,K)=ii*(kx(I)*F_exp_spec_terms_w2_x(I,J,K)+kz(J)*F_exp_spec_terms_w2_z(I,J,K))
END FORALL
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Cranck_Exp_v2_temp, Cranck_Exp_v2_temp )


Cranck_Exp_v2=Cranck_Exp_v2+Cranck_Exp_v2_temp

!! ------------------- Finishing off the adjoint v equation ------------------
! IMPLICIT PART IS NOW ADDED
A=0.0_DP
B=0.0_DP
C=0.0_DP
DO J = 1, NZ
  DO I = 1, NX
    DO K = K_Start, K_End
A(I,J,K) =  (V_bar_rkp1(I,J,K)+v_rkp1(I,J,K)+V_bar_rkp1(I,J,K-1)+v_rkp1(I,J,K-1))/(4.0_DP*DY(K-1))&
            - 2.0_DP/(Re*DYF(K-1)*DY(K-1))*mutot_rkp1(I,J,K-1)


C(I,J,K) = -(V_bar_rkp1(I,J,K+1)+v_rkp1(I,J,K+1)+V_bar_rkp1(I,J,K)+v_rkp1(I,J,K))/(4.0_DP*DY(K-1))&
           - 2.0_DP/(Re*DYF(K)*DY(K-1))*mutot_rkp1(I,J,K)

B(I,J,K) =1.0_DP - &
(V_bar_rkp1(I,J,K+1)+v_rkp1(I,J,K+1)-1.0_DP*(V_bar_rkp1(I,J,K-1)+v_rkp1(I,J,K-1)) &
)/(4.0_DP*DY(K-1))&
-(v(I,J,K+1) - v(I,J,K-1))/(2.0_DP*DY(K-1))&
+2.0_DP/(Re*DYF(K-1)*DY(K-1))*mutot_rkp1(I,J,K-1)+2.0_DP/(Re*DYF(K)*DY(K-1))*mutot_rkp1(I,J,K)

END DO
END DO
END DO

v2 = v2+ delta_t*(gamma(RK_step)*Exp_v2 + zeta(RK_step)*Exp_v2_m1 + &
                    (alpha(RK_step)/2.0_DP)*Cranck_Exp_v2)
Exp_v2_m1 = Exp_v2

CALL Adjoint_v_Boundary_Conditions(A,B,C,v2,NX,NY,NZ,v2_BC_Lower,v2_BC_Upper,v2_wall_Lower,v2_wall_Upper)
CALL Thomas_Matrix_Algorithm_real(A,B,C,v2,NX,NY,NZ)

alfa_t=-1.0_DP*alpha(RK_step)*delta_t
CALL Remove_Divergence(NX, NY, NZ, Lx, Lz, alfa_t, kx, kz, DY, DYF, plan_fwd, &
                       plan_bkd, v1, v2, v3, Q)


FORALL
Cranck_Exp_stau(I,J,K)=Cranck_Exp_stau(I,J,K) + &
&
-2.0_DP/Re(&
(mu_T_rkp1(I,J,K)*(ux_rkp1(I,J,K)+U_barx_rkp1(I,J,K))+mu_bar_T_rkp1(I,J,K)*ux_rkp1(I,J,K))*v1x(I,J,K) + &
&
(mu_T_rkp1(I,J,K)*((u_dbl_breve_rkp1(I,J,K+1)+U_bar_dbl_breve_rkp1(I,J,K+1)&
-u_dbl_breve_rkp1(I,J,K)-U_bar_dbl_breve_rkp1(I,J,K))/(2.0_DP*DYF(K))&
+0.25_DP*(vx_rkp1(I,J,K+1)+vx_rkp1(I,J,K)+V_barx_rkp1(I,J,K+1)+V_barx_rkp1(I,J,K)))+&
mu_bar_T_rkp1(I,J,K)*((u_dbl_breve_rkp1(I,J,K+1)-u_dbl_breve_rkp1(I,J,K))/(2.0_DP*DYF(K))&
+0.25_DP*(vx_rkp1(I,J,K+1)+vx_rkp1(I,J,K))))*(0.5_DP*(v2x(I,J,K+1)+v2x(I,J,K))&
+(v1_dbl_breve(I,J,K+1)-v1_dbl_breve(I,J,K))/DYF(K) )  +&
&
((mu_T_rkp1(I,J,K)*(uz_rkp1(I,J,K)+U_barz_rkp1(I,J,K)+wx_rkp1(I,J,K)+W_barx_rkp1(I,J,K))/2.0_DP+&
  mu_bar_T_rkp1(I,J,K)*(uz_rkp1(I,J,K)+wx_rkp1(I,J,K))/2.0_DP)*(v3x(I,J,K)+v1z(I,J,K)))+&
&
(mu_T_rkp1(I,J,K)*(wz_rkp1(I,J,K)+W_barz_rkp1(I,J,K))+mu_bar_T_rkp1(I,J,K)*wz_rkp1(I,J,K))*v3z(I,J,K) + &
&
(mu_T_rkp1(I,J,K)*((w_dbl_breve_rkp1(I,J,K+1)+W_bar_dbl_breve_rkp1(I,J,K+1)&
-w_dbl_breve_rkp1(I,J,K)-W_bar_dbl_breve_rkp1(I,J,K))/(2.0_DP*DYF(K))&
+0.25_DP*(vz_rkp1(I,J,K+1)+vz_rkp1(I,J,K)+V_barz_rkp1(I,J,K+1)+V_barz_rkp1(I,J,K)))+&
mu_bar_T_rkp1(I,J,K)*((w_dbl_breve_rkp1(I,J,K+1)-w_dbl_breve_rkp1(I,J,K))/(2.0_DP*DYF(K))&
+0.25_DP*(vz_rkp1(I,J,K+1)+vz_rkp1(I,J,K))))*(0.5_DP*(v2z_rkp1(I,J,K+1)+v2z_rkp1(I,J,K))&
+(v3_dbl_breve(I,J,K+1)-v3_dbl_breve(I,J,K))/DYF(K) )  +&
&
(mu_T_rkp1(I,J,K)*(v_rkp1(I,J,K+1)+V_bar_rkp1(I,J,K+1)-v_rkp1(I,J,K)-V_bar_rkp1(I,J,K))/DYF(K)+&
mu_bar_T_rkp1(I,J,K)*(v_rkp1(I,J,K+1)-v_rkp1(I,J,K))/DYF(K))*((v2(I,J,K+1)-v2(I,J,K))/DYF(K))&
)&
&
Ri*(n1*v1(I,J,K)+n2*v2(I,J,K)+n3*v3(I,J,K))

END FORALL

!! ------------------- Finishing off the adjoint stau equation ------------------
! IMPLICIT PART IS NOW ADDED
A=0.0_DP
B=0.0_DP
C=0.0_DP
DO J = 1, NZ
  DO I = 1, NX
    DO K = K_Start, K_End
A(I,J,K) =  (V_bar_rkp1(I,J,K)+v_rkp1(I,J,K))/(2.0_DP*DY(K-1))&
           -1.0_DP/(Re*Pr*DYF(K)*DY(K-1))
C(I,J,K) =  -(V_bar_rkp1(I,J,K+1)+v_rkp1(I,J,K+1))/(2.0_DP*DY(K))&
           -1.0_DP/(Re*Pr*DYF(K)*DY(K))
B(I,J,K) =1.0_DP - (V_bar_rkp1(I,J,K+1)+v_rkp1(I,J,K+1))*DYF(K+1)/(2.0_DP*DYF(K)*DY(K))&
          + (V_bar_rkp1(I,J,K)+v_rkp1(I,J,K))*DYF(K-1)/(2.0_DP*DYF(K)*DY(K-1)) &
          + 1.0_DP/(Re*PR*DYF(K))*(1.0_DP/DY(K)+1.0_DP/DY(K-1))
END DO
END DO
END DO

stau = stau+ delta_t*(gamma(RK_step)*Exp_stau + zeta(RK_step)*Exp_stau_m1 + &
                    (alpha(RK_step)/2.0_DP)*Cranck_Exp_stau)
Exp_stau_m1 = Exp_stau

CALL Adjoint_scalar_Boundary_Conditions(A,B,C,stau,NX,NY,NZ,stau_BC_Lower,stau_BC_Upper,stau_wall_Lower,stau_wall_Upper)
CALL Thomas_Matrix_Algorithm_real(A,B,C,stau,NX,NY,NZ)

END DO

END SUBROUTINE RK_Solver_Back
END MODULE Adjoint_Solvers
