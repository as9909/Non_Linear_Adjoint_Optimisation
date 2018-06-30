MODULE Adjoint_Solvers

! Make current step rk and the previous rkm1 instead of how it is presently: where current step is rkp1 and previous is rk

SUBROUTINE RK_Solver_Back



stau_dblebreve=0.0_DP
FORALL  (I=1:NX, J=1:NZ, K=K_start:K_end)
stau_dblebreve(I,J,K)=(stau(I,J,K)*DYF(K-1)+stau(I,J,K-1)*DYF(K))/(2*DY(K-1))
END FORALL
INTENT(IN):: TH_p_THbar, TH_bar
TH=TH_p_THbar-TH_bar ! Remove laminar temperature profile from the total real temperature
TH_rkp1=TH_p_THbar_rkp1-TH_bar_rkp1
THpTHB=TH_p_THbar+THB

CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, u, F_u )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, v, F_v )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, v_rkp1 ,F_v_rkp1 )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, w, F_w )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, TH, F_TH )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, TH_rkp1, F_TH_rkp1 )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, THBpTHB, F_THBpTHB )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, v1, F_v1 )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, v2, F_v2 )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, v3, F_v3 )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, Q, F_Q )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, tau, F_tau )

CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, U_bar, F_U_bar )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, V_bar, F_V_bar )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, W_bar, F_W_bar )

FORALL (I=1:NX/2+1,J=1:NZ,K=0:NY+1)
F_ux(I,J,K)=ii*kx(I)*F_u(I,J,K)
F_vx(I,J,K)=ii*kx(I)*F_v(I,J,K)
F_v_rkp1x(I,J,K) =ii*kx(I)*F_v_rkp1(I,J,K)
F_wx(I,J,K)=ii*kx(I)*F_w(I,J,K)
F_uz(I,J,K)=ii*kz(J)*F_u(I,J,K)
F_vz(I,J,K)=ii*kz(J)*F_v(I,J,K)
F_v_rkp1z(I,J,K) =ii*kz(J)*F_v_rkp1(I,J,K)
F_wz(I,J,K)=ii*kz(J)*F_w(I,J,K)

F_THBpTHBx(I,J,K)=ii*kx(I)*F_THpTHB(I,J,K)
F_THBpTHBz(I,J,K)=ii*kz(J)*F_THpTHB(I,J,K)

F_v1x(I,J,K)=ii*kx(I)*F_v1(I,J,K)
F_v2x(I,J,K)=ii*kx(I)*F_v2(I,J,K)
F_v3x(I,J,K)=ii*kx(I)*F_v3(I,J,K)
F_Qx(I,J,K) =ii*kx(I)*F_Q(I,J,K)
F_v1z(I,J,K)=ii*kz(J)*F_v1(I,J,K)
F_v2z(I,J,K)=ii*kz(J)*F_v2(I,J,K)
F_v3z(I,J,K)=ii*kz(J)*F_v3(I,J,K)
F_Qz(I,J,K) =ii*kz(J)*F_Q(I,J,K)

F_U_barx(I,J,K)=ii*kx(I)*F_U_bar(I,J,K)
F_U_barz(I,J,K)=ii*kz(J)*F_U_bar(I,J,K)
F_V_barx(I,J,K)=ii*kx(I)*F_V_bar(I,J,K)
F_V_barz(I,J,K)=ii*kz(J)*F_V_bar(I,J,K)
F_W_barx(I,J,K)=ii*kx(I)*F_W_bar(I,J,K)
F_W_barz(I,J,K)=ii*kz(J)*F_W_bar(I,J,K)

END FORALL

CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_ux, ux )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_vx, vx )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_v_rkp1x, v_rkp1x )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_wx, wx )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_uz, uz )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_vz, vz )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_v_rkp1z, v_rkp1z )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_wz, wz )

CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_U_barx, U_barx )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_U_barz, U_barz )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_V_barx, V_barx )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_V_barz, V_barz )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_W_barx, W_barx )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_W_barz, W_barz )

CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_THBpTHBx, THBpTHBx )
CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_THBpTHBz, THBpTHBz )

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

tau_Utot=tau*(u+U_bar)
tau_Wtot=tau*(w+W_bar)
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, tau_Utot, F_tau_Utot )
CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, tau_Wtot, F_tau_Wtot )

FORALL (I=1:NX/2+1,J=1:NZ,K=0:NY+1)
F_Exp_v1(I,J,K)=ii*(kx(I)*F_v1Utot_p_mutot_v1x(I,J,K)+kz(J)*F_v1Wtot_p_mutot_v1z_v3x(I,J,K))
F_Exp_v3(I,J,K)=ii*(kx(I)*F_v3Utot_p_mutot_v3x(I,J,K)+kz(J)*F_v3Wtot_p_mutot_v3z(I,J,K))
F_Exp_v2(I,J,K)=ii*(kx(I)*F_v2Utot_p_mutot_v2x(I,J,K)+kz(J)*F_v2Wtot_p_mutot_v2z(I,J,K))
F_Exp_stau(I,J,K)=ii*(kx(I)*F_tau_Utot(I,J,K)+kz(J)*F_tau_Wtot(I,J,K))-1.0_DP/(Re*Pr)*(kx(I)**2+kz(J)**2)*F_tau(I,J,K)
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
Exp_v2(I,J,K)=Exp_v2(I,J,K)-stau_dblebreve(I,J,K)*(THpTHB(I,J,K)-THpTHB(I,J,K-1))/(DY(K-1))
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
   -1.0_DP*A2*(mutot_rkp1_dbl_breve(I,J,K+1)*(u_rkp1(I,J,K+1)-u_rkp1(I,J,K))/DY(K) &
   - mutot_rkp1_dbl_breve(I,J,K)*(u_rkp1(I,J,K)-u_rkp1(I,J,K-1))/DY(K-1)))



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
                            -A2/Re*(mutot_dbl_breve(I,J,K)*vx(I,J,K)+mutot_rkp1_dbl_breve(I,J,K)*v_rkp1x(I,J,K))
exp_spec_terms_w2_z(I,J,K) =1.0_DP/Re*mutot_dbl_breve(I,J,K)*(v3(I,J,K+1)-v3(I,J,K-1))/(2.0_DP*DY(K-1)) &
                            -A2/Re*(mutot_dbl_breve(I,J,K)*vz(I,J,K)+mutot_rkp1_dbl_breve(I,J,K)*v_rkp1z(I,J,K))


Cranck_Exp_stau(I,J,K)= Cranck_Exp_stau(I,J,K)&
&
+(stau_dblebreve(I,J,K+1)*(V_bar(I,J,K+1)+v(I,J,K+1))&
                  -stau_dblebreve(I,J,K)*(V_bar(I,J,K)+v(I,J,K)))/DYF(K) &
&
+ (1.0_DP/(Re*Pr))*((tau(I,J,K+1)-tau(I,J,K))/DY(K)&
                    -(tau(I,J,K)-tau(I,J,K-1))/DY(K-1))/DYF(K) &
&
-2.0_DP*A2*Ri/(Re*Pr*(T_ref**2))*&
((TH(I,J,K+1)-TH(I,J,K))/DY(K)-(TH(I,J,K)-TH(I,J,K-1))/DY(K-1)&
&
+(TH_rkp1(I,J,K+1)-TH_rkp1(I,J,K))/DY(K)-(TH_rkp1(I,J,K)-TH_rkp1(I,J,K-1))/DY(K-1))/DYF(K) &
&
+A2*(mu_T(I,J,K)+mu_bar_T(I,J,K))/Re*(ux(I,J,K)**2+uy(I,J,K)**2+uz(I,J,K)**2&
                +vx(I,J,K)**2+vy(I,J,K)**2+vz(I,J,K)**2)&
                +wx(I,J,K)**2+wy(I,J,K)**2+wz(I,J,K)**2)&
&
+A2*(mu_T_rkp1(I,J,K)+mu_bar_T_rkp1(I,J,K))/Re*(ux_rkp1(I,J,K)**2+uy_rkp1(I,J,K)**2+uz_rkp1(I,J,K)**2&
                +vx_rkp1(I,J,K)**2+vy_rkp1(I,J,K)**2+vz_rkp1(I,J,K)**2)&
                +wx_rkp1(I,J,K)**2+wy_rkp1(I,J,K)**2+wz_rkp1(I,J,K)**2)&
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
               - mutot_rkp1_dbl_breve(I,J,K)/(Re*DYF(K)*DY(K-1))
    C(I,J,K) = -(V_bar_rkp1(I,J,K+1)+v_rkp1(I,J,K+1))/(2.0_DP*DYF(K))&
               - mutot_rkp1_dbl_breve(I,J,K+1)/(Re*DYF(K)*DY(K))
    B(I,J,K) =1.0_DP - (V_bar_rkp1(I,J,K+1)+v_rkp1(I,J,K+1))/(2.0_DP*DYF(K))&
              + (V_bar_rkp1(I,J,K)+v_rkp1(I,J,K))/(2.0_DP*DYF(K)) &
              + mutot_rkp1_dbl_breve(I,J,K)/(Re*DYF(K)*DY(K-1)) &
              + mutot_rkp1_dbl_breve(I,J,K+1)/(Re*DYF(K)*DY(K))
END DO
END DO
END DO

v1=v1+ delta_t*(gamma(RK_step)*Exp_v1 + zeta(RK_step)*Exp_v1_m1 + &
                    (alpha(RK_step)/2.0_DP)*Cranck_Exp_v1)
Exp_v1_m1=Exp_v1

CALL Adjoint_uw_Boundary_Conditions(A,B,C,v1,NX,NY,NZ,v1_BC_Lower,v1_BC_Upper)
CALL Thomas_Matrix_Algorithm_real(A,B,C,v1,NX,NY,NZ)

CALL physical_to_fourier_2D( plan_fwd, NX, NY, NZ, v1, F_v1 )
F_v1x(I,J,K)=ii*kx(I)*F_v1(I,J,K)
F_v1z(I,J,K)=ii*kz(J)*F_v1(I,J,K)
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
 -1.0_DP*A2*(mutot_rkp1_dbl_breve(I,J,K+1)*(w_rkp1(I,J,K+1)-w_rkp1(I,J,K))/DY(K) &
 - mutot_rkp1_dbl_breve(I,J,K)*(w_rkp1(I,J,K)-w_rkp1(I,J,K-1))/DY(K-1)))
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
           - mutot_rkp1_dbl_breve(I,J,K)/(Re*DYF(K)*DY(K-1))
C(I,J,K) = -(V_bar_rkp1(I,J,K+1)+v_rkp1(I,J,K+1))/(2.0_DP*DYF(K))&
           - mutot_rkp1_dbl_breve(I,J,K+1)/(Re*DYF(K)*DY(K))
B(I,J,K) =1.0_DP - (V_bar_rkp1(I,J,K+1)+v_rkp1(I,J,K+1))/(2.0_DP*DYF(K))&
          + (V_bar_rkp1(I,J,K)+v_rkp1(I,J,K))/(2.0_DP*DYF(K)) &
          + mutot_rkp1_dbl_breve(I,J,K)/(Re*DYF(K)*DY(K-1)) &
          + mutot_rkp1_dbl_breve(I,J,K+1)/(Re*DYF(K)*DY(K))
END DO
END DO
END DO

v3=v3+ delta_t*(gamma(RK_step)*Exp_v3 + zeta(RK_step)*Exp_v3_m1 + &
                    (alpha(RK_step)/2.0_DP)*Cranck_Exp_v3)
Exp_v3_m1=Exp_v3

CALL Adjoint_uw_Boundary_Conditions(A,B,C,v3,NX,NY,NZ,v3_BC_Lower,v3_BC_Upper)
CALL Thomas_Matrix_Algorithm_real(A,B,C,v3,NX,NY,NZ)

! ------


CALL fourier_to_physical_2D( plan_bkd, NX, NY, NZ, F_Cranck_Exp_v2, Cranck_Exp_v2_temp )

FORALL (I=1:NX, J=1:NZ, K=1:NY+1)

Cranck_Exp_v2(I,J,K) = Cranck_Exp_v2(I,J,K) &
                      +(v1(I,J,K)*DYF(K)+v1(I,J,K-1)*DYF(K-1))*(u(I,J,K)-u(I,J,K-1))/(2.0_DP*(DY(K-1)**2)) &
                      +(v3(I,J,K)*DYF(K)+v3(I,J,K-1)*DYF(K-1))*(w(I,J,K)-w(I,J,K-1))/(2.0_DP*(DY(K-1)**2)) &
                      + 2.0_DP*(Q(I,J,K)-Q(I,J,K-1))/DY(K-1)

exp_spec_terms_w2_x(I,J,K) =exp_spec_terms_w2_x(I,J,K)+1.0_DP/Re*mutot_rkp1_dbl_breve(I,J,K)*(v1(I,J,K+1)-v1(I,J,K-1))/(2.0_DP*DY(K-1))
exp_spec_terms_w2_z(I,J,K) =exp_spec_terms_w2_z(I,J,K)+1.0_DP/Re*mutot_rkp1_dbl_breve(I,J,K)*(v3(I,J,K+1)-v3(I,J,K-1))/(2.0_DP*DY(K-1))

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

v2=v2+ delta_t*(gamma(RK_step)*Exp_v2 + zeta(RK_step)*Exp_v2_m1 + &
                    (alpha(RK_step)/2.0_DP)*Cranck_Exp_v2)
Exp_v2_m1=Exp_v2

CALL Adjoint_uw_Boundary_Conditions(A,B,C,v2,NX,NY,NZ,v3_BC_Lower,v3_BC_Upper)
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

stau=stau+ delta_t*(gamma(RK_step)*Exp_stau + zeta(RK_step)*Exp_stau_m1 + &
                    (alpha(RK_step)/2.0_DP)*Cranck_Exp_stau)
Exp_stau_m1=Exp_stau

CALL Adjoint_uw_Boundary_Conditions(A,B,C,stau,NX,NY,NZ,stau_BC_Lower,stau_BC_Upper)
CALL Thomas_Matrix_Algorithm_real(A,B,C,stau,NX,NY,NZ)


END SUBROUTINE RK_Solver_Back
END MODULE Adjoint_Solvers
