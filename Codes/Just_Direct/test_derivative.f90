PROGRAM TEST
USE Grid_Definition
IMPLICIT NONE
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14), NX=20, NY=32, NZ=32
REAL(KIND=DP), PARAMETER :: pi = 4.0_DP*ATAN(1.0_DP), &
 	U_bulk=1.0_DP, &
  Kick_ini_vel=0.0_DP, Lx=2.0_DP*pi, Ly=2.0_DP, Lz=2.0_DP*pi, Stretch_y=0.75_DP,&
  n1=0.0_DP, n2=1.0_DP, n3=0.0_DP, Kick_ini_temp_fluct=0.0_DP

REAL(KIND=DP), DIMENSION(1:NX) :: GX
REAL(KIND=DP), DIMENSION(0:NY+1) :: GY, GYF
REAL(KIND=DP), DIMENSION(1:NZ) :: GZ, kz
REAL(KIND=DP), DIMENSION(0:NY) :: DY, DYF
REAL(KIND=DP), DIMENSION(0:NY+1) :: U,DU2_DY2,DY2_theory, errr
INTEGER :: K
call xyz_grid(Lx, Ly, Lz, Stretch_y, NX, NY, NZ, GX, GY, GYF, GZ, DY, DYF)

    DO K = 0, NY+1
U(K)=(1.0_DP-(GYF(K))**3)
DY2_theory(K)=-6.0_DP*(GYF(K)**1)
END DO
DY2_theory(0)=0.0_DP
DY2_theory(1)=0.0_DP
DY2_theory(NY)=0.0_DP
DY2_theory(NY+1)=0.0_DP

DU2_DY2=0.0_DP

DO K = 2, NY-1
DU2_DY2(K)= ((U(K+1)-U(K))/DY(K)-(U(K)-U(K-1))/DY(K-1))/((DY(K)+DY(K-1))/2.0_DP)
errr(K)=-2.0_DP*(DY(K)-DY(K-1))
END DO 

print *, maxval(abs(DU2_DY2-DY2_theory))
print *, maxval(DYF)/maxval(abs(DU2_DY2-DY2_theory))
print *, maxval(DYF)**2
print *, maxval(abs(errr))/maxval(abs(DU2_DY2-DY2_theory))
!print *, du2_dy2
END PROGRAM TEST
