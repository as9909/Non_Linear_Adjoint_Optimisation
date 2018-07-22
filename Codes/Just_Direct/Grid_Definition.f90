MODULE Grid_Definition
! ----------------------- Module for generating grid ---------------------------
! ---------------- It consists of the following subroutines: -------------------
! 1. xyz_grid: To generate x, y and z grid arrays (1 D arrays)
! 2. full_xyz_grid: Make grid arrays 3 D

! -------------------------------- 1. xyz_grid ---------------------------------
! Input : Nx, Ny, Nz, Lx, Ly, Lz, stretch_y
! Output : GX(x base grid), GY(y base grid), GYF(y fractional grid),
!          GZ(z base grid), DY(Delta y base grid), DYF(Delta y fractional grid)
CONTAINS
SUBROUTINE xyz_grid(Lx, Ly, Lz, Stretch_y, NX, NY, NZ, GX, GY, GYF, GZ, DY, DYF)
IMPLICIT NONE
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: NX, NY, NZ
INTEGER :: I
REAL(KIND=DP), INTENT(IN) :: Lx, Ly, Lz, Stretch_y
REAL(KIND=DP), DIMENSION(1:NX), INTENT(OUT) :: GX
REAL(KIND=DP), DIMENSION(0:NY+1), INTENT(OUT) :: GY, GYF
REAL(KIND=DP), DIMENSION(0:NY), INTENT(OUT) :: DY, DYF
REAL(KIND=DP), DIMENSION(1:NZ), INTENT(OUT) :: GZ
REAL(KIND=DP) :: fact, shift, GY_lower, GY_upper, GYF_lower, GYF_upper
DO I=1,NX
  GX(I)=(I*Lx)/NX
END DO
DO I=1,NZ
  GZ(I)=(I*Lz)/NZ
END DO
DO I=1,NY+1
  GY(I)=Ly/2.0_DP*TANH(Stretch_y*((2.0_DP*(I-1))/NY-1.0_DP))/TANH(Stretch_y)
END DO
DO I=1,NY
  GYF(I)=(GY(I)+GY(I+1))/2.0_DP
END DO
GY_lower=GY(1)
GY_upper=GY(NY+1)
GYF_lower=GYF(1)
GYF_upper=GYF(NY)
fact= (GY_upper-GY_lower)/(GYF_upper-GYF_lower)
GY=GY*fact
GYF=GYF*fact
GYF_lower=GYF(1)
GYF_upper=GYF(NY)
shift=GYF_lower-GY_lower
GYF=GYF-shift
GY=GY-shift
!GYF(0)=2.0_DP*GYF(1)-GYF(2)
!GYF(NY+1)=2.0_DP*GYF(NY)-GYF(NY-1)
GY(0)=2.0_DP*GY(1)-GY(2)
GYF(0)=(GY(0)+GY(1))/2.0_DP
GYF(NY+1)=2.0_DP*GYF(NY)-GYF(NY-1)
DO I=0,NY
  DY(I)=GYF(I+1)-GYF(I)
  DYF(I)=GY(I+1)-GY(I)
END DO
END SUBROUTINE xyz_grid
! ------------------------------------------------------------------------------

! ----------------------------- 2. full_xyz_grid -------------------------------
! -------------------------- AS OF NOW IT IS NOT USED --------------------------
! Input : Nx, Ny, Nz, and 1 D grids- GX, GY, GYF and GZ
! Output : 3D grids- GX_base, GY_base, GZ_base and GY_frac
SUBROUTINE full_xyz_grid(NX, NY, NZ, GX,GY,GYF,GZ, GX_base, GY_base, GZ_base,&
                                                                      GY_frac)
IMPLICIT NONE
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: NX, NY, NZ
INTEGER :: I, J, K
REAL(KIND=DP), DIMENSION(1:NX), INTENT(IN) :: GX
REAL(KIND=DP), DIMENSION(0:NY+1), INTENT(IN) :: GY, GYF
REAL(KIND=DP), DIMENSION(1+NZ), INTENT(IN) :: GZ
REAL(KIND=DP), DIMENSION(1:NX,1:NZ,0:NY+1), INTENT(OUT) :: GX_base, GY_base,&
                                                           GZ_base, GY_frac
FORALL (I=1:NX,J=1:NZ,K=0:NY+1)
GX_base(I,J,K)=GX(I)
GZ_base(I,J,K)=GZ(J)
GY_base(I,J,K)=GY(K)
GY_frac(I,J,K)=GYF(K)
END FORALL
END SUBROUTINE full_xyz_grid
END MODULE Grid_Definition
