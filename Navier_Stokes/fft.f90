MODULE Fourier_Spectral
! ------------------------------------ Module Collecting all the subroutines for
! spectral differentiation in  x and z direction -------------------------------
! ---------------- It consists of the following subroutines: -------------------
! 1. FFT_Initialise: To initialise the plans for forward and backward FFTs
! 2. create_wavenumbers: create wavenumbers in x and z directions
! 3. create_wavenumbers: create 1 D arrays of wavenumbers in x and z directions
! 4. physical_to_fourier_2D: simply taking 2D fourier transform
! 5. fourier_to_physical_2D: simply taking inverse 2D fourier transform
! 6. dfdx: 1st derivative
! 7. df2dx2: 2nd z derivative
! 8. df2dxy: xy derivative
! 9. FFT_destroy: To destroy the plans
CONTAINS
! -------------- 1. Subroutine for initialising the FFTW routines --------------
! This is used to initialise the plans for forward and backward FFTW in x and z
! directions. It uses FFTW_MEASURE such that most optimal FFTW algorithms are
! chosen.
! Input: NX and NZ
! Output: plan_fwd and plan_bkd
SUBROUTINE FFT_Initialise(NX, NZ, plan_fwd, plan_bkd)
USE, INTRINSIC :: iso_c_binding
IMPLICIT NONE
include 'fftw3.f03'
INTEGER, INTENT(IN) :: NX, NZ
type(C_PTR), INTENT(OUT) :: plan_fwd, plan_bkd
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
REAL(KIND=DP), DIMENSION(NX,NZ) :: in_x
COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(NX/2+1,NZ) :: out_x
plan_fwd = fftw_plan_dft_r2c_2d (NZ, NX, in_x, out_x, FFTW_MEASURE)
plan_bkd = fftw_plan_dft_c2r_2d (NZ, NX, out_x, in_x, FFTW_MEASURE)
END SUBROUTINE FFT_Initialise

! -------------- 2. Subroutine for creating wavenumbers for FFTW ---------------
! -------------------------- AS OF NOW IT IS NOT USED --------------------------
! This is used to create the x and z create_wavenumbers
! Input: NX, NZ, Lx and Lz
! Output: kx and kz
SUBROUTINE create_wavenumbers ( NX, NZ, Lx, Lz, kx, kz)
IMPLICIT NONE
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER :: I, J
REAL(KIND=DP), PARAMETER :: pi=4.0_DP*ATAN(1.0_DP)
INTEGER, INTENT(IN) :: NX, NZ
REAL(KIND=DP), INTENT(IN) :: Lx, Lz
REAL(KIND=DP), DIMENSION(NX/2+1,NZ), INTENT(OUT) :: kx, kz
FORALL (I  = 1:NX/2, J = 1:NZ)
  kx(I,J) = 2.0_DP*pi/Lx*real(I-1)
END FORALL
kx(NX/2+1,:) = 0.0_DP
FORALL (I  = 1:NX/2, J = 1:NZ/2+1)
  kz(I,J) = 2.0_DP*pi/Lz*real(J-1)
END FORALL
FORALL (I  = 1:NX/2, J = NZ/2+2:NZ)
  kz(I,J) = -1*kz(I,NZ-J+2)
END FORALL
kz(NX/2+1,:) = kz(1,:)
END SUBROUTINE create_wavenumbers

! ------------- 3. Subroutine for creating 1D wavenumbers for FFTW -------------
! This is used to create the x and z create_wavenumbers
! Input: NX, NZ, Lx and Lz
! Output: kx and kz
SUBROUTINE create_wavenumbers_1D ( NX, NZ, Lx, Lz, kx, kz)
IMPLICIT NONE
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER :: I, J
REAL(KIND=DP), PARAMETER :: pi=4.0_DP*ATAN(1.0_DP)
INTEGER, INTENT(IN) :: NX, NZ
REAL(KIND=DP), INTENT(IN) :: Lx, Lz
REAL(KIND=DP), DIMENSION(1:NX/2+1), INTENT(OUT) :: kx
REAL(KIND=DP), DIMENSION(1:NZ), INTENT(OUT) :: kz
FORALL (I  = 1:NX/2)
  kx(I) = 2.0_DP*pi/Lx*real(I-1)
END FORALL
kx(NX/2+1) = 0.0_DP
DO J = 1,NZ/2+1
  kz(J) = 2.0_DP*pi/Lz*real(J-1)
END DO
DO J = NZ/2+2,NZ
  kz(J) = -1.0_DP*kz(NZ-J+2)
END DO
END SUBROUTINE create_wavenumbers_1D


! -------- 4. Subroutine for performing fast fourier transform in x and y:------
! --------------------------------- Cu(kx,kz,y)=fft(u(x,z,y)) ------------------
! Input: plan_fwd, fft_u, NX, NY, NZ
! Output: u
SUBROUTINE physical_to_fourier_2D( plan_fwd, NX, NY, NZ, u, fft_u )
USE, INTRINSIC :: iso_c_binding
IMPLICIT NONE
include 'fftw3.f03'
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: NX, NY, NZ
INTEGER :: K
type(C_PTR), INTENT(IN) :: plan_fwd
REAL(KIND=DP), DIMENSION(1:NX,1:NZ,0:NY+1), INTENT(INOUT) :: u
COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(1:NX/2+1,1:NZ,0:NY+1), INTENT(OUT) :: fft_u
DO K=0,Ny+1
  CALL fftw_execute_dft_r2c (plan_fwd, u(:,:,K), fft_u(:,:,K))
END DO
END SUBROUTINE physical_to_fourier_2D

! --- 5. Subroutine for performing inverse fast fourier transform in x and y:---
! --------------------------------- u(x,z,y)=ifft(Cu(kx,kz,y)) ----------------------------------
! Input: plan_bkd, fft_u, NX, NY, NZ
! Output: u
SUBROUTINE fourier_to_physical_2D( plan_bkd, NX, NY, NZ, fft_u, u )
USE, INTRINSIC :: iso_c_binding
IMPLICIT NONE
include 'fftw3.f03'
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
INTEGER, INTENT(IN) :: NX, NY, NZ
INTEGER :: K
type(C_PTR), INTENT(IN) :: plan_bkd
REAL(KIND=DP), DIMENSION(1:NX,1:NZ,0:NY+1), INTENT(OUT) :: u
COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(1:NX/2+1,1:NZ,0:NY+1), INTENT(INOUT) :: fft_u
DO K=0,Ny+1
  CALL fftw_execute_dft_c2r (plan_bkd, fft_u(:,:,K), u(:,:,K))
  u(:,:,K) = u(:,:,K)/size(u(:,:,K))
END DO
END SUBROUTINE fourier_to_physical_2D

! ----------- 6. Subroutine for calculating x (or `y`) derivative spectrally -----
! -------------------------- AS OF NOW IT IS NOT USED --------------------------
! This calculates the x derivative of the function f using the formula:
! ---------------------- du/dx=ifft(sqrt(-1)*k*fft(u)) -----------------------
! Input: plan_fwd, plan_bkd, NX, NZ, k, u
! Output: u_x
SUBROUTINE dfdx( plan_fwd, plan_bkd, NX, NZ, k, u, u_x)
USE, INTRINSIC :: iso_c_binding
IMPLICIT NONE
include 'fftw3.f03'
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
COMPLEX(C_DOUBLE_COMPLEX),PARAMETER :: ii=(0.d0, 1.d0)
type(C_PTR), INTENT(IN) :: plan_fwd, plan_bkd
INTEGER, INTENT(IN) :: NX, NZ
REAL(KIND=DP), DIMENSION(NX,NZ), INTENT(INOUT) :: u
REAL(KIND=DP), DIMENSION(NX,NZ), INTENT(OUT) :: u_x
COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(NX/2+1,NZ) :: fft_u
REAL(KIND=DP), DIMENSION(NX/2+1,NZ), INTENT(IN) :: k
CALL fftw_execute_dft_r2c (plan_fwd, u, fft_u)
fft_u = ii * k * fft_u
CALL fftw_execute_dft_c2r(plan_bkd, fft_u, u_x)
u_x = u_x/size(u_x)
END SUBROUTINE dfdx


! ----- 7. Subroutine for calculating 2nd x (or z) derivative spectrally -------
! -------------------------- AS OF NOW IT IS NOT USED --------------------------
! This caluclates the 2nd x derivative of the function f using the formula:
! ----------------------- d^2u/dx^2=ifft(-1*k^2*fft(u)) -----------------------
! Input: plan_fwd, plan_bkd, NX, NZ, k, u
! Output: u_xx
SUBROUTINE df2dx2( plan_fwd, plan_bkd, NX, NZ, k, u, u_xx)
USE, INTRINSIC :: iso_c_binding
IMPLICIT NONE
include 'fftw3.f03'
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
type(C_PTR), INTENT(IN) :: plan_fwd, plan_bkd
INTEGER, INTENT(IN) :: NX, NZ
REAL(KIND=DP), DIMENSION(NX,NZ), INTENT(INOUT) :: u
REAL(KIND=DP), DIMENSION(NX,NZ), INTENT(OUT) :: u_xx
COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(NX/2+1,NZ) :: fft_u
REAL(KIND=DP), DIMENSION(NX/2+1,NZ), INTENT(IN) :: k
CALL fftw_execute_dft_r2c (plan_fwd, u, fft_u)
fft_u = -1*fft_u*k**2
CALL fftw_execute_dft_c2r(plan_bkd, fft_u, u_xx)
u_xx = u_xx/size(u_xx)
END SUBROUTINE df2dx2


! ----- 8. Subroutine for calculating cross derivative spectrally -------
! -------------------------- AS OF NOW IT IS NOT USED --------------------------
! This caluclates the cross x- derivative of the function f using the formula:
! ----------------------- d^2u/dxy=ifft(-1*k_x,k_z*fft(u)) -----------------------
! Input: plan_fwd, plan_bkd, NX, NZ, k_x, k_z, u
! Output: u_xz
SUBROUTINE df2dxy( plan_fwd, plan_bkd, NX, NZ, k_x, k_z, u, u_xz)
USE, INTRINSIC :: iso_c_binding
IMPLICIT NONE
include 'fftw3.f03'
INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
type(C_PTR), INTENT(IN) :: plan_fwd, plan_bkd
INTEGER, INTENT(IN) :: NX, NZ
REAL(KIND=DP), DIMENSION(NX,NZ), INTENT(INOUT) :: u
REAL(KIND=DP), DIMENSION(NX,NZ), INTENT(OUT) :: u_xz
COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(NX/2+1,NZ) :: fft_u
REAL(KIND=DP), DIMENSION(NX/2+1,NZ), INTENT(IN) :: k_x, k_z
CALL fftw_execute_dft_r2c (plan_fwd, u, fft_u)
fft_u = -1*fft_u*k_x*k_z
CALL fftw_execute_dft_c2r(plan_bkd, fft_u, u_xz)
u_xz = u_xz/size(u_xz)
END SUBROUTINE df2dxy

! --------------- 9. Subroutine for destroying FFTW plans ----------------------
SUBROUTINE FFT_destroy(plan_fwd, plan_bkd)
USE, INTRINSIC :: iso_c_binding
IMPLICIT NONE
include 'fftw3.f03'
type(C_PTR), INTENT(IN) :: plan_fwd, plan_bkd
call fftw_destroy_plan(plan_fwd)
call fftw_destroy_plan(plan_bkd)
END SUBROUTINE FFT_destroy

END MODULE Fourier_Spectral
