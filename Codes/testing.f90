PROGRAM Test_Filename
IMPLICIT NONE

INTEGER :: I=5, J, i1, i2, i3, i4
REAL :: K=9.0, L
REAL(KIND=16) :: Rand_num
INTEGER, DIMENSION(:), ALLOCATABLE :: seed
REAL(KIND=16), DIMENSION(10,10,10,3):: Y, Z
CHARACTER(LEN=17) :: filename


write(filename,'(A,I3.3,A)') 'U_total_0_',I,'.bin'

CALL RANDOM_SEED

DO  i1=1,10
DO i2=1,10
DO i3=1,10
DO i4=1,3
CALL RANDOM_NUMBER(Rand_num)
Y(i1,i2,i3,i4)=Rand_num
END DO
END DO
END DO
END DO

open(unit=11,file=filename,form='unformatted',status='replace')
write(11) ((((Y(i1,i2,i3,i4), &
     i1=1,10),i2=1,10),i3=1,10),i4=1,3)
close(11)
! print *,K
L=0.0
open(unit=12,file=filename,form='unformatted',status='old')
read(12) ((((Z(i1,i2,i3,i4), &
     i1=1,10),i2=1,10),i3=1,10),i4=1,3)
close(12)

print *, Y(1,1,1,1), Z(1,1,1,1),Y(2,1,1,1), Z(2,1,1,1)
open(unit=13,file='tt.dat',status='replace')
write(13,*) ((((Z(i1,i2,i3,i4), &
     i1=1,10),i2=1,10),i3=1,10),i4=1,3)
close(13)


END PROGRAM Test_Filename
