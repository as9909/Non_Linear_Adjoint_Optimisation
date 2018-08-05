PROGRAM ThomasTest

  INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14), Ny=10
INTEGER :: I
REAL(KIND=DP), DIMENSION(0:Ny+1)::A,B,C, A_old,B_old,C_old, Check,D,D_old
CALL RANDOM_SEED

DO I=0,Ny
  CALL RANDOM_NUMBER(Rand_num)
  A(I)=Rand_num
    CALL RANDOM_NUMBER(Rand_num)
    B(I)=Rand_num
      CALL RANDOM_NUMBER(Rand_num)
      C(I)=Rand_num
        CALL RANDOM_NUMBER(Rand_num)
        D(I)=Rand_num
END DO
A(0)=0
C(NY+1)=0
A_old=A
B_old=B
C_old=C
D_old=D
C(0) = C(0)/B(0)
D(0) = D(0)/B(0)
DO I = 1,NY
  C(I) =  C(I)/(B(I)-A(I)*C(I-1))
  D(I) =  (D(I)-A(I)*D(I-1))/(B(I)-A(I)*C(I-1))
END DO
 D(NY+1) =  (D(NY+1)-A(NY+1)*D(NY))/(B(NY+1)&
          -A(NY+1)*C(NY))
DO I = NY,0,-1
  D(I) = D(I)-C(I)*D(I+1)
END DO
DO I=0,NY+1
  Check(I)=A_old(I)*D(I-1)+B_old(I)*D(I)+C_old(I)*D(I+1)-D_Old(I)
END DO
print*, maxval(Check)
END PROGRAM ThomasTest
