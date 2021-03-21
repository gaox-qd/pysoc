      REAL FUNCTION CLEBSCH(J1,J2,M1,M2,J,M)
C
C     Purpose: Calculate Clebsch-Gordan coefficients using the 
C              Racah formula.
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'fileio.h'
      INCLUDE 'math.h'
C
      INTEGER :: J1,J2,M1,M2,J,M
      INTEGER :: T1,T2,T3,T4,T5,KMIN,KMAX,PHASE
      INTEGER :: K
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      CLEBSCH = 0.0 
C
      IF ((ABS(M1) > J1).OR.(ABS(M2) > J2)
     $                  .OR.(ABS(M) > J)
     $                  .OR.(J1 < 0)
     $                  .OR.(J2 < 0)
     $                  .OR.(J < 0)
     $                  .OR.(ABS(J1-J2) > J)
     $                  .OR.(J > J1+J2)
     $                  .OR.(M1+M2 /= M)) RETURN
C
      IF ((MAX(J1,J2)+J) > MAXFAC) THEN
        WRITE(OUT,1000) 
        WRITE(OUT,2000) MAXFAC,MAX(J1,J2)+J
      ENDIF
C
      T1 = (J1 - M1)/2
      T2 = (J - J2 + M1)/2
      T3 = (J2 + M2)/2
      T4 = (J - J1 - M2)/2
      T5 = (J1 + J2 - J)/2
C
      IF ((T1*2 /= J1-M1).OR.(T3*2 /= J2+M2)
     $                   .OR.(T5*2 /= J1+J2-J)) RETURN
C
      KMIN = MAX(MAX(-T2,-T4),0)
      KMAX = MIN(MIN(T1,T3),T5)
C
      PHASE = 1
      IF ((KMIN/2)*2 /= KMIN) PHASE = -1
C
      DO K=KMIN,KMAX
        CLEBSCH = CLEBSCH + PHASE/(FAC(T1-K)*FAC(T2+K)*
     $            FAC(T3-K)*FAC(T4+K)*FAC(K)*FAC(T5-K))
        PHASE = -PHASE
      END DO
C
      IF (KMIN > KMAX) CLEBSCH = 1.0 
      CLEBSCH = CLEBSCH*SQRT(FAC(T5)*FAC((J1+J-J2)/2)*
     $          FAC((J2+J-J1)/2)/FAC((J1+J2+J)/2+1)*(J+1)*
     $          FAC((J1+M1)/2)*FAC(T1)*FAC(T3)*FAC((J2-M2)/2)*
     $          FAC((J+M)/2)*FAC((J-M)/2))
C
C     
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
 1000 FORMAT(/,T2,'*** CLEBSCH: INCREASE THE NUMBER OF FACTORIALS ***')
 2000 FORMAT(/,T5,'ARE NEEDED ',I5,' INSTEAD OF ',I5) 
C
C     O________________________________________________________________O
C
      END FUNCTION CLEBSCH
