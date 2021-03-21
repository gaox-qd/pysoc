      SUBROUTINE CRFTABLE
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'math.h'
C
      INTEGER, PARAMETER :: ITMAX = 50, NMAX = MXLT+6
      INTEGER :: J,K,N,ITAB
      REAL :: FACTOR,PROD,SUMTERM,SUMT,TERM,EPS,T
      REAL, DIMENSION(ITMAX + 10) :: R
C
C     ------------------------------------------------------------------
C
      EPS = EPSILON(0.0)
C      
      DO ITAB=1,MAXITAB
        T = FLOAT(ITAB)*0.1
        FTABLE(ITAB,0:NMAX) = 0.0
C     
        R(ITMAX+10) = 0.0
        DO J=ITMAX+9,1,-1
          R(J) = -T/((4.0*J + 2.0) - T*R(J+1))
        END DO
C
        FACTOR = 2.0*SINH(0.5*T)*EXP(-0.5*T)/T
C
        DO N=0,NMAX
C
C     ___ Initialize the iteration ___
C      
          SUMT = FACTOR/(2.0*N + 1.0)
          TERM = 1.0
C
C     ___ Start the summation and recursion ___
C      
          DO K=1,ITMAX
            TERM = TERM*(2.0*N - 2.0*K + 1.0)/(2.0*N + 2.0*K + 1.0)
C
C     ___ Product of Bessel function quotients ___
C
            PROD = 1.0
            DO J=1,K
              PROD = PROD*R(J)
            END DO
C
            SUMTERM = FACTOR*TERM*PROD*(2.0*K + 1.0)/(2.0*N + 1.0)
C
C     ___ Convergence test ___
C
            IF (ABS(SUMTERM) < EPS) THEN
C
              EXIT
C
            ELSE IF (K == ITMAX) THEN
C
              STOP 'MAXIMUM NUMBER OF ITERATIONS REACHED IN CRFTABLE'
C
            ELSE
C      
              SUMT = SUMT + SUMTERM
C
            END IF
          END DO
C
          FTABLE(ITAB,N) = SUMT
C
        END DO
      END DO
C
C     ------------------------------------------------------------------
C
      END SUBROUTINE CRFTABLE
