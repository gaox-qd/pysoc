      SUBROUTINE CPUTIME 
C
C     *****************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'fileio.h'
      INCLUDE 'moleculeD.h'
C
      REAL(KIND = 4) :: ETIME
      REAL(KIND = 4), DIMENSION(2) :: TIMEVEC
      REAL :: TIMER
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      TIMEVEC(1) = 0.0; TIMEVEC(2) = 0.0
      TIMER = ETIME(TIMEVEC)
      CPUT = MAX(TIMER,TOLTIME)
C
C     O________________________________________________________________O
C
      END SUBROUTINE CPUTIME
