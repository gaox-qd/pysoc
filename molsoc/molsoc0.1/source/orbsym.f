      SUBROUTINE ORBSYM
C
C     *****************************************************************
C
      IMPLICIT NONE
      INCLUDE 'parameter.h'
      INCLUDE 'fileio.h'
      INCLUDE 'molecule.h'
      INCLUDE 'moleculeD.h'
C
      INTEGER :: IS,L
      INTEGER, DIMENSION(10) :: POINT
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Initialization ___
C
      POINT = 0
C
      DO IS=1,NSHL
        L = PSHELL(IS,2)
        IF (L == 0) THEN
          ORSYM(LLSTO(IS)) = 'S'
        ELSE IF (L == 1) THEN
          POINT(1) = LLSTO(IS)
          POINT(2) = LLSTO(IS) + 1
          POINT(3) = ULSTO(IS)
          ORSYM(POINT(1)) = 'PX'
          ORSYM(POINT(2)) = 'PY'
          ORSYM(POINT(3)) = 'PZ'
        ELSE IF (L == 2) THEN
          POINT(1) = LLSTO(IS)
          POINT(2) = LLSTO(IS) + 1
          POINT(3) = LLSTO(IS) + 2
          POINT(4) = LLSTO(IS) + 3
          POINT(5) = LLSTO(IS) + 4
          POINT(6) = ULSTO(IS)
          ORSYM(POINT(1)) = 'DXX'
          ORSYM(POINT(2)) = 'DXY'
          ORSYM(POINT(3)) = 'DXZ'
          ORSYM(POINT(4)) = 'DYY'
          ORSYM(POINT(5)) = 'DYZ'
          ORSYM(POINT(6)) = 'DZZ'
        ELSE IF (L == 3) THEN
          POINT(1) = LLSTO(IS)
          POINT(2) = LLSTO(IS) + 1
          POINT(3) = LLSTO(IS) + 2
          POINT(4) = LLSTO(IS) + 3
          POINT(5) = LLSTO(IS) + 4
          POINT(6) = LLSTO(IS) + 5
          POINT(7) = LLSTO(IS) + 6
          POINT(8) = LLSTO(IS) + 7
          POINT(9) = LLSTO(IS) + 8
          POINT(10) = ULSTO(IS)
          ORSYM(POINT(1))  = 'FXXX'
          ORSYM(POINT(2))  = 'FXXY'
          ORSYM(POINT(3))  = 'FXXZ'
          ORSYM(POINT(4))  = 'FXYY'
          ORSYM(POINT(5))  = 'FXYZ'
          ORSYM(POINT(6))  = 'FXZZ'
          ORSYM(POINT(7))  = 'FYYY'
          ORSYM(POINT(8))  = 'FYYZ'
          ORSYM(POINT(9))  = 'FYZZ'
          ORSYM(POINT(10)) = 'FZZZ'
        ELSE 
          WRITE(OUT,*) 'ORBSYM: L CANNOT BE > 3'
          STOP
        END IF
      END DO
C
      POINT = 0
C
      DO IS=1,NSHL
        L = PSHELL(IS,2)
        IF (L == 0) THEN
          ORSYMS(SLL(IS)) = 's'
        ELSE IF (L == 1) THEN
          POINT(1) = SLL(IS)
          POINT(2) = SLL(IS) + 1
          POINT(3) = SUL(IS)
          ORSYMS(POINT(1)) = 'py'
          ORSYMS(POINT(2)) = 'pz'
          ORSYMS(POINT(3)) = 'px'
        ELSE IF (L == 2) THEN
          POINT(1) = SLL(IS)
          POINT(2) = SLL(IS) + 1
          POINT(3) = SLL(IS) + 2
          POINT(4) = SLL(IS) + 3
          POINT(5) = SUL(IS)
          ORSYMS(POINT(1)) = 'd-2'
          ORSYMS(POINT(2)) = 'd-1'
          ORSYMS(POINT(3)) = 'd 0'
          ORSYMS(POINT(4)) = 'd+1'
          ORSYMS(POINT(5)) = 'd+2'
        ELSE IF (L == 3) THEN
          POINT(1) = SLL(IS)
          POINT(2) = SLL(IS) + 1
          POINT(3) = SLL(IS) + 2
          POINT(4) = SLL(IS) + 3
          POINT(5) = SLL(IS) + 4
          POINT(6) = SLL(IS) + 5
          POINT(7) = SUL(IS)
          ORSYMS(POINT(1)) = 'f-3'
          ORSYMS(POINT(2)) = 'f-2'
          ORSYMS(POINT(3)) = 'f-1'
          ORSYMS(POINT(4)) = 'f 0'
          ORSYMS(POINT(5)) = 'f+1'
          ORSYMS(POINT(6)) = 'f+2'
          ORSYMS(POINT(7)) = 'f+3'
        ELSE
          WRITE(OUT,*) 'ORBSYM: L CANNOT BE > 3'
        END IF
      END DO
C
C     O_______________________________________________________________O
C
      END SUBROUTINE ORBSYM
