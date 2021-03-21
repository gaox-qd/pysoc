      SUBROUTINE READINP
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'fileio.h'
      INCLUDE 'molecule.h'
      INCLUDE 'moleculeD.h'
C
      CHARACTER :: LINE*80,REVMA*1,REVMB*1
C
      INTEGER :: IATOM,I,ILINE
C
      LOGICAL :: SUCCESS
C
      CHARACTER :: UCC*7
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Read input line and keywords ___
C
      ILINE = 0
   50 CONTINUE
      ILINE = ILINE + 1
      READ(INP,"(A80)",ERR=9999) LINE
      IF (LINE(1:1) == '#') GOTO 50
      SUCCESS = .FALSE.
      CALL FINKEY(LINE,SUCCESS)
      IF (.NOT.SUCCESS) THEN 
        BACKSPACE(INP)
        ILINE = ILINE - 1
      END IF
C
C   70 CONTINUE
C      ILINE = ILINE + 1
C      READ(INP,"(A80)",ERR=9999) LINE
C      IF (LINE(1:1) == '#') GOTO 70
C
C     ___ Read molecule specification ___
C
C      BACKSPACE(INP)
C      READ(INP,*,ERR=9999) CHARGE,MULTF,REVMA,MULTS,REVMB
C
C      DISCIN = MULTS - MULTF
C      DIRECTA = .TRUE.
C      DIRECTB = .TRUE.
C
C      IF (DISCIN == -2) THEN
C        WRITE(OUT,500) 
C        STOP
C      END IF
C
C      IF (UCC(REVMA(1:1)) == 'D') THEN
C        DIRECTA = .TRUE.
C      ELSE IF (UCC(REVMA(1:1)) == 'R') THEN
C        DIRECTA = .FALSE.
C      ELSE
C        WRITE(OUT,1000) 
C        STOP
C      END IF
C
C      IF (UCC(REVMB(1:1)) == 'D') THEN
C        DIRECTB = .TRUE.
C      ELSE IF (UCC(REVMB(1:1)) == 'R') THEN
C        DIRECTB = .FALSE.
C      ELSE
C        WRITE(OUT,2000) 
C        STOP
C      END IF
C
C      CALL ERRORS
C
C     ___ Read geometry ___
C
  100 CONTINUE
      ILINE = ILINE + 1
      READ(INP,"(A80)",ERR=9999) LINE
      IF (LINE(1:1) == '#') GOTO 100
      SUCCESS = .FALSE.
      CALL FINSTR(LINE,SUCCESS,'END')
      IF (SUCCESS) GOTO 200
      BACKSPACE(INP)
C
      NAT = NAT + 1
      IF (NAT.GT.MXATM) THEN 
        WRITE(OUT,3000) 
        STOP
      END IF
      READ(INP,*,ERR=9999) ELEMENT(NAT),(C(I,NAT),I=1,3),SCAL(NAT)
      CALL ATOMDATA(NAT,ILINE)
      GOTO 100
C
  200 CONTINUE
      NATOM = NAT 
C
C     ___ Transform coordinates in bohr ___
C
      DO IATOM=1,NAT
        C(:,IATOM) = C(:,IATOM)*TRANS
      END DO
C
      RETURN
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
 9999 CONTINUE
      WRITE(OUT,"(/,'ERROR READING INPUT LINE ',I5)") ILINE
      STOP
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Format statements ___
C
  500 FORMAT(T2,'READINP: UNKNOWN DISCOINCIDENCE OPTION ',/,
     $          'COUPLING MUST BE BETWEEN LOW SPIN AND HIGH SPIN')
 1000 FORMAT(T2,'READINP: UNKNOWN COUPLING ORBITAL MODE OF THE FIRST ',
     $          'STATE')
 2000 FORMAT(T2,'READINP: UNKNOWN COUPLING ORBITAL MODE OF THE SECOND ',
     $          'STATE')
 3000 FORMAT(T2,'NUMBER OF ATOMS EXCEEDED.',/,
     $       T2,'INCREASE THE PARAMETER MXATM IN parameter.h AND',/,
     $       T2,'RECOMPILE AGAIN') 
C
C     O________________________________________________________________O
C
      END SUBROUTINE READINP
      
