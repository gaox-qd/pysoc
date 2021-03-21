      SUBROUTINE FINELE(LINE,SUCCESS,STRING,INAT,IATOM,CLINE) 
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
      CHARACTER*(*) :: LINE,EL*2
C
      INTEGER :: I,ILINE,ILOW,ISTART,INU,LIM,IATOM,RNUM,CNUM,CLINE
C
      INTEGER, DIMENSION(80) :: LENSTR,LLSTR,ULSTR
C
      INTEGER, DIMENSION(NATOM) :: INAT
C
      LOGICAL :: SUCCESS
C
      CHARACTER*(*) :: UCC*7,STRING,IFOR*4,CHECKL*25,CHECKS*4,LINE1*1
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Find keywords decoding the line ___
C
      ILINE = 1
      ILOW = 0
      LENSTR = 0
      LLSTR = 0
      INAT = 0
      SUCCESS = .FALSE.
      CHECKL = 'ABCDEFGHIYJKLMNOPQRSTUVWZ'
      CHECKS = '#?!*'
C
      DO I=1,LEN(LINE)
        IF (INDEX(' ',LINE(I:I)) == 0) THEN
          ILOW = ILOW + 1
          LENSTR(ILINE) = LENSTR(ILINE) + 1
        ELSE
          ULSTR(ILINE) = LLSTR(ILINE) + LENSTR(ILINE) - 1
          ILINE = ILINE + 1
          ILOW = 0
        END IF
        IF (ILOW == 1) LLSTR(ILINE) = I
      END DO
C
      RNUM = 0
      CNUM = 0
      DO I=1,LEN(LINE)
        LINE1 = UCC(LINE(I:I))
        IF (INDEX('.',LINE(I:I)) /= 0) RNUM = RNUM + 1
        IF ((INDEX(CHECKL,LINE1) /= 0).OR.
     $      (INDEX(CHECKS,LINE(I:I)) /= 0)) CNUM = CNUM + 1
      END DO
      IF (RNUM /= 0) GOTO 300
      IF (CNUM /= 0) GOTO 200
C
      LIM = 0
  100 CONTINUE
      LIM = LIM + 1
      BACKSPACE(BASIS)
      READ(BASIS,*,ERR=9999) (INAT(I),I=1,LIM)
      IF (INAT(LIM) /= 0) GOTO 100
      DO I=1,LIM
        IF (INAT(I) == IATOM) SUCCESS = .TRUE.
      END DO
C
  200 CONTINUE
C
      DO I=1,ILINE
        IF (LENSTR(I) /= 0) THEN
          IF (UCC(LINE(LLSTR(I):ULSTR(I))) == UCC(STRING(1:2))) THEN
             SUCCESS = .TRUE.
             GOTO 300
          END IF
        END IF
      END DO
C
  300 CONTINUE
      RETURN 
C
 9999 CONTINUE
      WRITE(OUT,1000) CLINE
      STOP
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Format statements ___
C
 1000 FORMAT(/,T2,'FINELE: ERROR READING BASIS SET INPUT LINE',I5)
C
C     O________________________________________________________________O
C
      END SUBROUTINE FINELE
      
