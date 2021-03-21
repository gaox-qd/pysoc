      SUBROUTINE GENOP
C
C     *******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
C
      INCLUDE 'fileio.h'
      INCLUDE 'molecule.h'
      INCLUDE 'moleculeD.h'
C
      INTEGER :: AX,AY,AZ,IATOM,IS,ISLAT,LSTO,MXLB,NUM,ICOU
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Initialize pointers ___
C
      GLL = 0
      GUL = 0
      LLS = 0
      ULS = 0
      LLSTO = 0
      ULSTO = 0
C
C     ___ Initialize counters ___
C
      ISLAT = 0
      IATOM = 1
      NUM = 0
      ICOU = 1
      MXLB = 0
C
C     ___ Build pointers ___
C
      DO IS=1,NSHL
        IF (IS == 1) THEN 
          GLL(1) = 1
        ELSE
          GLL(IS) = GUL(IS-1) + 1
        END IF
        GUL(IS) = GLL(IS) + PSHELL(IS,3) - 1
        LSTO = PSHELL(IS,2)
        IATOM = PSHELL(IS,1)
        MXLB = MAX(MXLB,LSTO)
        IF (IATOM == 1) THEN
          LLS(1) = 1
          ULS(1) = ULS(1) + 1
        ELSE
          LLS(IATOM) = ULS(IATOM-1) + 1
          ULS(IATOM) = IS
        END IF
C
        DO AX=LSTO,0,-1
          DO AY=LSTO-AX,0,-1
            AZ = LSTO - AX - AY
            ISLAT = ISLAT + 1
            IF (ISLAT > NSTO) GOTO 10000
            IF (ICOU == IS) THEN
              NUM = NUM + 1
            ELSE
              ICOU = IS
              NUM = 1
            END IF
          END DO
        END DO
        IF (IS == 1) THEN
          LLSTO(1) = 1
        ELSE
          LLSTO(IS) = ULSTO(IS-1) + 1
        END IF
        ULSTO(IS) = LLSTO(IS) + NUM - 1
      END DO
C
      DSHL = ((MXLB + 1)*(MXLB + 2))/2
C
      RETURN
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Error handling ___
C
10000 CONTINUE
      WRITE(OUT,"(/,T2,'ERROR FOR ORBITAL POINTER')")
      STOP
C
C     O________________________________________________________________O
C
      END SUBROUTINE GENOP
