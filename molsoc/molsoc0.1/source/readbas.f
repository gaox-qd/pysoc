      SUBROUTINE READBAS
C
C     ******************************************************************
C     ***                    Read basis set                          ***
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'fileio.h'
      INCLUDE 'molecule.h'
      INCLUDE 'moleculeD.h'
C
      CHARACTER :: ELEM*2,LINE*80,ORS*1
      INTEGER :: IATOM,ILINE,IAT
      REAL :: SCALZET
C
      INTEGER ALLOCATION
      INTEGER, DIMENSION(:), ALLOCATABLE :: ILIN,INAT
C
      CHARACTER :: UCC*7
C
      LOGICAL :: SUCCESS
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Allocate local fields ___
C
      ALLOCATE(ILIN(NATOM),INAT(NATOM),STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,"(/,T2,'READBAS: ALLOCATION FAILED')")
        STOP
      END IF
C
C     ___ Find atomic basis set line ___
C
      DO IATOM=1,NATOM
C
        ILINE = 0
        ILIN(IATOM) = 0
C
   50   CONTINUE
        ILINE = ILINE + 1
        READ(BASIS,"(A80)",ERR=9999) LINE
        IF (LINE(1:1) == '#') GOTO 50
        CALL FINELE(LINE,SUCCESS,ELEMENT(IATOM),INAT,IATOM,ILINE)
        IF (SUCCESS) ILIN(IATOM) = ILINE
        CALL FINSTR(LINE,SUCCESS,'END')
        IF (SUCCESS) THEN
          IF (ILIN(IATOM) == 0) THEN
            WRITE(OUT,1000) UCC(ELEMENT(IATOM))
            STOP
          END IF
        END IF
        IF (ILIN(IATOM) == 0) GOTO 50
C
        REWIND(BASIS)
C
      END DO
C
C     ___ Read basis set ___
C
      DO IATOM=1,NAT
C
        ILINE = 0
C
  100   CONTINUE
        ILINE = ILINE + 1
        READ(BASIS,"(A80)",ERR=9999) LINE
        IF (ILIN(IATOM) /= ILINE) GOTO 100
C
C     ___ Loop over the shells ___
C
  300   CONTINUE
C
C       ___ Read shell pointers and GTF ___
C
        CALL READSHGF(IATOM,ILINE)
C
        ILINE = ILINE + 1
        READ(BASIS,"(A80)",ERR=9999) LINE
        CALL FINSTR(LINE,SUCCESS,'***')
        IF (SUCCESS) GOTO 200
        CALL FINSTR(LINE,SUCCESS,'END')
        IF (SUCCESS) GOTO 200
        BACKSPACE(BASIS)
        ILINE = ILINE - 1
        GOTO 300
C
  200   CONTINUE
C
        REWIND(BASIS)
C
      END DO
C
C     ___ Deallocate local fields ___
C
      DEALLOCATE(ILIN,INAT,STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,"(/,T2,'READBAS: DEALLOCATION FAILED')")
        STOP
      END IF
C
      RETURN
C
 9999 CONTINUE
      WRITE(OUT,"('ERROR READING BASIS SET INPUT LINE',I5)") ILINE
      STOP
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Format statements ___
C
 1000 FORMAT(/,T2,'READBAS: CANNOT FIND BASIS SET FOR THE ELEMENT ',A2)
C
C     O________________________________________________________________O
C
      END SUBROUTINE READBAS
