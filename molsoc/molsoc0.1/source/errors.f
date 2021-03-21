      SUBROUTINE ERRORS 
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
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Write error messages and exit ___ 
C
      IF (DISCIN == 2) THEN
        IF (((.NOT.DIRECTA).AND.(DIRECTB)).OR.
     $      ((DIRECTA).AND.(.NOT.DIRECTB))) THEN
          IF ((MULTF /= 1).AND.(MULTS /= 3)) THEN
            WRITE(OUT,1000)
            STOP
          END IF
        END IF
      END IF
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     *** Format statements ***
C
 1000 FORMAT(/,T2,'ERRONEOUS COMBINATION OF THE COUPLE OF IDENTIFIER O',
     $       /,T2,'FOR DIRECT OR REVERSE PROCEDURE')
C
C     O________________________________________________________________O
C
      END SUBROUTINE ERRORS
