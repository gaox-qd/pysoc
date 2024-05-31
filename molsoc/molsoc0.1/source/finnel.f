      SUBROUTINE FINNEL
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'molecule.h'
      INCLUDE 'moleculeD.h'
C
      INTEGER :: IATOM
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Calculate number of electrons ___
C
      NELEC = 0
      DO IATOM=1,NATOM
        NELEC = NELEC + INT(ZATOM(IATOM))
      END DO
      NELEC = NELEC - CHARGE
C
C     ___ Set threshold for integrals ___
C
      INTHRESH2 = 1.0E-10
      INTHRESH = 1.0E-10
C
C     O________________________________________________________________O
C
      END SUBROUTINE FINNEL
