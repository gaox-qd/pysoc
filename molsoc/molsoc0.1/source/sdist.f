      SUBROUTINE SDIST
C
C     ******************************************************************
C     ***      Calculate the square of the distance matrix           ***
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
C
      INCLUDE 'molecule.h'
      INCLUDE 'moleculeD.h'
C
      INTEGER :: I,IATOM,JATOM
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Calculate the square of the distance matrix ___
C
      DO IATOM=1,NATOM
        SDISTANCE(IATOM,IATOM) = 0.0
        DO JATOM=IATOM+1,NATOM
          SDISTANCE(IATOM,JATOM) = (C(1,IATOM) - C(1,JATOM))**2 + 
     $                             (C(2,IATOM) - C(2,JATOM))**2 +
     $                             (C(3,IATOM) - C(3,JATOM))**2
          SDISTANCE(JATOM,IATOM) = SDISTANCE(IATOM,JATOM)
        END DO
      END DO  
C
C     O________________________________________________________________O
C
      END SUBROUTINE SDIST
