      SUBROUTINE TURBORD(MATRIX,DMAT,OPTION)
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
C
      INCLUDE 'fileio.h'
      INCLUDE 'molecule.h'
      INCLUDE 'moleculeD.h'
C
      CHARACTER*(*) :: OPTION
C
      INTEGER :: IS,L,IMO,DMAT
      INTEGER, DIMENSION(10) :: POINT
C
      REAL, DIMENSION(DMAT,DMAT) :: MATRIX,NEWMAT
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      NEWMAT = 0.0
      POINT = 0
C
      IF (OPTION.EQ.'SPHERICAL') THEN
C
        DO IS=1,NSHL
          L = PSHELL(IS,2)
          IF (L == 0) THEN
            DO IMO=1,NSPH
              NEWMAT(SLL(IS),IMO) = MATRIX(SLL(IS),IMO)
            END DO
          ELSE IF (L == 1) THEN
            POINT(1) = SLL(IS)
            POINT(2) = SLL(IS) + 1
            POINT(3) = SUL(IS)
            DO IMO=1,NSPH
              NEWMAT(POINT(3),IMO) = MATRIX(POINT(1),IMO)
              NEWMAT(POINT(1),IMO) = MATRIX(POINT(2),IMO)
              NEWMAT(POINT(2),IMO) = MATRIX(POINT(3),IMO)
            END DO
          ELSE IF (L == 2) THEN
            POINT(1) = SLL(IS)
            POINT(2) = SLL(IS) + 1
            POINT(3) = SLL(IS) + 2
            POINT(4) = SLL(IS) + 3
            POINT(5) = SUL(IS)
            DO IMO=1,NSPH
              NEWMAT(POINT(3),IMO) = MATRIX(POINT(1),IMO)
              NEWMAT(POINT(4),IMO) = MATRIX(POINT(2),IMO)
              NEWMAT(POINT(2),IMO) = MATRIX(POINT(3),IMO)
              NEWMAT(POINT(1),IMO) = MATRIX(POINT(4),IMO)
              NEWMAT(POINT(5),IMO) = MATRIX(POINT(5),IMO)
            END DO
          ELSE IF (L == 3) THEN
            WRITE(OUT,*) 'F CARTESIAN CONVERSION AOs IS NOT AVAILABLE' 
            STOP
          ELSE 
            WRITE(OUT,*) 'TURBORD: NO ANGULAR MOMENTUM FOUND'
            STOP
          END IF
        END DO
C
      ELSE IF (OPTION == 'CARTESIAN') THEN
C
        WRITE(OUT,*) 'CARTESIAN CONVERSION AOs IS NOT AVAILABLE'
        STOP
C
      ELSE 
C
        WRITE(OUT,*) 'TURBORD: UNKNOWN OPTION'
        STOP
C
      END IF
C
      MATRIX = NEWMAT
C
C     O________________________________________________________________O
C
      END SUBROUTINE TURBORD
