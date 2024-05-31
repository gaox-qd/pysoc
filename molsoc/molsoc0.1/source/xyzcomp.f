      SUBROUTINE XYZCOMP 
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'molecule.h'
      INCLUDE 'moleculeD.h'
      INCLUDE 'fileio.h'
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C
      ENER(1) = 0.5*ENER(1); ENERD(1) = 0.5*ENERD(1)
      ENER(2) = 0.5*ENER(2); ENERD(2) = 0.5*ENERD(2)
      ENER(3) = 0.5*ENER(3); ENERD(3) = 0.5*ENERD(3)
C
      IF ((MULTF == 1).AND.(MULTS == 1)) THEN
      WRITE(OUT,1000)
      WRITE(OUT,3000) ENERD(1),ENER(1),
     $                ENERD(2),ENER(2),
     $                ENERD(3),ENER(3)
C
      IF (SOCTWOM) THEN
        WRITE(OUT,2000)
        WRITE(OUT,3000) ENERTWOAU(1),ENERTWO(1),
     $                  ENERTWOAU(2),ENERTWO(2),
     $                  ENERTWOAU(3),ENERTWO(3)
      END IF
C
      WRITE(OUT,4000)
      WRITE(OUT,3000) ENERD(1)+ENERTWOAU(1),
     $                ENER(1)+ENERTWO(1),
     $                ENERD(2)+ENERTWOAU(2),
     $                ENER(2)+ENERTWO(2),
     $                ENERD(3)+ENERTWOAU(3),
     $                ENER(3)+ENERTWO(3)
C
      ELSE
C     IF ((MULTF == 1).AND.(MULTS == 1)) GOTO 100
C
C     ___ Calculate the matrix elements ___
C
        CALL MATRELEM
C
      END IF
  100 CONTINUE
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Format statements ___
C
 1000 FORMAT(/,T2,'ONE ELECTRON SOC CONTRIBUTION:')
 2000 FORMAT(/,T2,'TWO ELECTRON SOC CONTRIBUTION:')
 3000 FORMAT(T2,42('-'),
     $       /,T2,15X,'A.U.',15X,'CM-1',/,
     $       /,T2,'X ',2F20.9,
     $       /,T2,'Y ',2F20.9,
     $       /,T2,'Z ',2F20.9)
 4000 FORMAT(/,T2,'TOTAL:')
C
C     O________________________________________________________________O
C
      END SUBROUTINE XYZCOMP
