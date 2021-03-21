      SUBROUTINE ORTHO(S,CA,DMAT,DIM)
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      INTEGER :: DMAT,DIM
      INTEGER :: IMAT,JMAT,KMAT,LMAT
C
      REAL :: PROJEC,SUM
      REAL, DIMENSION(DMAT,DMAT) :: S,CA
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Cram-Schmidt Orthogonalization ___
C
      DO IMAT=1,DIM
C
        DO JMAT=1,IMAT-1
          PROJEC = 0.0
          DO KMAT=1,DMAT
            DO LMAT=1,DMAT
              PROJEC = PROJEC + CA(LMAT,JMAT)*CA(KMAT,IMAT)*S(LMAT,KMAT)
            END DO
          END DO
          DO KMAT=1,DMAT
            CA(KMAT,IMAT) = CA(KMAT,IMAT) - PROJEC*CA(KMAT,JMAT)
          END DO
        END DO
C
        SUM = 0.0
        DO JMAT=1,DMAT
          DO KMAT=1,DMAT
            SUM = SUM + CA(JMAT,IMAT)*CA(KMAT,IMAT)*S(KMAT,JMAT)
          END DO
        END DO
        SUM = 1.0/SQRT(SUM)
        DO KMAT=1,DMAT
          CA(KMAT,IMAT) = CA(KMAT,IMAT)*SUM
        END DO
C
      END DO
C
C     O_______________________________________________________________O
C
      END SUBROUTINE ORTHO
