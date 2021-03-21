      SUBROUTINE MULDENS(COE,DMAT,LFOMO,OCCN,OPTION1)
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
      CHARACTER :: OPTION1*(*),STRING*34
C
      INTEGER :: IMAT,JMAT,DMAT,IMO,LFOMO
C
      REAL MUL,NELECTR,MC
      REAL, DIMENSION(DMAT,DMAT) :: COE
      REAL, DIMENSION(MXSTO) :: OCCN
C
      INTEGER :: ALLOCATION
      REAL, DIMENSION(:,:), ALLOCATABLE :: S,P
C
      CHARACTER :: UCC*7
      REAL :: TRAATB
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      IF (OPTION1 == 'FIRST') THEN
        STRING = ' BEFORE SVD (BIORTHOGONALIZATION) '
      ELSE IF (OPTION1 == 'SECOND') THEN
        STRING = ' AFTER SVD (BIORTHOGONALIZATION) '
      END IF
C
C     *** Allocate local fields ***
C
      ALLOCATE(S(DMAT,DMAT),P(DMAT,DMAT),STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN 
        WRITE(OUT,*) 'MULDENS: ALLOCATION FAILED'
        STOP
      END IF
C
C     *** Read Overlap matrix from tape ***
C
      CALL IOSCR(S,DMAT**2,'READ','ZERO',4)
C
      IF (GRAMS) THEN
        CALL ORTHO(S,COE,DMAT,NSPH)
        IF (OPTION1 == 'FIRST') THEN
          WRITE(OUT,"(/,T2,'GRAM-SCHMIDT ORTHOGONALIZATION')")
        END IF
      END IF
C
C     *** Calculate spin-density matrix ***
C
      P = 0.0
      DO IMAT=1,DMAT
        DO JMAT=IMAT,DMAT
          DO IMO=1,LFOMO
            P(IMAT,JMAT) = P(IMAT,JMAT) + 2.0*OCCN(IMO)*
     $                     COE(IMAT,IMO)*COE(JMAT,IMO)
          END DO
        END DO
      END DO
C
      DO IMAT=1,DMAT
        P(IMAT,IMAT) = 0.5*P(IMAT,IMAT)
      END DO
C
C     *** Calculate the number of sigma electrons ***
C
      NELECTR = TRAATB(P,S,DMAT,DMAT)
C
C     *** Write the number of spin electrons ***
C
      IF ((OPTION1 == 'FIRST').OR.(OPTION1 == 'SECOND')) THEN 
        WRITE(OUT,"(//,T2,3('*'),A34,3('*'),/)") STRING
      END IF
C
      WRITE(OUT,"(T2,'NUMBER OF ELECTRONS =', F20.9)") NELECTR
C
C     *** Deallocate local fields ***
C
      DEALLOCATE(S,P,STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN 
        WRITE(OUT,*) 'MULDENS: DEALLOCATION FAILED'
        STOP
      END IF
C
C     O________________________________________________________________O
C
      END SUBROUTINE MULDENS 
