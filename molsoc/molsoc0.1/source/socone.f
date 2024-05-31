      SUBROUTINE SOCONE
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'fileio.h'
      INCLUDE 'moleculeD.h'
      INCLUDE 'molecule.h'
C
      INTEGER :: IATOM,ISTO,JSTO,U
      REAL :: CPUTOLD
C
      INTEGER :: ALLOCATION
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: MATRIX,NMATRIX
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Allocate local fields ___
C
      ALLOCATE(MATRIX(NSTO,NSTO,3),NMATRIX(NSTO,NSTO,3),STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'SOCONE: ALLOCATION FAILED'
        STOP
      END IF
C
C     ___ Get CPU time ___
C
      CALL CPUTIME
      CPUTOLD = CPUT
C
C     ___ Calculate SOC matrix ___
C
      NMATRIX = 0.0
      DO IATOM=1,NATOM
        CALL SOCINT(MATRIX,C(:,IATOM))
        NMATRIX = NMATRIX + ZEFF(IATOM)*MATRIX
      END DO
C
C     ___ Get CPU time ___
C
      CALL CPUTIME
      WRITE(OUT,2000) CPUT - CPUTOLD
C
C     ___ Store SOC matrix on tape ___
C
      DO U=1,3
        CALL IOSCR(NMATRIX(:,:,U),NSTO**2,'WRITE','ZERO',U)
      END DO
C
      WRITE(SOCA,999)
      DO U=1,3
        IF (U == 1) WRITE(SOCA,*) 'X COMPONENT'
        IF (U == 2) WRITE(SOCA,*) 'Y COMPONENT'
        IF (U == 3) WRITE(SOCA,*) 'Z COMPONENT'
        DO ISTO=1,NSTO
          DO JSTO=1,NSTO
            WRITE(SOCA,1100) JSTO,ORSYM(JSTO),
     $                      ISTO,ORSYM(ISTO),NMATRIX(JSTO,ISTO,U)
          END DO
        END DO
      END DO
C
C     ___ Deallocate local fields ___
C
      DEALLOCATE(MATRIX,NMATRIX,STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'SOCONE: DEALLOCATION FAILED'
        STOP
      END IF 
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Format statements ___
C
  999 FORMAT(//,T2,'ONE-ELECTRON SOC INTEGRALS <PHI|SOC|PHI>',
     $       /,T2,45('-'))
 1000 FORMAT(I7,4X,I7,10X,F20.9)
 1100 FORMAT(I5,A5,3X,I5,A5,F20.9)
 2000 FORMAT(//,T2,'ONE-ELECTRON SOC INTEGRALS CALCULATED IN',F12.6,
     $          X,'SEC.') 
C
C     O________________________________________________________________O
C
      END SUBROUTINE SOCONE
