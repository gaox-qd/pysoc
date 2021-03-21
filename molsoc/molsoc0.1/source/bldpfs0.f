      SUBROUTINE BLDPFS0(CA1,CB1,CA2,CB2,HOMOA,HOMOAS,HOMOB,DMAT)
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
      REAL, PARAMETER :: TOLS = 1.0E-10
C
      INTEGER :: DMAT,HOMOA,HOMOB,IMAT,JMAT,HOMOAM,HOMOAS
      INTEGER :: IMO,IMOA,IMOB
C
      REAL :: PRODAM,PRODBM,SMOP
      REAL, DIMENSION(DMAT,DMAT) :: CA1,CA2,CB1,CB2
C
      INTEGER :: ALLOCATION
      REAL, DIMENSION(:,:), ALLOCATABLE :: SMOA,SMOB,PABA,PABB
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      IF (.NOT.SOCTWOM) RETURN
C
C     ___ Allocate local fields ___
C
      ALLOCATE(SMOA(DMAT,DMAT),SMOB(DMAT,DMAT),PABA(DMAT,DMAT),
     $         PABB(DMAT,DMAT),STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'BLDPFS0: ALLOCATION FAILED'
        STOP
      END IF
C
      IF (HOMOA /= HOMOAS) THEN
        WRITE(OUT,*) 'BLDPFS0: HOMOAF NE HOMOAS'
        STOP
      END IF
C
      HOMOAM = HOMOA - 1
C
      CALL IOSCR(SMOA,DMAT**2,'READ','ZERO',5)
      CALL IOSCR(SMOB,DMAT**2,'READ','ZERO',6)
C
      PRODAT = 1.0
      IF (HOMOAM > 0) THEN
        DO IMOA=1,HOMOAM
          SMOP = SMOA(IMOA,IMOA) 
          IF (ABS(SMOP) <= TOLS) SMOP = 0.0
          PRODAT = PRODAT*SMOP
        END DO
      END IF
C
      PRODBT = 1.0
      IF (HOMOB > 0) THEN
        DO IMOB=1,HOMOB
          SMOP = SMOB(IMOB,IMOB)
          IF (ABS(SMOP) <= TOLS) SMOP = 0.0
          PRODBT = PRODBT*SMOP
        END DO
      END IF
C
      PABA = 0.0
      IF (HOMOAM > 0) THEN
        DO IMO=1,HOMOAM
          PRODAM = 1.0
          DO IMOA=1,HOMOAM
            IF (IMOA /= IMO) THEN
              SMOP = SMOA(IMOA,IMOA)
              IF (ABS(SMOP) <= TOLS) SMOP = 0.0
              PRODAM = PRODAM*SMOP 
            END IF
          END DO
          DO IMAT=1,DMAT
            DO JMAT=1,DMAT
              PABA(IMAT,JMAT) = PABA(IMAT,JMAT) + 
     $                          CA1(IMAT,IMO)*CA2(JMAT,IMO)*PRODAM
            END DO
          END DO
        END DO
      ELSE
        WRITE(OUT,*) 'BLDPFS0: HOMOA-1 LE ZERO'
        STOP
      END IF
C
      PABA(:,:) = PABA(:,:)*PRODBT
C
      PABB = 0.0
      IF (HOMOB > 0) THEN
        DO IMO=1,HOMOB
          PRODBM = 1.0
          DO IMOB=1,HOMOB
            IF (IMOB /= IMO) THEN
              SMOP = SMOB(IMOB,IMOB)
              IF (ABS(SMOP) <= TOLS) SMOP = 0.0
              PRODBM = PRODBM*SMOP
            END IF
          END DO
          DO IMAT=1,DMAT
            DO JMAT=1,DMAT
              PABB(IMAT,JMAT) = PABB(IMAT,JMAT) + 
     $                          CB1(IMAT,IMO)*CB2(JMAT,IMO)*PRODBM
            END DO
          END DO
        END DO
      ELSE 
        WRITE(OUT,*) 'BLDPFS: HOMOB LE ZERO'
        STOP
      END IF
C
      PABB(:,:) = PABB(:,:)*PRODAT
C
C     ___ Write density matrices on tape ___
C
      CALL IOSCR(PABA,DMAT**2,'WRITE','PMAT',1)
      CALL IOSCR(PABB,DMAT**2,'WRITE','PMAT',2)
      CALL IOSCR(CA1(:,HOMOA),DMAT,'WRITE','PMAT',3)
      CALL IOSCR(CA2(:,HOMOAS),DMAT,'WRITE','PMAT',4)
C
C     ___ Deallocate local fields ___
C
      DEALLOCATE(SMOA,SMOB,PABA,PABB,STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN 
        WRITE(OUT,*) 'BLDPFS0: DEALLOCATION FAILED'
        STOP
      END IF
C
C     O________________________________________________________________O
C
      END SUBROUTINE BLDPFS0
