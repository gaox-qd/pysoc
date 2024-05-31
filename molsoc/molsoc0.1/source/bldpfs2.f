      SUBROUTINE BLDPFS2(CA1,CB1,CA2,CB2,HOMOA1,HOMOB1,HOMOA2,HOMOB2,
     $                   DMAT)
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
      INTEGER :: DMAT,HOMOA1,HOMOB1,HOMOA2,HOMOB2,IMAT,JMAT
      INTEGER :: IMO,IMOA,IMOB
C
      REAL :: PRODAM,PRODBM,SMOP
      REAL, DIMENSION(DMAT,DMAT) :: CA1,CA2,CB1,CB2
C
      INTEGER ALLOCATION
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
        WRITE(OUT,*) 'BLDPFS2: ALLOCATION FAILED'
        STOP
      END IF
C
      CALL IOSCR(SMOA,DMAT**2,'READ','ZERO',5)
      CALL IOSCR(SMOB,DMAT**2,'READ','ZERO',6)
C
      PRODAT = 1.0
      IF (HOMOA1.GT.0) THEN
        DO IMOA=1,HOMOA1
          SMOP = SMOA(IMOA,IMOA)
          IF (ABS(SMOP).LE.TOLS) SMOP = 0.0
          PRODAT = PRODAT*SMOP
        END DO
      END IF
C
      PRODBT = 1.0
      IF (HOMOB2.GT.0) THEN
        DO IMOB=1,HOMOB2
          SMOP = SMOB(IMOB,IMOB)
          IF (ABS(SMOP).LE.TOLS) SMOP = 0.0
          PRODBT = PRODBT*SMOP
        END DO
      END IF
C
      PABA = 0.0
      IF (HOMOA1.GT.0) THEN
        DO IMO=1,HOMOA1
          PRODAM = 1.0
          DO IMOA=1,HOMOA1
            IF (IMOA.NE.IMO) THEN
              SMOP = SMOA(IMOA,IMOA)
              IF (ABS(SMOP).LE.TOLS) SMOP = 0.0
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
        WRITE(OUT,*) 'BLDPFS2: HOMOA1 LE ZERO'
        STOP
      END IF
C
      PABA(:,:) = PABA(:,:)*PRODBT
C
      PABB = 0.0
      IF (HOMOB2 > 0) THEN
        DO IMO=1,HOMOB2
          PRODBM = 1.0
          DO IMOB=1,HOMOB2
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
        WRITE(OUT,*) 'BLDPFS2: HOMOB2 LE ZERO'
        STOP
      END IF
C
      PABB(:,:) = PABB(:,:)*PRODAT
C
C     ___ Write MOs overlap on tape ___
C
      CALL IOSCR(PABA,DMAT**2,'WRITE','PMAT',1)
      CALL IOSCR(PABB,DMAT**2,'WRITE','PMAT',2)
      CALL IOSCR(CB1(:,HOMOB1),DMAT,'WRITE','PMAT',3)
      CALL IOSCR(CA2(:,HOMOA2),DMAT,'WRITE','PMAT',4)
C
C     ___ Deallocate local fields ___
C
      DEALLOCATE(SMOA,SMOB,PABA,PABB,STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN 
        WRITE(OUT,*) 'BLDPFS2: DEALLOCATION FAILED'
        STOP
      END IF
C
C     O________________________________________________________________O
C
      END SUBROUTINE BLDPFS2
