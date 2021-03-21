      SUBROUTINE SOCFST0
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
      INTEGER :: IATOM,JATOM,KATOM,LATOM,ISHL,JSHL,KSHL,LSHL,LA,LB,LC,
     $           LD,IMO,IMOA,IMOB
C
      REAL :: CPUTOLD
      REAL, DIMENSION(3) :: TWOSOA,TWOSOB
      REAL, DIMENSION(NSTO) :: CA1,CA2
C
      INTEGER :: ALLOCATION
      REAL, DIMENSION(:,:), ALLOCATABLE :: PABA,PABB
      REAL, DIMENSION(:,:,:,:,:), ALLOCATABLE :: BLOCKS
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Calculate MO SOC TWO electron integrals ___
C
      IF (.NOT.SOCTWOM) RETURN
C
C     ___ Allocate local fields ___
C
      ALLOCATE(PABA(NSTO,NSTO),PABB(NSTO,NSTO),
     $         BLOCKS(DSHL,DSHL,DSHL,DSHL,3),STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN 
        WRITE(OUT,*) 'SOCFST0: ALLOCATION FAILED'
        STOP
      END IF
C
      CALL IOSCR(PABA,NSTO**2,'READ','PMAT',1)
      CALL IOSCR(PABB,NSTO**2,'READ','PMAT',2)
      CALL IOSCR(CA1,NSTO,'READ','PMAT',3)
      CALL IOSCR(CA2,NSTO,'READ','PMAT',4)
C
C     ___ Get CPU time ___
C
      CALL CPUTIME
      CPUTOLD = CPUT
C
      ENERTWO = 0.0
      ENERTWOAU = 0.0
C
      BLOCKS = 0.0
      TWOSOA = 0.0
      TWOSOB = 0.0
C
      DO IATOM=1,NATOM
        DO JATOM=IATOM,NATOM
          CALL ATOMDIS(IATOM,JATOM,RAB,RSQAB)
          AAT(:) = C(:,IATOM); BAT(:) = C(:,JATOM)
          DO ISHL=LLS(IATOM),ULS(IATOM)
            DO JSHL=MAX(ISHL,LLS(JATOM)),ULS(JATOM)
              DO KATOM=1,NATOM
                DO LATOM=KATOM,NATOM
                  CALL ATOMDIS(KATOM,LATOM,RCD,RSQCD)
                  CAT(:) = C(:,KATOM); DAT(:) = C(:,LATOM)
                  DO KSHL=LLS(KATOM),ULS(KATOM)
                    DO LSHL=MAX(KSHL,LLS(LATOM)),ULS(LATOM)
                      LA = PSHELL(ISHL,2)
                      LB = PSHELL(JSHL,2)
                      LC = PSHELL(KSHL,2)
                      LD = PSHELL(LSHL,2)
                      CALL RR2SOCIN(LA,LB,LC,LD,ISHL,JSHL,KSHL,LSHL,
     $                              BLOCKS)
                      CALL BLDSOC0(ISHL,JSHL,KSHL,LSHL,TWOSOA,TWOSOB,
     $                             CA1,CA2,PABA,PABB,BLOCKS)
C 
                    END DO
                  END DO
                END DO
              END DO 
            END DO
          END DO
        END DO
      END DO
C
      ENERTWO = TWOSOA + TWOSOB
      ENERTWOAU = -0.5*FINESTRC**2*ENERTWO*1.0/2.0
      ENERTWO = -CMM1*0.5*FINESTRC**2*ENERTWO*1.0/2.0 
C
C     ___ Get CPU time ___
C
      CALL CPUTIME
      WRITE(OUT,1000) CPUT - CPUTOLD
C
C     ___ Deallocate local fields ___
C
      DEALLOCATE(PABA,PABB,BLOCKS,STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN 
        WRITE(OUT,*) 'SOCFST0: DEALLOCATION FAILED'
        STOP
      END IF
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Format statements ___
C
 1000 FORMAT(/,T2,'TWO-ELECTRON SOC INTEGRALS CALCULATED IN',
     $         F12.6,X,'SEC.')
C
C     O________________________________________________________________O
C
      END SUBROUTINE SOCFST0
