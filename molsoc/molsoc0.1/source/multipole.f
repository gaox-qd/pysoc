      SUBROUTINE MULTIPOLE(CA1,CB1,CA2,CB2,HOMOA1,HOMOB1,HOMOA2,HOMOB2,
     $                     DMAT)
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
      REAL, PARAMETER :: TOL = 1.0E-10
C
      INTEGER :: DMAT,DMOM
      INTEGER :: IMAT,JMAT,IMO,IATOM,U,ILU,UX,UY,UZ,I
      INTEGER :: HOMOA1,HOMOB1,HOMOA2,HOMOB2
C
      REAL, DIMENSION(DMAT,DMAT) :: CA1,CB1,CA2,CB2
C
      INTEGER :: ALLOCATION
      REAL, DIMENSION(:), ALLOCATABLE :: EM1,EM2,NM,DEB
      REAL, DIMENSION(:,:), ALLOCATABLE :: MATRIX,P
C
      REAL :: TRAATB
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      IF (LM == 0) RETURN
C
      DMOM = (LM**3 - LM)/6 + LM**2 + 2*LM + 1
C
C     *** Allocate local fields ***
C
      ALLOCATE(MATRIX(DMAT,DMAT),EM1(2:DMOM),EM2(2:DMOM),
     $         NM(2:DMOM),P(DMAT,DMAT),DEB(LM),STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'MULTIPOLE: ALLOCATION FAILED'
        STOP
      END IF
C
      DEB(1) = DEBYE
      IF (LM > 1) THEN
        DO I=2,LM
          DEB(I) = DEB(I-1)/BOHR
        END DO
      END IF
C
      DO U=2,DMOM
C
C     *** Read moment integrals from tape ***
C
        CALL IOSCR(MATRIX,NSTO**2,'READ','MOME',U)
C
C     *** Initialization ***
C
        EM1(U) = 0.0
        EM2(U) = 0.0
        NM(U) = 0.0
C
C     *** Alpha moments I state ***
C
        P = 0.0
        DO IMAT=1,DMAT
          DO JMAT=IMAT,DMAT
            DO IMO=1,HOMOA1
              P(IMAT,JMAT) = P(IMAT,JMAT) + 2.0*
     $                       CA1(IMAT,IMO)*CA1(JMAT,IMO)
            END DO
          END DO
        END DO
C
        DO IMAT=1,DMAT
          P(IMAT,IMAT) = 0.5*P(IMAT,IMAT)
        END DO
C
        EM1(U) = -TRAATB(P,MATRIX,DMAT,DMAT)
C
C     *** Beta moments I state ***
C
        P = 0.0
        DO IMAT=1,DMAT
          DO JMAT=IMAT,DMAT
            DO IMO=1,HOMOB1
              P(IMAT,JMAT) = P(IMAT,JMAT) + 2.0*
     $                       CB1(IMAT,IMO)*CB1(JMAT,IMO)
            END DO
          END DO
        END DO
C
        DO IMAT=1,DMAT
          P(IMAT,IMAT) = 0.5*P(IMAT,IMAT)
        END DO
C
        EM1(U) = EM1(U) - TRAATB(P,MATRIX,DMAT,DMAT)
C
C     *** Alpha moments II state ***
C
        P = 0.0
        DO IMAT=1,DMAT
          DO JMAT=IMAT,DMAT
            DO IMO=1,HOMOA2
              P(IMAT,JMAT) = P(IMAT,JMAT) + 2.0*
     $                       CA2(IMAT,IMO)*CA2(JMAT,IMO)
            END DO
          END DO
        END DO
C
        DO IMAT=1,DMAT
          P(IMAT,IMAT) = 0.5*P(IMAT,IMAT)
        END DO
C
        EM2(U) = -TRAATB(P,MATRIX,DMAT,DMAT)
C
C     *** Beta moments II state ***
C
        P = 0.0
        DO IMAT=1,DMAT
          DO JMAT=IMAT,DMAT
            DO IMO=1,HOMOB2
              P(IMAT,JMAT) = P(IMAT,JMAT) + 2.0*
     $                       CB2(IMAT,IMO)*CB2(JMAT,IMO)
            END DO
          END DO
        END DO
C
        DO IMAT=1,DMAT
          P(IMAT,IMAT) = 0.5*P(IMAT,IMAT)
        END DO
C
        EM2(U) = EM2(U) - TRAATB(P,MATRIX,DMAT,DMAT)
C
      END DO
C
      DO ILU=1,LM
        DO UX=ILU,0,-1
          DO UY=ILU-UX,0,-1
            UZ = ILU - UX - UY
            DO IATOM=1,NATOM
              NM(COP(UX,UY,UZ)) = NM(COP(UX,UY,UZ)) + ZATOM(IATOM)*
     $                            C(1,IATOM)**UX*
     $                            C(2,IATOM)**UY*
     $                            C(3,IATOM)**UZ
            END DO
          END DO
        END DO
      END DO
C
      DO ILU=1,LM
        DO UX=ILU,0,-1
          DO UY=ILU-UX,0,-1
            UZ = ILU - UX - UY
            NM(COP(UX,UY,UZ)) = NM(COP(UX,UY,UZ))*DEB(ILU)
            EM1(COP(UX,UY,UZ)) = EM1(COP(UX,UY,UZ))*DEB(ILU)
            EM2(COP(UX,UY,UZ)) = EM2(COP(UX,UY,UZ))*DEB(ILU)
          END DO
        END DO
      END DO
C
      WRITE(OUT,1000)
      WRITE(OUT,2000) 
      WRITE(OUT,3000)
      DO U=2,DMOM
        IF (U.EQ.2) WRITE(OUT,"(T2,'DIPOLE')")
        IF (U.EQ.5) WRITE(OUT,"(T2,'QUADRUPOLE')")
        IF (U.EQ.11) WRITE(OUT,"(T2,'OCTUPOLE')")
        WRITE(OUT,4000) MOMSYM(U),NM(U),EM1(U),EM1(U)+NM(U)
      END DO
      WRITE(OUT,5000)
      WRITE(OUT,2000) 
      WRITE(OUT,3000)
      DO U=2,DMOM
        IF (U.EQ.2) WRITE(OUT,"(T2,'DIPOLE')")
        IF (U.EQ.5) WRITE(OUT,"(T2,'QUADRUPOLE')")
        IF (U.EQ.11) WRITE(OUT,"(T2,'OCTUPOLE')")
        WRITE(OUT,4000) MOMSYM(U),NM(U),EM2(U),EM2(U)+NM(U)
      END DO
C    
C     *** Deallocate local fields **
C
      DEALLOCATE(MATRIX,EM1,EM2,NM,P,DEB,STAT=ALLOCATION)
      IF (ALLOCATION /= 0) THEN
        WRITE(OUT,*) 'MULTIPOLE: DEALLOCATION FAILED'
        STOP
      END IF
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     *** Format statements ***
C
 1000 FORMAT(/,T2,'MOMENTS IN DEBYES OF THE I STATE:')
 2000 FORMAT(T2,45("-"))
 3000 FORMAT(T2,10X,'NUCLEAR',5X,'ELECTRONIC',6X,'TOTAL',/)
 4000 FORMAT(T2,A3,3F14.6)
 5000 FORMAT(/,T2,'MOMENTS IN DEBYES OF THE II STATE:')
C
C     O________________________________________________________________O
C
      END SUBROUTINE MULTIPOLE 
