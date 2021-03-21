      SUBROUTINE WRITOUT
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
      INTEGER :: IATOM,I,ISHELL,IG,NCO
      REAL :: TRANS1,TRANS2
C
      CHARACTER :: UCC*7
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      WRITE(OUT,1000) CHARGE,MULTF,MULTS,NELEC
C
      WRITE(OUT,1100)
      WRITE(OUT,2100)
      DO IATOM=1,NAT
        WRITE(OUT,1200) UCC(ELEMENT(IATOM)(1:2)),
     $                  (C(I,IATOM)/BOHR,I=1,3),
     $                  ZATOM(IATOM),ZEFF(IATOM)
      END DO
      WRITE(OUT,1300)
      WRITE(OUT,2100)
      DO IATOM=1,NAT
        WRITE(OUT,1200) UCC(ELEMENT(IATOM)(1:2)),
     $                  (C(I,IATOM),I=1,3),
     $                  ZATOM(IATOM),ZEFF(IATOM)
      END DO
C
      WRITE(OUT,1400) NAT,NSHL,NGTF,NSTO
      WRITE(OUT,"(//,T2,'BASIS SET',/)")
C
C     ___ Write basis sets ___
C
      DO IATOM=1,NATOM
        WRITE(OUT,"(T2,A2)") UCC(ELEMENT(IATOM)(1:2))
        DO ISHELL=LLS(IATOM),ULS(IATOM)
C          print *, "ISHELL", IATOM, ISHELL, LLS(IATOM),ULS(IATOM)
          NCO = ULSTO(ISHELL) - LLSTO(ISHELL) + 1
          WRITE(OUT,1900) PSHELL(ISHELL,1),PSHELL(ISHELL,2),
     $                    PSHELL(ISHELL,3),NCO
C          print *, "shell", PSHELL(ISHELL,1),PSHELL(ISHELL,2), 
C     $                    PSHELL(ISHELL,3),NCO
          DO IG=GLL(ISHELL),GUL(ISHELL)
            WRITE(OUT,2000) ZET(IG),GCC(IG)
C            print *, "basis", ZET(IG),GCC(IG)
          END DO
        END DO
      END DO
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Format statements ___
C
 1000 FORMAT(/,T2,'MOLECULAR CHARGE                 = ',I5,
     $       /,T2,'MULTIPLICITY OF THE FIRST STATE  = ',I5,
     $       /,T2,'MULTIPLICITY OF THE SECOND STATE = ',I5,
     $       /,T2,'NUMBER OF ELECTRONS              = ',I5)
 1100 FORMAT(//,T2,'GEOMETRY IN ANGSTROM',/)
 1200 FORMAT(T3,A2,3X,3F12.6,F6.1,F8.3)
 1300 FORMAT(/,T2,'GEOMETRY IN BOHR',/)
 1400 FORMAT(//,T2,'NUMBER OF ATOMS                  = ',I6,
     $        /,T2,'NUMBER OF SHELLS                 = ',I6,
     $        /,T2,'NUMBER OF GTFs                   = ',I6,
     $        /,T2,'NUMBER OF STOs                   = ',I6)
 1900 FORMAT(T2,4I5)
 2000 FORMAT(T5,2F20.9)
 2100 FORMAT(T2,'ATOM',10X,'X',11X,'Y',11X,'Z',6X,'Z_A',3X,'Z_eff',/)
C
C     O_______________________________________________________________O
C
      END SUBROUTINE WRITOUT
