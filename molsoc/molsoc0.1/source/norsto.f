      SUBROUTINE NORSTO
C
C     ******************************************************************
C     ***           Normalization of Slater type orbitals.           ***
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'math.h'
      INCLUDE 'molecule.h'
      INCLUDE 'moleculeD.h'
C
      INTEGER :: AX,AY,AZ,IATOM,IG,IS,ISLAT,JG,L
      REAL :: FACT,NOR
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      NCSTO = 0.0
      ISLAT = 0
C
      IF (TDTB) THEN
        NCSTO = 1.0
      ELSE
        DO IATOM=1,NATOM
          DO IS=LLS(IATOM),ULS(IATOM)
            L = PSHELL(IS,2)
            FACT = PI**1.5/2**L
            NOR = 0.0
            DO IG=GLL(IS),GUL(IS)
              DO JG=GLL(IS),GUL(IS)
                NOR = NOR + GCC(IG)*GCC(JG)*FACT/
     $                      (ZET(IG) + ZET(JG))**(L+1.5)
              END DO
            END DO  
C
            DO AX=L,0,-1
              DO AY=L-AX,0,-1
                AZ = L - AX - AY
                ISLAT = ISLAT + 1
                NCSTO(ISLAT) = 1.0/SQRT(NOR*DFAC(2*AX-1)*DFAC(2*AY-1)*
     $                                                 DFAC(2*AZ-1))
              END DO
            END DO
C
          END DO  
        END DO  
      END IF
C
C     O________________________________________________________________O
C
      END SUBROUTINE NORSTO
