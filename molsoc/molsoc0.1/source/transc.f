      SUBROUTINE TRANSC
C
C     *****************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
C
      INCLUDE 'math.h'
      INCLUDE 'moleculeD.h'
C
      INTEGER :: I,IC,IS,J,K,L,LX,LY,LZ,M,ABM
C
      REAL :: FITERM,STERM,TTERM,FOTERM,STOT 
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      DO L=0,MXL
        IC = 0
        DO LX=L,0,-1
          DO LY=L-LX,0,-1
            LZ = L - LX - LY
            IC = IC + 1
            DO M=-L,L
              ABM = ABS(M)
              IS = L + M + 1
              IF ((ABM.GT.L).OR.(LZ.LT.0)) GOTO 10
              FITERM = SQRT(FAC(2*LX)*FAC(2*LY)*FAC(2*LZ)*FAC(L)*
     $                      FAC(L-ABM)/(FAC(2*L)*FAC(LX)*FAC(LY)*
     $                      FAC(LZ)*FAC(L+ABM)))/((2**L)*FAC(L))
              J = (LX + LY - ABM)/2
              STOT = 0.0
              IF ((2*J).EQ.(LX + LY - ABM)) THEN
                DO I=0,(L-ABM)/2
                  STERM = 0.0
                  IF ((J.GE.0).AND.(J.LE.I)) THEN
                    STERM = (FAC(L)/(FAC(I)*FAC(L-I)))*
     $                      (FAC(2*L-2*I)*((-1)**I)/FAC(L-ABM-2*I))*
     $                      (FAC(I)/(FAC(J)*FAC(I-J)))
                    DO K=0,J
                      TTERM = 0.0
                      IF (((LX-2*K).GE.0).AND.((LX-2*K).LE.ABM)) THEN
                        TTERM = (FAC(J)/(FAC(K)*FAC(J-K)))*(FAC(ABM)/
     $                          (FAC(LX-2*K)*FAC(ABM+2*K-LX)))
                        FOTERM = 0.0
                        IF ((M.EQ.0).AND.((MOD(LX,2).EQ.0)))
     $                  FOTERM = (-1)**(K - LX/2)
                        IF (((M.GT.0).AND.
     $                      (MOD(ABS(ABM-LX),2).EQ.0)).OR.
     $                      ((M.LT.0).AND.(MOD(ABS(ABM-LX),2).EQ.1)))
     $                  FOTERM = (-1)**((2*K+ABM-LX)/2)*SQRT(2.0)
                        STOT = STOT + STERM*TTERM*FOTERM
                      END IF
                    END DO
                  END IF
                END DO
              END IF
              TRMSC(IS,IC,L) = FITERM*STOT
   10         CONTINUE
            END DO
          END DO
        END DO
      END DO
C
C     O________________________________________________________________O
C
      END SUBROUTINE TRANSC
