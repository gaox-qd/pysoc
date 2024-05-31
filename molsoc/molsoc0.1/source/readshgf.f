      SUBROUTINE READSHGF(ATOM,ILINE)
C
C     ******************************************************************
C     ***        Read shell pointers and gaussian functions         ***
C     ******************************************************************
C
      IMPLICIT NONE
C
      INCLUDE 'parameter.h'
      INCLUDE 'fileio.h'
      INCLUDE 'molecule.h'
      INCLUDE 'moleculeD.h'
C
      CHARACTER :: ORS*2
C
      INTEGER :: PSH,ILINE,IPSHA,IPSHB,ATOM,IG,NGTF1,NGTF2
      INTEGER :: ICA,ICB,AX,AY,AZ,I,J
C
      REAL :: SCALZET
C
      CHARACTER :: UCC*7
C
      LOGICAL :: TWOSH
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
      TWOSH = .FALSE.
C
      ILINE = ILINE + 1
      READ(BASIS,*,ERR=9990) ORS,PSH,SCALZET
C      print *, "shell_info", ORS,PSH,SCALZET
C
      IF (UCC(ORS(1:2)) == 'S ') THEN
        IPSHA = 0
        NSHL = NSHL + 1
      ELSE IF (UCC(ORS(1:2)) == 'P ') THEN
        IPSHA = 1
        NSHL = NSHL + 1
      ELSE IF (UCC(ORS(1:2)) == 'D ') THEN
        IPSHA = 2
        NSHL = NSHL + 1
      ELSE IF (UCC(ORS(1:2)) == 'F ') THEN
        IPSHA = 3
        NSHL = NSHL + 1
      ELSE IF (UCC(ORS(1:2)) == 'SP') THEN
        IPSHA = 0
        IPSHB = 1
        NSHL = NSHL + 2
        TWOSH = .TRUE.
      ELSE
        GOTO 9991
      END IF
C
      IF (TWOSH) THEN
        IF (NSHL > MXSH) GOTO 9992 
C
C     ___ Set shell pointers ___
C
        PSHELL(NSHL-1,1) = ATOM
        PSHELL(NSHL,1) = ATOM
        PSHELL(NSHL-1,2) = IPSHA
        PSHELL(NSHL,2) = IPSHB
        PSHELL(NSHL-1,3) = PSH
        PSHELL(NSHL,3) = PSH
C
        ICA = 0
        DO I=IPSHA,0,-1
          DO J=IPSHA-I,0,-1
            ICA = ICA + 1
          END DO
        END DO
C
        ICB = 0
        DO I=IPSHB,0,-1
          DO AY=IPSHB-I,0,-1
            ICB = ICB + 1
          END DO
        END DO
C
        NSTO = NSTO + ICA + ICB 
C
C      ___ Loop over primitive gaussian functions ___
C
        NGTF1 = NGTF
        NGTF2 = NGTF + PSH
        IF ((NGTF+2*PSH) > MXGTF) GOTO 9993
        DO IG=1,PSH
          NGTF1 = NGTF1 + 1
          NGTF2 = NGTF2 + 1
          ILINE = ILINE + 1
          READ(BASIS,*,ERR=9990) ZET(NGTF1),GCC(NGTF1),GCC(NGTF2)
C          print *, "base1", ZET(NGTF1),GCC(NGTF1),GCC(NGTF2)
          ZET(NGTF1) = ZET(NGTF1)*SCALZET**2
          ZET(NGTF2) = ZET(NGTF1)
        END DO
        NGTF = NGTF + 2*PSH
C
      ELSE
C
        IF (NSHL > MXSH) GOTO 9992
C
C     ___ Set shell pointers ___
C
        PSHELL(NSHL,1) = ATOM
        PSHELL(NSHL,2) = IPSHA
        PSHELL(NSHL,3) = PSH
C        print *, "PSHELL", PSHELL(NSHL,1), PSHELL(NSHL,2), 
C     $              PSHELL(NSHL,3)
C
        ICA = 0
        DO I=IPSHA,0,-1
          DO J=IPSHA-I,0,-1
            ICA = ICA + 1
          END DO
        END DO
C
        NSTO = NSTO + ICA
C
C      ___ Loop over primitive gaussian functions ___
C
        DO IG=1,PSH
          NGTF = NGTF + 1
          IF (NGTF > MXGTF) GOTO 9993
          ILINE = ILINE + 1
          READ(BASIS,*,ERR=9990) ZET(NGTF),GCC(NGTF)
C          print *, "base2", ZET(NGTF),GCC(NGTF)
          ZET(NGTF) = ZET(NGTF)*SCALZET**2
        END DO
C
      END IF
C
      RETURN
C
 9990 CONTINUE
      WRITE(OUT,1000) ILINE
      STOP
C
 9991 CONTINUE
      WRITE(OUT,2000)
      STOP
C
 9992 CONTINUE
      WRITE(OUT,3000)
      STOP
C
 9993 CONTINUE
      WRITE(OUT,4000)
      STOP
C
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Format statements ___
C
 1000 FORMAT(/,T2,'READSHGF: ERROR READING BASIS SET INPUT LINE',I5)
C
 2000 FORMAT(/,T2,'MAXIMUM L QUANTUM NUMBER FOR ORBITAL BASIS CANNOT ',
     $            'BE LARGER THAN 3')
C
 3000 FORMAT(/,T2,'NUMBER OF SHELLS EXCEEDED.',
     $       /,T2,'INCREASE THE PARAMETER MXSH IN parameter.h AND',
     $       /,T2,'RECOMPILE AGAIN')
C
 4000 FORMAT(/,T2,'NUMBER OF GTFs EXCEEDED.',
     $       /,T2,'INCREASE THE PARAMETER MXGTF IN parameter.h AND',
     $       /,T2,'RECOMPILE AGAIN')
C
C     O________________________________________________________________O
C
      END SUBROUTINE READSHGF

