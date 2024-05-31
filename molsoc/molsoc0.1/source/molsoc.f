C     O_______________--------------oOOOo--------------________________O
                               PROGRAM MolSOC       
C      ****************************************************************
C       ****************************     ****************************
C       --------------------------------------------------------------
C      ***     Calculate Spin-Orbit Coupling between preoptimized   ***
C     ***                         states                             ***
C      ---------------------------------------------------------------- 
C                                  AUTHOR:
C                           Sandro Giuseppe Chiodo
C     O~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~O
C
C     ___ Open permanent files ___
C
      CALL FILEIO('OPEN PERMANENT FILES')
C
C     ___ Write copyright ___
C
      CALL COPYR
C
C     ___ Get date and time ___
C
      CALL GETDATE
C
C     ___ Initialization ___
C
      CALL INITIALIZE
C
C     ___ Set orbital pointers ___
C
      CALL POINTER
C
C     ___ Read input data ___
C
      CALL READINP
C
C     ___ Write input header ___
C
      CALL WRTHEAD
C
C     ___ Open MO permanent files ___
C
      CALL FILEIO('OPEN MOS PERMANENT FILES')
C
C     ___ Read Gaussian basis set ___
C
      CALL READBAS
C
C     ___ Generate spherical orbital pointer ___
C
      CALL SHELLOP
C
C     ___ Calculate number of electrons ___
C
      CALL FINNEL
C
C     ___ Check if charge and multiplicities are consistent ___
C
C      CALL CHKMULT
C
C     ___ Calculate the square of the distance matrix ___
C
      CALL SDIST
C
C     ___ Initialization of mathematical constants ___
C
      CALL MATHCON
      CALL TRANSC
      CALL CRFTABLE
C
C     ___ Generate orbital pointer ___
C
      CALL GENOP
C
C     ___ Write outuput informations ___
C
      CALL WRITOUT
C
C     ___ Normalize Gaussian-type functions ___
C
      CALL NORGTF
C
C     ___ Normalize slater type orbitals ___
C
      CALL NORSTO
C
C     ___ Generate orbital symmetry label ___
C
      CALL ORBSYM
C
C     ___ Generate moment symmetry label ___
C
      CALL MOMESYM
C
C     ___ Calculate one-electron SOC integrals ___
C
      CALL SOCONE
C
C     ___ Calculate overlap integrals ___
C
      CALL BLDOVER
C
C     ___ Calculate moments integrals ___
C
      CALL MOMENT
C
C     ___ Calculate spin-orbit coupling ___
C
C      CALL BLDSOC
C
C     ___ Get date and time ___
C
      CALL GETDATE
C
C     ___ Close temporary files ___
C
      CALL FILEIO('CLOSE ALL FILES')
C
C     O_______________--------------oOOOo--------------________________O
C      ***                 End of PROGRAM MolSOC                    ***
C     ******************************************************************
C
      END PROGRAM MolSOC
