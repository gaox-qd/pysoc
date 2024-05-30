# PySOC - Calculation of Spin-Orbit Coupling

**Now updated to Python3.0!**

PySOC is a program, library, and interface to calculate spin-orbit coupling (SOC) between singlet ground states/excited states and triplet excited states.

PySOC currently supports excited state calculations with Gaussian and DFTB+.

## Basic features

 - Evaluation of spin-orbit coupling elements between singlet and triplet states  
 - Python scripts + FORTRAN  
 - Interfaced to third-party quantum chemistry packages, such as Gaussian 09 and DFTB+  
 - Based on Casida’s wave functions in LR-TDDFT, TDA, TDDFTB  
 - Breit-Pauli spin-orbit Hamiltonian with effective charge approximation  

## Dependencies

```console
$ pip install cclib tabulate periodictable scipy
```

## Pre-calculation

Before spin-orbit coupling can be calculated, it is first necessary to perform an excited states calculation with either the Gaussian or DFTB+ calculation programs.
The next section will briefly cover the necessary options to perform SOC calculations with each,
but the reader is expected to be familiar with the basic operation of the QM program of choice.

> [!NOTE]
> Shells higher than the f shell are not currently supported by PySOC. The user should ensure that the basis set used in the QM program does not contain g (or higher) shells.

### Gaussian

Time-dependant DFT excited states calculations both with and without the Tamm–Dancoff approximation are supported for Gaussian versions 09 and 16.
Older versions of Gaussian and alternative calculation methods (CIS etc.) may additionally be supported, but have not been tested.

PySOC requires both the log file (.log/.out) and read/write binary file (.rwf) from the calculation as input.
The .rwf file is normally deleted automatically by Gaussian after the calculation finishes. You can prevent this behaviour by adding the following line to the Gaussian input file:

```
%Rwf="FILE.rwf"
```
Where FILE can be changed as appropriate.

Additionally, PySOC requires that the basis set be printed in the calculation output. This can be requested with the following options:

```
6D 10F GFInput
```

A full example Gaussian input file for formaldehyde would look like this:

```
%mem=1GB
%chk=formaldehyde.chk
%rwf=formaldehyde.rwf
# TD(50-50,nstates=5) wB97XD/TZVP 6D 10F GFInput

test

0 1

C         -0.131829      -0.000001      -0.000286
O          1.065288       0.000001       0.000090
H         -0.718439       0.939705       0.000097
H         -0.718441      -0.939705       0.000136


```

### DFTB+

> [!NOTE]
> This section is under revision.

> [!NOTE]  
> PySOC requires Gaussian type orbitals (GTOs) to perform the SOC calculation, yet DFTB+ uses Slater-type orbitals (STOs).
> Thus the parameter set used for the DFTB+ calculation must be fitted to GTOs for use with PySOC.
> PySOC contains a fitted set for the mio-1-1 parameter set which will be used by default, alternative fitted sets can be specified with the `--fitted_basis` option.

An example dftb_in.hsd file is given below:

```
Geometry = GenFormat {
  <<< "ch2o.gen"
}

Driver = {}

Hamiltonian = DFTB {
  SCC = Yes
  SCCTolerance = 1e-10  # Very tight for test purposes only
  MaxAngularMomentum = {
    H  = "s"
    C  = "p"
    O  = "p"
  }
  SlaterKosterFiles = Type2FileNames {
    Prefix = "/fsnfs/users/xinggao/work/gsh/thiothymine/gtsh/test_python/tddftb/sk/mio-1-1/"
    Separator = "-"
    Suffix = ".skf"
  }
  LinearResponse { 
    NrOfExcitations = 10
    StateOfInterest = 0
    Symmetry = both
    HubbardDerivatives {
      H  = 0.347100   0.491900 
      C  = 0.341975   0.387425 
      O  = 0.467490   0.523300
    }
    WriteTransitions = Yes
    WriteTransitionDipole = Yes
    WriteXplusY = Yes
  }
}

Options {
  WriteEigenvectors = Yes
  WriteHS = No
}

ParserOptions {
  ParserVersion = 4
}
```

After performing the DFTB+ calculation, set WriteHS = Yes and run the calculation a second time.

## Usage

Once the QC/TB calculation has completed, SOC can be calculated using the `pysoc` command.

`pysoc` has only one mandatory argument: the calculation output file. For Gaussian, this is the main calculation .log file, while for DFTB+ this is the .xyz geometry file.

For example, to calculate SOC for a file in the current directory with the name ch2o.log:

```console
$ pysoc formaldehyde.log
```

The rwf file, which is also required, will be guessed based on the name of the given log file. The rwf file can also be specified
explicitly:

```console
$ pysoc formaldehyde.log --rwf-file gaussian.rwf
```

To calculate SOC from a completed DFTB+ calculation:

```console
$ pysoc formaldehyde.xyz
```

Once the calculation is complete, the SOC table will be printed:

```console
$ pysoc ch2o.log
 Singlet    Triplet     RSS (cm-1)    +1 (cm-1)    0 (cm-1)    -1 (cm-1)
---------  ---------  ------------  -----------  ----------  -----------
   S0         T1           60.7419      42.9510      0.0077      42.9510
   S0         T2            0.0194       0.0137      0.0000       0.0137
   S0         T3           10.6234       0.0133     10.6234       0.0133
   S0         T4           59.8873      42.3467      0.0012      42.3467
   S0         T5           11.7332       0.0013     11.7332       0.0013
   S1         T1            0.0008       0.0006      0.0000       0.0006
   S1         T2           44.3144      31.3350      0.0224      31.3350
   S1         T3            8.6280       6.1009      0.0002       6.1009
   S1         T4           50.7211       0.0151     50.7211       0.0151
   S1         T5            5.5317       3.9115      0.0001       3.9115
   S2         T1            7.3854       5.2223      0.0001       5.2223
   S2         T2            0.2532       0.0008      0.2531       0.0008
   S2         T3            0.0003       0.0002      0.0000       0.0002
   S2         T4            0.2432       0.1720      0.0013       0.1720
 ...
```

Here, each line indicates SOC between one singlet and one triplet state, which are given by the ‘Singlet’ and ‘Triplet’ columns respectively.
The last three columns each contain the calculated spin-orbit coupling for the triplet sub-state with quantum number +1, 0, or -1.
The root-sum-square for the three triplet sub-states is presented in the RSS column.

Alternatively, a comma-separated values (CSV) format can be requested with the ‘-c’ option, which can be more easily imported into a spreadsheet for analysis:

```console
$ pysoc ch2o.log -c
Singlet,Triplet,RSS (cm-1),+1 (cm-1),0 (cm-1),-1 (cm-1)
S0,T1,60.7419154885397,42.95102,0.00769,42.95102
S0,T2,0.019431296920174937,0.01374,1e-05,0.01374
S0,T3,10.623396701112124,0.01332,10.62338,0.01332
S0,T4,59.887319900989056,42.34673,0.00124,42.34673
S0,T5,11.733160146260683,0.00131,11.73316,0.00131
...
```

To write to a file instead of to the screen, use standard Linux file redirection with the ‘>’ or ‘>>’ characters (to overwrite or append to an existing file respectively):

```console
$ pysoc ch2o.log -c > SOC.csv
```

### Intermediate Files
PySOC generates a number of intermediate files during operation. By default, these are deleted at the end of the calculation and only the final SOC values are displayed.
To override this behaviour and keep these intermediate files, use the `-o` (output) option with a path to a directory where the intermediate files should be stored.
One may use `-o ./` to write the intermediate files to the current directory. The `-o` option has no effect on the final SOC values:

```console
$ pysoc ch2o.log -o ./
 Singlet    Triplet     RSS (cm-1)    +1 (cm-1)    0 (cm-1)    -1 (cm-1)
---------  ---------  ------------  -----------  ----------  -----------
   S0         T1           60.7419      42.9510      0.0077      42.9510
   S0         T2            0.0194       0.0137      0.0000       0.0137
   S0         T3           10.6234       0.0133     10.6234       0.0133
   S0         T4           59.8873      42.3467      0.0012      42.3467
   S0         T5           11.7332       0.0013     11.7332       0.0013
...
```

## Additional documentation

1. code: pysoc.tar.gz; short tutorial: [doc/pysoc.pdf](doc/pysoc.pdf)  
1. Even more better tutorial written in Chinese by sobereva Lu Tian:http://bbs.keinsci.com/thread-9442-1-1.html  
1. For Gaussian 16 by ggdh:http://bbs.keinsci.com/thread-19813-1-1.html  

## Reference and Citation

Evaluation of Spin-Orbit Couplings with Linear-Response Time-Dependent Density Functional Methods  
Xing Gao, Shuming Bai, Daniele Fazzi, Thomas Niehaus, Mario Barbatti, and Walter Thiel  
J. Chem. Theory Comput., 2017, 13 (2), pp 515–524

DOI: 10.1021/acs.jctc.6b00915  

## Authorship

PySOC was originally written by Xing Gao et al. for python 2.x (J. Chem. Theory Comput. 13 (2017) 515–524).

[MolSOC](molsoc/), which is used for the calculation of atomic integrals, was originally written by Sandro Giuseppe Chiodo et al. (Computer Physics Communications 185 (2014) 676–683).

PySOC was re-written for Python 3.x by Oliver S. Lee.