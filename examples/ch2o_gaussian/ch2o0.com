%mem=1GB
%chk=gaussian.chk
%rwf=gaussian.rwf
# td(50-50,nstates=5) wB97XD/TZVP 6D 10F nosymm GFInput
                 
test
                 
0 1
C         -0.131829      -0.000001      -0.000286
O          1.065288       0.000001       0.000090
H         -0.718439       0.939705       0.000097
H         -0.718441      -0.939705       0.000136

