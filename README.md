# scf-cpp
SCF program written for a M.Sc. course in quantum chemistry.

Usage: ```./scf -i <path-to-input> [-rhf -uhf] [-geom -bas]```

Giving an input file is required. The other options are optional. The default is a RHF single-point calculation.
Bugs can arise when trying to run a RHF calculation with unpaired electrons since this is not handled in the input parser.

scf-dscan is a version of the normal SCF program used to automatically scan electron densities. It will scan through each bond if there are multiple atoms and through the atom if there is only one.

Usage: ```./scf-dscan -i <path-to-input>```

It calculates the electron densities by running a RHF single-point. 

All output is printed to stdout, so it should be piped to a file if you want to keep the information.
