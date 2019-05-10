non-magnetic bulk NiO

Tried to use Sirius to generate G-vectors, but got our-of-bound error before finishing init0, check elk.out and slurm-*.out. 
It is out of bound by 1, which looks like mis-match of array indexing between C++ and Fortran. 
But currently have no clue to fix it. 
