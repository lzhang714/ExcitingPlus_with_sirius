
-------------------------------------------------------------------------------------------------
module environment on NERSC Cori (Cray XC40): 

zhang714@cori07:~/excitingplus_with_sirius> module list
Currently Loaded Modulefiles:
  1) modules/3.2.10.6                                 12) job/2.2.3-6.0.7.0_44.1__g6c4e934.ari
  2) nsg/1.2.0                                        13) dvs/2.7_2.2.117-6.0.7.1_9.2__gf817677
  3) intel/18.0.1.163                                 14) alps/6.6.43-6.0.7.0_26.4__ga796da3.ari
  4) craype-network-aries                             15) rca/2.2.18-6.0.7.0_33.3__g2aa4f39.ari
  5) craype/2.5.15                                    16) atp/2.1.3
  6) udreg/2.3.2-6.0.7.0_33.18__g5196236.ari          17) PrgEnv-intel/6.0.4
  7) ugni/6.0.14.0-6.0.7.0_23.1__gea11d3d.ari         18) craype-haswell
  8) pmi/5.0.14                                       19) cray-mpich/7.7.3
  9) dmapp/7.1.1-6.0.7.0_34.3__g5a674e0.ari           20) altd/2.0
 10) gni-headers/5.0.12.0-6.0.7.0_24.1__g3b1768f.ari  21) cray-hdf5/1.10.2.0
 11) xpmem/2.2.15-6.0.7.1_5.11__g7549d06.ari


--------------------------------------------------------------------------------------------------
codes involving sirius: 

./src/modsirius.f90            // newly created
./src/modmain.f90              // added some variables at the end
./src/addons/mod_mpi_grid.f90  // added some variables and a subroutine called "gatherir"
./src/energy.f90
./src/genapwfr.f90
./src/gencfun.f90
./src/gengvec.f90
./src/genlofr.f90   
./src/genveffig.f90
./src/gndstate.f90             // greatly simplified 
./src/hmlrad_sirius.f90        // newly created, it is same as the hmlint.f90 in EXCITING.
./src/init0.f90
./src/init1.f90
./src/olprad.f90
./src/poteff.f90


