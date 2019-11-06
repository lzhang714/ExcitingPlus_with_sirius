! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! 
! DESCRIPTION:
!   Computes the self-consistent Kohn-Sham ground-state. General information is written to {\tt INFO.OUT}. 
!   1st- and 2nd- variational eigenvalues, eigenvectors and occupancies are written to unformatted files: 
!   {\tt EVALFV.OUT}, {\tt EVALSV.OUT}, {\tt EVECFV.OUT}, {\tt EVECSV.OUT} and {\tt OCCSV.OUT}. 
! 
! 2018 Dec. L.Zhang @ UFL, simplified and modified to interface to the SIRIUS library
! 

subroutine gndstate_900

  use modmain
  use mod_wannier
  use mod_sic
  use mod_seceqn
  use mod_addons

  implicit none

  ! local variables
  logical exist
  integer ik,is,ia,ias,ist,i,ikloc,idm
  integer n,nwork
  real(8) dv,etp,de,timetot
  
  ! allocatable arrays
  real(8), allocatable :: v(:)
  real(8), allocatable :: work(:)
  ! real(8), allocatable :: evalfv(:,:,:)        ! renamed to evalfvloc, declaration has been moved to addons/mod_addons.f90
  character(100) :: label
  
  call timer_start(t_init,reset=.true.)  
  call init0                                 ! sirius calls inside
  call init1                                 ! sirius calls inside
  call timer_stop(t_init)
   
  if (.not.mpi_grid_in()) return
  
  wproc = mpi_grid_root()  
  
  ! allocate arrays for eigenvalues/vectors
  ! "tsveqn" assigned .true. in mod_addons.f90, we never changed.
  if (tsveqn) then
    allocate(evalfvloc(nstfv,nspnfv,nkptloc))          ! renamed to evalfvloc
    allocate(evecfvloc(nmatmax,nstfv,nspnfv,nkptloc))
    allocate(evecsvloc(nstsv,nstsv,nkptloc))
  else
    allocate(evecfdloc(nspinor*nmatmax,nstsv,nkptloc))
  endif

  ! write info and open files
  if (wproc) then
    ! write the real and reciprocal lattice vectors to file
    call writelat
    ! write interatomic distances to file
    call writeiad(.false.)
    ! write symmetry matrices to file
    call writesym
    ! output the k-point set to file
    call writekpts
    ! write lattice vectors and atomic positions to file
    call writegeom(.false.)
    ! open INFO.OUT file
    open(ioinfo,file='INFO'//trim(filext),action='WRITE',form='FORMATTED')
    ! open TOTENERGY.OUT
    open(61,file='TOTENERGY'//trim(filext),action='WRITE',form='FORMATTED')
    ! open FERMIDOS.OUT
    open(62,file='FERMIDOS'//trim(filext),action='WRITE',form='FORMATTED')
    ! open MOMENT.OUT if required
    if (spinpol) open(63,file='MOMENT'//trim(filext),action='WRITE', form='FORMATTED')
    ! open RMSDVEFF.OUT
    open(65,file='RMSDVEFF'//trim(filext),action='WRITE',form='FORMATTED')
    ! open DTOTENERGY.OUT
    open(66,file='DTOTENERGY'//trim(filext),action='WRITE',form='FORMATTED')
    ! write out general information to INFO.OUT
    call writeinfo(ioinfo)
    write(ioinfo,*)
    write(ioinfo,'("initialization time : ",F10.2," sec.")')timer_get_value(t_init)
    call writenn
    call writegshells
  endif
  
  ! size of mixing vector
  n = lmmaxvr * nrmtmax * natmtot + ngrtot
  if (spinpol) n = n*(1+ndmag)
  ! allocate mixing arrays
  allocate(v(n))
  ! determine the size of the mixer work array
  nwork=-1
  call mixerifc(mixtype,n,v,dv,nwork,v)
  allocate(work(nwork))

  
  ! -----------------------------------------------------------------! 
  ! initialize or read the charge density and potentials from file   ! 
  ! -----------------------------------------------------------------! 
  
  if (wproc) write(ioinfo,*)
  
  ! ---------------------- task selection 
  if ((task.eq.1).or.(sic.and.isclsic.gt.1)) then

      call readstate
      if (wproc) write(ioinfo,'(" Potential has been read in from STATE.OUT ")')
      if (autolinengy) call readfermi

  else if (task.eq.200) then

      call phveff
      if (wproc) write(ioinfo,'(" Supercell potential constructed from STATE.OUT ")')

  else

      call timer_start(t_rhoinit,reset=.true.)
      call allatoms 
      do ias=1,natmtot
        is=ias2is(ias)
        do ist=1,spnst(is)
          evalcr(ist,ias)=speval(ist,is)
        enddo
      enddo
      !
      ! check that we can run the whole ground state; used to debug the EP/SIRIUS interface;
      ! 
      if (use_sirius_library.and.sirius_run_full_scf) then
        ! set the free atomic desnity for each atom type
        ! we can do it only here after we called allatoms() to solve the isolated atom
        !   and generate free-atom density sprho
        do is=1,nspecies
          label = adjustl(spfname(is))
          call sirius_add_atom_type_radial_function(sctx, string(trim(label)), string("ae_rho"),&
                                                    &sprho(1, is), spnr(is))
        enddo
        write(*,*)' SIRIUS full scf started! '
        call sirius_find_ground_state(gs_handler, potential_tol=epspot, energy_tol=epsengy, niter=maxscl,&
                                     &save_state=bool(.false.))
        write(*,*)' SIRIUS full scf done! '
      else
        ! origianl EP rhoinit should be "OR"ed with sirius full scf run
        call rhoinit  
      endif
      
      call poteff      ! sirius calls 
      !call genveffig   ! EXCITING
      call timer_stop(t_rhoinit)

      if (wproc) then
        write(ioinfo,'("Density and potential initialised from atomic data")')
        write(ioinfo,'(" done in : ",F10.2," sec.")')timer_get_value(t_rhoinit)
        write(ioinfo,*)
        do is=1,nspecies
          write(ioinfo,'("Species : ",I4," (",A,")")') is,trim(spsymb(is))
          do ist=1,spnst(is)
            write(ioinfo,'("  n : ",I1,"   l : ",I1,"   k : ",I1,"   energy : ",G18.10)')&
            &spn(ist,is),spl(ist,is),spk(ist,is),speval(ist,is)
          enddo
          write(ioinfo,*)
        enddo
      endif
  
  end if
  ! ---------------------- task selection 
  
  if (wproc) call flushifc(ioinfo)
  
  tstop=.false.   
  tlast=.false.   
  etp=0.d0


!moved the chunk below to the task selection above,
!  !
!  ! check that we can run the whole ground state; used to debug the EP/SIRIUS interface
!  ! 
!  if (use_sirius_library.and.sirius_run_full_scf) then
!    ! set the free atomic desnity for each atom type
!    ! we can do it only here after we called allatoms() to solve the isolated atom
!    !   and generate free-atom density sprho
!    do is=1,nspecies
!      label = adjustl(spfname(is))
!      call sirius_add_atom_type_radial_function(sctx, string(trim(label)), string("ae_rho"),&
!                                                &sprho(1, is), spnr(is))
!    enddo
!    call sirius_find_ground_state(gs_handler, bool(.false.))
!    write(*,*)'SIRIUS scf done!'
!  endif
!
!the "sirius_find_ground_state" will do the following, 
!
!/* @fortran begin function void sirius_find_ground_state        Find the ground state
!   @fortran argument in required void* gs_handler               Handler of the ground state
!   @fortran argument in optional bool  save__                   boolean variable indicating if we want to save the ground state
!   @fortran end 
!   */
!   
!void sirius_find_ground_state(void* const* gs_handler__, bool const *save__)
!{
!    GET_GS(gs_handler__)
!    auto& ctx = gs.ctx();
!    auto& inp = ctx.parameters_input();
!    gs.initial_state();
!
!    bool save{false};
!    if (save__ != nullptr) {
!        save = *save__;
!    }
!
!    auto result = gs.find(inp.potential_tol_, inp.energy_tol_, inp.num_dft_iter_, save);
!}
  
  ! ================================================================== ! 
  !                               SCF loop                             ! 
  ! ================================================================== ! 

  if (wproc) then
    write(ioinfo,*)
    write(ioinfo,'("+------------------------------+")')
    write(ioinfo,'("| Self-consistent loop started |")')
    write(ioinfo,'("+------------------------------+")')
  endif

  iscl=0
  do iscl=1,maxscl

    ! reset timers and check max-iter
    call timer_reset(0)
    call timer_start(t_iter_tot)
    if (wproc) then
      write(ioinfo,*)
      write(ioinfo,'("+-------------------------+")')
      write(ioinfo,'("| Iteration number : ",I4," |")') iscl
      write(ioinfo,'("+-------------------------+")')
      call flushifc(ioinfo)
    endif
    if (iscl.ge.maxscl) then
      if (wproc) then
        write(ioinfo,*)
        write(ioinfo,'("Reached self-consistent loops maximum")')
      endif
      tlast=.true.
    end if
    if (wproc) call flushifc(ioinfo)

    ! Fourier transform effective potential to G-space
    !
    ! will pass pw coeff to sirius if (use_sirius_library.and.pass_veffig_to_sirius)
    !
    call genveffig                                        ! sirius calls inside

    ! set the new spherical potential in SIRIUS
    if (use_sirius_library.and.update_atomic_pot) then
#ifdef _SIRIUS_
      call sirius_update_atomic_potential(gs_handler)
#else
      stop sirius_error
#endif
    endif
    
    
    ! linengy + genapwfr + genlofr + getufr + genufrp calls were originally grouped in addons/genradf.f90. 
    ! it's broken down here just for clarity. 
        
    ! generate the core wavefunctions and densities
    call gencore                              
    ! find the new linearisation energies
    call linengy                              
    ! write out the linearisation energies
    if (wproc) call writelinen
    ! generate the APW radial functions
    call genapwfr                                 ! sirius call inside
    ! generate the local-orbital radial functions
    call genlofr                                  ! sirius call inside
    ! collect radial mt functions, check mod_addons
    call getufr
    ! compute the product of radial functions, check mod_addons
    call genufrp
    ! compute the overlap radial integrals
    call olprad                                   ! sirius call inside
    ! compute the Hamiltonian radial integrals
    ! hmlrad_sirius.f90 is newly created, following the hmlint.f90 of Exciting
    if (use_sirius_library.and.pass_hmlrad_to_sirius) then
        call hmlrad_sirius                        ! sirius calls inside 
    else 
        call hmlrad
    endif    
    ! generate effective magntic field integrals for full diagonalization 
    call genbeff
    ! generate muffin-tin effective magnetic fields and s.o. coupling functions 
    call genbeffmt
    
    evalsv(:,:)=0.d0
    if (sic) sic_evalsum=0.d0

    ! ========================================= loop k-points
    if (use_sirius_library.and.use_sirius_eigen_states) then

#ifdef _SIRIUS_
        !call sirius_update_atomic_potential
        !call sirius_generate_radial_integrals
        !call sirius_generate_potential_pw_coefs
        call sirius_find_eigen_states(gs_handler, ks_handler, logical(.false.,kind=c_bool))     ! the actual diagonalisation
        do ik = 1, nkptloc
             if (ndmag.eq.0.or.ndmag.eq.3) then
               call sirius_get_band_energies(ks_handler, ik, 0, evalsv(1, ik))   ! get eigenvalues from SIRIUS
             else
               call sirius_get_band_energies(ks_handler, ik, 0, evalsv(1, ik))
               call sirius_get_band_energies(ks_handler, ik, 1, evalsv(nstfv+1, ik))
             endif
             !TODO: get eigen vectors from SIRIUS. do it here? or at the end of scf iteration (after converge, before writing STATE.OUT)?
        enddo
#else
        stop sirius_error
#endif
        !TODO: add wannier downfolding here, like that in seceqn.f90
        !if (wannier) then
        !  call wan_gencsv_aux(ikloc,evecfv=evecfv,evecsv=evecsv)
        !  if (wann_add_poco) then
        !    call wann_seceqn(ikloc,evecsv)
        !    call wan_gencsv_aux(ikloc,evecfv=evecfv,evecsv=evecsv)
        !  endif
        !endif
        call mpi_grid_reduce(evalsv(1,1),nstsv*nkpt,dims=(/dim_k/),all=.true.)    
        
    else 
        
        do ikloc=1,nkptloc
          call seceqn( ikloc, evalfvloc(1,1,ikloc), evecfvloc(1,1,1,ikloc), evecsvloc(1,1,ikloc) )        ! renamed to evalfvloc
        enddo
        call mpi_grid_reduce(evalsv(1,1),nstsv*nkpt,dims=(/dim_k/),all=.true.)
        
    endif
    ! ========================================= loop k-points
    
    ! find occsv(nstsv,nkpt) and E-Fermi, then mpi-bcast occsv
    if (wproc) then
      call occupy
      if (autoswidth) then
        write(ioinfo,*)
        write(ioinfo,'("New smearing width : ",G18.10)') swidth
      end if    
      call writeeval        ! write eigenvalues and occupation numbers to file EIGVAL.OUT
      call writefermi       ! write Fermi energy to file EFERMI.OUT
      call flushifc(ioinfo)
    endif
    call mpi_grid_bcast(swidth)
    call mpi_grid_bcast(occsv(1,1),nstsv*nkpt)

    ! wannier occupancy and energy matrices
    if (wannier) call wann_ene_occ
    
    ! density and magnetisation
    if (.false.) then             ! if (use_sirius_library.and.use_sirius_density) then
!#ifdef _SIRIUS_
!        do ik = 1, nkpt
!            if (ndmag.eq.0.or.ndmag.eq.3) then
!              call sirius_set_band_occupancies(ks_handler, ik, 0, occsv(1, ik))   ! pass occsv to sirius
!            else
!              call sirius_set_band_occupancies(ks_handler, ik, 0, occsv(1, ik))
!              call sirius_set_band_occupancies(ks_handler, ik, 1, occsv(nstfv+1, ik))
!            endif
!        enddo
!        ! let sirius generate charge density
!        call sirius_generate_density(gs_handler)
!        ! TODO: get the charge density from SIRIUS. do it here or at the end of scf iteration (after converge, before writing STATE.OUT)?
!        call gatherir(rhoir(1))
!        do idm = 1, ndmag
!          call gatherir(magir(1, idm))
!        enddo
!#else
!        stop sirius_error
!#endif
    else
          call rhomag
    endif
    
    if (sic) call mpi_grid_reduce(sic_evalsum,dims=(/dim_k/))

    ! LDA+U commented out for the moment. 
    ! LDA+U
    !if (iscl.gt.1.and.ldapu.ne.0) then  
    !  ! generate the LDA+U density matrix
    !  call gendmatrsh
    !  ! generate the LDA+U potential matrix
    !  call genvmatlu
    !  ! write the LDA+U matrices to file
    !  if (wproc) call writeldapu
    !  ! calculate and write tensor moments to file
    !  if (tmomlu.and.wproc) then
    !    call tensmom(67)
    !    call flushifc(67)
    !  end if
    !end if
    
    ! compute the effective potential
    call poteff                                      ! sirius calls

    ! pack interstitial and muffin-tin effective potential and field into one array
    call mixpack(.true.,n,v)
    ! mix in the old potential and field with the new
    call mixerifc(mixtype,n,v,dv,nwork,work)
    ! unpack potential and field
    call mixpack(.false.,n,v)

    ! add the fixed spin moment effect field
    if (fixspin.ne.0) call fsmfield

    ! reduce the external magnetic fields if required
    if (reducebf.lt.1.d0) then
      bfieldc(:)=bfieldc(:)*reducebf
      bfcmt(:,:,:)=bfcmt(:,:,:)*reducebf
    end if

    ! compute the energy components
    call energy                                      ! sirius 
    call timer_stop(t_iter_tot)

    ! output energy components
    if (wproc) then
      call writeengy(ioinfo)
      write(ioinfo,*)
      write(ioinfo,'("Density of states at Fermi energy : ",G18.10)') fermidos
      write(ioinfo,'(" (states/Hartree/unit cell)")')
      write(ioinfo,'("Estimated band gap (eV) : ",G18.10)') bandgap(1)*ha2ev
      write(ioinfo,'(" from k-point ",I6," to k-point ",I6)') ikgap(1),ikgap(2)
      write(ioinfo,'("Estimated direct band gap (eV) : ",G18.10)') bandgap(2)*ha2ev
      write(ioinfo,'(" at k-point ",I6)') ikgap(3)
      call flushifc(ioinfo)
      ! write total energy to TOTENERGY.OUT and flush
      write(61,'(G22.12)') engytot
      call flushifc(61)
      ! write DOS at Fermi energy to FERMIDOS.OUT and flush
      write(62,'(G18.10)') fermidos
      call flushifc(62)
      ! output charges and moments
      call writechg(60)
      ! write total moment to MOMENT.OUT and flush
      if (spinpol) then
        write(63,'(3G18.10)') momtot(1:ndmag)
        call flushifc(63)
      end if
      ! output effective fields for fixed spin moment calculations
      if (fixspin.ne.0) call writefsm(60)
      ! check for WRITE file
      inquire(file='WRITE',exist=exist)
      if (exist) then
        write(ioinfo,*)
        write(ioinfo,'("WRITE file exists - writing STATE.OUT")')
        call writestate
        open(50,file='WRITE')
        close(50,status='DELETE')
      end if
      ! write STATE.OUT file if required
      if (nwrite.ge.1) then
        if (mod(iscl,nwrite).eq.0) then
          call writestate
          write(ioinfo,*)
          write(ioinfo,'("Wrote STATE.OUT")')
        end if
      end if
    endif ! wproc

    ! exit self-consistent loop if last iteration is complete
    call mpi_grid_bcast(tlast)
    if (tlast) goto 20

    ! check for convergence
    if (wproc) then
      if (iscl.ge.2) then
        write(ioinfo,*)
        write(ioinfo,'("RMS change in effective potential (target) : ",G18.10," (",G18.10,")")') dv,epspot
        write(65,'(G18.10)') dv
        call flushifc(65)
        de=abs(engytot-etp)
        write(60,'("Absolute change in total energy (target)   : ",G18.10," (",G18.10,")")') de,epsengy
        write(66,'(G18.10)') de
        call flushifc(66)
        if ((dv.lt.epspot).and.(de.lt.epsengy)) then
          write(ioinfo,*)
          write(ioinfo,'("Convergence targets achieved")')
          tlast=.true.
        end if
      end if
      etp=engytot
      if (xctype(1).lt.0) then
        write(ioinfo,*)
        write(ioinfo,'("Magnitude of OEP residual : ",G18.10)') resoep
      end if
      ! check for STOP file
      inquire(file='STOP',exist=exist)
      if (exist) then
        write(ioinfo,*)
        write(ioinfo,'("STOP file exists - stopping self-consistent loop")')
        tstop=.true.
        tlast=.true.
        open(50,file='STOP')
        close(50,status='DELETE')
      end if
      ! output the current total CPU time
      timetot=timeinit+timemat+timefv+timesv+timerho+timepot+timefor
      write(ioinfo,*)
      write(ioinfo,'("Time (CPU seconds) : ",F12.2)') timetot
      write(ioinfo,*)
      write(ioinfo,'("iteration time (seconds)              : ",F12.2)' ) timer_get_value(t_iter_tot)
      write(ioinfo,'("  radial setup (linen, APW, radint)   : ",3F12.2)') timer_get_value(t_lin_en),&
        &timer_get_value(t_apw_rad),timer_get_value(t_hbo_rad)
      write(ioinfo,'("  total for secular equation          : ",F12.2)' ) timer_get_value(t_seceqn)
      write(ioinfo,'("  first-variational                   : ",F12.2)' ) timer_get_value(t_seceqnfv)
      write(ioinfo,'("  setup                               : ",F12.2)' ) timer_get_value(t_seceqnfv_setup)
      write(ioinfo,'("        setup H (total, MT, IT)             : ",3F12.2)')&
      &timer_get_value(t_seceqnfv_setup_h),&
      &timer_get_value(t_seceqnfv_setup_h_mt),&
      &timer_get_value(t_seceqnfv_setup_h_it)
      write(ioinfo,'("        setup O (total, MT, IT)             : ",3F12.2)')&
      &timer_get_value(t_seceqnfv_setup_o),&
      &timer_get_value(t_seceqnfv_setup_o_mt),&
      &timer_get_value(t_seceqnfv_setup_o_it)
      write(ioinfo,'("      diagonalization                       : ",F12.2)')&
      &timer_get_value(t_seceqnfv_diag)
      write(ioinfo,'("    second-variational                      : ",F12.2)')&
      &timer_get_value(t_seceqnsv)
      write(ioinfo,'("      setup (total, MT, IT)                 : ",3F12.2)')&
      &timer_get_value(t_seceqnsv_setup),&
      &timer_get_value(t_seceqnsv_setup_mt),&
      &timer_get_value(t_seceqnsv_setup_it)
      write(ioinfo,'("      diagonalization                       : ",F12.2)')&
      &timer_get_value(t_seceqnsv_diag)
      write(ioinfo,'("  total for charge and magnetization        : ",F12.2)')&
      &timer_get_value(t_rho_mag_tot)
      write(ioinfo,'("    k-point summation (total, MT, IT)       : ",3F12.2)')&
      &timer_get_value(t_rho_mag_sum),&
      &timer_get_value(t_rho_mag_mt),&
      &timer_get_value(t_rho_mag_it)
      if (texactrho) then
      write(ioinfo,'("      wave-function setup                   : ",F12.2)')&
        &timer_get_value(t_rho_wf)
      write(ioinfo,'("      convert to r-mesh                     : ",F12.2)')&
        &timer_get_value(t_rho_mag_conv)
      endif
      write(ioinfo,'("    symmetrization                          : ",F12.2)')&
      &timer_get_value(t_rho_mag_sym)
      write(ioinfo,'("  potential (total, XC, Ha)                 : ",3F12.2)')&
      &timer_get_value(t_pot),&
      &timer_get_value(t_pot_xc),&
      &timer_get_value(t_pot_ha)
      write(ioinfo,'("  density matrix setup                      : ",F12.2)')&
      &timer_get_value(t_dmat)
      if (sic) then
      write(ioinfo,'("  sic_genfvprj                              : ",F12.2)')&
        &timer_get_value(t_sic_genfvprj)
      write(ioinfo,'("  sic_hunif                                 : ",F12.2)')&
        &timer_get_value(t_sic_hunif)
      endif
    endif ! wproc
  
  enddo ! iscl


  ! TODO: get eigen vectors from SIRIUS. do it here? 
  ! TODO: get the charge density from SIRIUS. do it here?
  
  
  20 continue
  if (wproc) then
    write(ioinfo,*)
    write(ioinfo,'("+------------------------------+")')
    write(ioinfo,'("| Self-consistent loop stopped |")')
    write(ioinfo,'("+------------------------------+")')
    call flushifc(ioinfo)
    ! write density and potentials to file only if maxscl > 1
    if (maxscl.gt.1) then
      call writestate
      write(ioinfo,*)
      write(ioinfo,'("Wrote STATE.OUT")')
      call flushifc(ioinfo)
    end if
  endif
  
  ! ================================================================== ! 
  !                               SCF loop                             ! 
  ! ================================================================== ! 

  write(*,*)'the whole gndstate scf loop done!'  
  
  ! write eigenvalues/vectors and occupancies to file
  if (mpi_grid_side(dims=(/dim_k/))) then
    do i=0,mpi_grid_dim_size(dim_k)-1
      if (mpi_grid_dim_pos(dim_k).eq.i) then
        do ikloc=1,nkptloc
          ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
          if (tsveqn) then
            call putevalfv(ik,evalfvloc(1,1,ikloc))        ! renamed to evalfvloc
            call putevecfv(ik,evecfvloc(1,1,1,ikloc))
            call putevecsv(ik,evecsvloc(1,1,ikloc))
          else
            call putevecfd(ikloc,evecfdloc(1,1,ikloc))
          endif
          call putevalsv(ik,evalsv(1,ik))
          call putoccsv(ik,occsv(1,ik))
        end do
      end if
      call mpi_grid_barrier(dims=(/dim_k/))
    end do
  endif
  call mpi_grid_bcast(tstop)
  
  ! compute force
  if (tforce) then
    call force
    if (wproc) then
      call writeforce(ioinfo)
      call flushifc(ioinfo)
    endif
  end if

  ! output timing information
  if (wproc) then
    write(ioinfo,*)
    write(ioinfo,'("Timings (CPU seconds) :")')
    write(ioinfo,'(" initialisation",T40,": ",F12.2)') timeinit
    write(ioinfo,'(" Hamiltonian and overlap matrix set up",T40,": ",F12.2)') timemat
    write(ioinfo,'(" first-variational secular equation",T40,": ",F12.2)') timefv
    if (spinpol) then
      write(ioinfo,'(" second-variational calculation",T40,": ",F12.2)') timesv
    endif
    write(ioinfo,'(" charge density calculation",T40,": ",F12.2)') timerho
    write(ioinfo,'(" potential calculation",T40,": ",F12.2)') timepot
    if (tforce) then
      write(ioinfo,'(" force calculation",T40,": ",F12.2)') timefor
    endif
    timetot = timeinit + timemat + timefv + timesv + timerho + timepot + timefor
    write(ioinfo,'(" total",T40,": ",F12.2)') timetot
    write(ioinfo,*)
    write(ioinfo,'("+-----------------------------+")')
    write(ioinfo,'("| Elk version ",I1.1,".",I1.1,".",I3.3," stopped |")') version
    write(ioinfo,'("+-----------------------------+")')
    close(ioinfo)           ! close the INFO.OUT file
    close(61)               ! close the TOTENERGY.OUT file
    close(62)               ! close the FERMIDOS.OUT file
    if (spinpol) close(63)  ! close the MOMENT.OUT file
    close(65)               ! close the RMSDVEFF.OUT file
    close(66)               ! close the DTOTENERGY.OUT file
  endif

  deallocate(v,work)

  if (tsveqn) then
    deallocate(evalfvloc)    ! renamed to evalfvloc
    deallocate(evecfvloc)
    deallocate(evecsvloc)
  else
    deallocate(evecfdloc)
  endif

  ! finalize SIRIUS at the end
  if (use_sirius_library.and.use_sirius_init) then
#ifdef _SIRIUS_
    call sirius_free_handler(gs_handler)
    call sirius_free_handler(ks_handler)
    call sirius_free_handler(sctx)
    call sirius_finalize(logical(.false.,kind=c_bool))
#endif
  endif

  return

end subroutine 





    ! ===================================================================================== !
    ! some notes from reading EXCITING-SIRIUS and EP, just to help clarify code structures. !
    ! ===================================================================================== !
    
    
    
    ! ------------------------------------------------------------------------------------------------------- 
    ! ALL about Exciting-Sirius:
    ! 
    ! 1. hmlrad is NOT used anywhere in Exciting except:
    !    hmlrad and hmlint are both used in hartfolk.f90, and 
    !    hmlrad is also used in src_rdmft/rdmft.f90. 
    ! 
    ! 2. the variables haa/halo/hlolo are ONLY used in hmlrad.f90 in Exciting, they seems obseleted. 
    !
    ! 3. Exciting used hmlint many times:
    !    ./bandstr.f90:      Call hmlint
    !    ./effmass.f90:      Call hmlint
    !    ./fermisurf.f90:    call hmlint
    !    ./hartfock.f90:     Call hmlint
    !    ./scf_cycle.f90:    call hmlint
    !    ./seceqnsv2.f90:    call hmlint
    ! 
    ! 4. the variables haaij/haloij/hloloij are used in many places: 
    !    forcek.f90
    !    forcek2.f90 
    !    src_eigensystem/hmlint.f90
    !    src_eigensystem/hamiltonandoverlapsetup.f90
    !    src_iterative_solver/davidson.f90
    !    src_iterative_solver/hapwsapw.f90
    ! 
    ! 5. the variables haaintegrals/halointegrals/hlolointegrals are 
    !    local variables in hmlint.f90, they ONLY appear in hmlint.f90, 
    !    and they are used to calculate haaij/haloij/hloloij.
    ! 
    ! 6. apwi is matching coeffi. defined in matchapwi.f90
    !    EP has the original apwalm defined in matchf90
    ! 
    ! 7. overlap and Hamiltonian integral arrays have same and different naming: 
    !    EP:       overlap:     oalo, ololo; 
    !              Hamiltonian: haa,  hloa, hlolo; 
    !    Exciting: overlap:     oalo, ololo, h1aa, h1loa, h1lolo; 
    !              Hamiltonian: haaij,  hloaij, hloloij (haa/halo/hlolo abandent);
    !
    !    EP: haa/hloa/hlolo are calculated in hmlrad.f90, 
    !        hmlaa/hmlalo/hmllolo.f90 use haa/hloa/hlolo to calculate Hamiltonian matrix,
    !        seceqnfv.f90 calls hmlaa/hmlalo/hmllolo to setup the total Hamiltonian matrix, then diagonalize.
    !        seceqnfv.f90 calls "sethml" and "setovl" do the same setup when it's not packed format, then diagonalize. 
    !        "sethml" and "setovl" are defined in src/addons/mod_seceqn.f90
    !        (Exciting does not have "sethml" and "setovl")
    !
    ! 8. In EP, an IMPORTANT step in sethml is to multiply apwalm with haa/hloa/hlolo, 
    !    which is done in src/addons/mod_seceqn.f90
    !    In Exciting, the similar job is done in src/src_eigensystem/hamiltonandoverlapsetup.f90, 
    !    using "apwi" and haaij/haloij/hloloij
    !
    !    addons/mod_seceqn.f90 needs to be studied !!!
    !
    ! 9. In Exciting, olprad.f90: 
    !    h1aa    is for valence-relativity-iora;  
    !    h1loa   is for valence-relativity-iora,  oalo is the normal one;
    !    h1lolo  is for valence-relativity-iora, ololo is the normal one; 
    ! 
    ! 10. EP::src/addons/mod_seceqn.f90    vs.    EX::src/src_eigensystem/hamint_simplified.f90
    !                            gntyry    vs.    listgnt,indgnt 
    !                                             (but only used for angular integrals, NOT the radial integrals haa/halo/hlolointegrals)
    !     
    !     It seems mod_seceqn.f90 puts angular AND radial integrals together and calculates three parts aa/alo/lolo, 
    !     while hamint_simplified.f90 separately calculates aa/alo/lolo for radial (haa/halo/hlolointegrals), 
    !     and then calculates aa/alo/lolo for angular (haa/halo/hloloij = listgnt * haa/halo/hlolointegrals). 
    !     haa/halo/hloloij (plus oalo/ololo/h1aa/h1lo/h1lolo) are then used in hamiltonandoverlapsetup.f90. 
    !     (remember hamiltonandoverlapsetup is called in seceqnfv.f90 for non-Davidson diagonalisation)
    


    
    ! ------------------------------------------------------------------------------------------------------- 
    ! 
    ! in Exciting, in scf_cycle.f90, it calls seceqn(ik, evalfv, evecfv, evecsv).
    ! 
    ! EX::seceqn.f90 calls seceqnfv, then calls seceqnss/seceqnsv/secqnsv2. 
    ! 
    ! EX::seceqnfv.f90 1. sets up system%hamilton and system%overlap, for Davidson or not
    !
    !                     condition ".ne.Davidson" or "constructHS" ---> hamiltonandoverlapsetup 
    !                     (i.e. hamiltonandoverlapsetup explicitly sets H and O, it's not for Davidson)
    !                     
    !                  2. diagonalise with options Lapack/Davidson/Arpack, 
    !
    !                     Lapack   ---> solvewithlapack(system,nstfv,evecfv,evalfv), defined in 
    !                                   src/src_eigensystem/modfvsystem.f90 and it's LAPACK 3.0 call to zhpgvx.
    !
    !                     Davidson ---> singularcomponents(system,nstfv,evecfv,evalfv,ik)
    !                                   davidson(system,nstfv,evecfv,evalfv,ik)
    ! 
    !                     Arpack   ---> iterativearpacksecequn(system,nstfv,evecfv,evalfv,InvertMethod,sigma,arpackseed(:,ik))
    !                     
    !                  3. apply valence relativity correction to evecfv, 
    ! 
    ! 
    ! EP::seceqn.f90 has very similar structure except 
    !                 (1) it has wann stuff at the end
    !                 (2) it has seceqnfv/seceqnit instead of seceqnfv
    ! 
    ! EP::seceqnfv.f90 1. sets up Hamiltonian and Overlap, for packed or not-packed
    !                  2. diagonalise with LAPACK 3.0 call to zhpgvx. 
    ! 
    !                  "packed" is set to .false. in seceqnfv.f90, never changed anywhere else !!
    !                  
    !                  "packed" will use hmlaa/hmlalo/hmllolo/hmlistl
    !                                    olpaa/olpalo/olplolo/olpistl                  
    !
    !                  "not-packed" will call sethml and setovl, which are defined in 
    !                  src/addons/mod_seceqn.f90.
    !
    !
    !









