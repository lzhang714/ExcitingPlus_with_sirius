! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine init0

  use modmain
  use modxcifc
  use modldapu
#ifdef _SIRIUS_
  use mod_sirius
#endif

  !DESCRIPTION:
  !  Performs basic consistency checks as well as allocating and initialising
  !  global variables not dependent on the $k$-point set. 
  !
  !REVISION HISTORY:
  !  Created January 2004 (JKD)
  !  Modified March 2019, to interface to SIRIUS library

  implicit none

  ! local variables
  integer is,ia,ia1,ia2,ias,ncls,io,jo,ilo,jlo
  integer ist,l,m,lm,iv(3)
  !real(8) sum
  real(8) ts0,ts1
  logical autoenu
  integer, allocatable :: icls(:)
  
  ! zero timing variables
  timeinit=0.d0
  timemat=0.d0
  timefv=0.d0
  timesv=0.d0
  timerho=0.d0
  timepot=0.d0
  timefor=0.d0
  call timesec(ts0)

!------------------------------------!
!     angular momentum variables     !
!------------------------------------!

  lmmaxvr=(lmaxvr+1)**2
  lmmaxapw=(lmaxapw+1)**2
  if (lmaxvr.gt.lmaxapw) then
    write(*,*)
    write(*,'("Error(init0): lmaxvr > lmaxapw : ",2I8)') lmaxvr,lmaxapw
    write(*,*)
    stop
  end if
  ! index to (l,m) pairs
  if (allocated(idxlm)) deallocate(idxlm)
  allocate(idxlm(0:50,-50:50))
  lm=0
  do l=0,50
    do m=-l,l
      lm=lm+1
      idxlm(l,m)=lm
    end do
  end do
  ! array of i**l values
  if (allocated(zil)) deallocate(zil)
  allocate(zil(0:50))
  do l=0,50
    zil(l)=zi**l
  end do

!------------------------------------!
!     index to atoms and species     !
!------------------------------------!

  ! check if the system is an isolated molecule
  if (molecule) then
    primcell=.false.
    tshift=.false.
  end if
  ! find primitive cell if required
  if (primcell) call findprim
  natmmax=0
  ias=0
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=ias+1
      idxas(ia,is)=ias
    end do
    ! maximum number of atoms over all species
    natmmax=max(natmmax,natoms(is))
  end do
  ! total number of atoms
  natmtot=ias

!------------------------!
!     spin variables     !
!------------------------!

  if (spinsprl) then
    spinpol=.true.
    spinorb=.false.
    select case(task)
    case(51,52,53,61,62,63,120,121)
      write(*,*)
      write(*,'("Error(init0): spin-spirals do not work with task ",I4)') task
      write(*,*)
      stop
    end select
    if (xctype(1).lt.0) then
      write(*,*)
      write(*,'("Error(init0): spin-spirals do not work with the OEP method")')
      write(*,*)
      stop
    end if
  end if
  ! spin-orbit coupling or fixed spin moment implies spin-polarised calculation  
  if ((spinorb).or.(fixspin.ne.0).or.(spinsprl)) spinpol=.true.
  ! number of spinor components and maximum allowed occupancy  
  if (spinpol) then
    nspinor=2
    occmax=1.d0
  else
    nspinor=1
    occmax=2.d0
  end if
  ! number of spin-dependent first-variational functions per state  
  if (spinsprl) then
    nspnfv=2
  else
    nspnfv=1
  end if
  ! spin-polarised calculations require second-variational eigenvectors  
  if (spinpol) tevecsv=.true.
  ! Hartree-Fock/RDMFT requires second-variational eigenvectors
  if ((task.eq.5).or.(task.eq.6).or.(task.eq.300)) tevecsv=.true.
  ! get exchange-correlation functional data
  call getxcdata(xctype,xcdescr,xcspin,xcgrad)
  if ((spinpol).and.(xcspin.eq.0)) then
    write(*,*)
    write(*,'("Error(init0): requested spin-polarised run with spin-unpolarised exchange-correlation functional")')
    write(*,*)
    stop
  end if
  ! set the magnetic fields to the initial values
  bfieldc(:)=bfieldc0(:)
  bfcmt(:,:,:)=bfcmt0(:,:,:)
  ! check for collinearity in the z-direction and set the dimension of the magnetisation and exchange-correlation vector fields
  if (spinpol) then
    ndmag=1
    if ((abs(bfieldc(1)).gt.epslat).or.(abs(bfieldc(2)).gt.epslat)) ndmag=3
    do is=1,nspecies
      do ia=1,natoms(is)
        if ((abs(bfcmt(1,ia,is)).gt.epslat).or.(abs(bfcmt(2,ia,is)).gt.epslat)) ndmag=3
      end do
    end do
    ! source-free fields and spin-spirals are non-collinear in general
    if ((nosource).or.(spinsprl)) ndmag=3
    ! spin-orbit coupling is non-collinear in general
    if (spinorb) ndmag=3
  else
    ndmag=0
  end if
  ! spin-polarised cores
  if (.not.spinpol) spincore=.false.
  ! set fixed spin moment effective field to zero
  bfsmc(:)=0.d0
  ! set muffin-tin FSM fields to zero
  bfsmcmt(:,:,:)=0.d0

  !----------------------------------!
  !     crystal structure set up     !
  !----------------------------------!

  ! generate the reciprocal lattice vectors and unit cell volume
  call reciplat
  ! compute the inverse of the lattice vector matrix
  call r3minv(avec,ainv)
  ! compute the inverse of the reciprocal vector matrix
  call r3minv(bvec,binv)
  ! Cartesian coordinates of the spin-spiral vector
  call r3mv(bvec,vqlss,vqcss)
  do is=1,nspecies
    do ia=1,natoms(is)
      ! map atomic lattice coordinates to [0,1) if not in molecule mode
      if (.not.molecule) call r3frac(epslat,atposl(:,ia,is),iv)
      ! determine atomic Cartesian coordinates
      call r3mv(avec,atposl(:,ia,is),atposc(:,ia,is))
    end do
  end do
  ! automatically determine the muffin-tin radii if required
  if (autormt) call autoradmt
  ! check for overlapping muffin-tins
  call checkmt

  !-------------------------------!  
  !     vector fields E and A     !  ! EXCITING does not have this part
  !-------------------------------!

  efieldpol=.false.
  if ((abs(efieldc(1)).gt.epslat).or.(abs(efieldc(2)).gt.epslat).or.(abs(efieldc(3)).gt.epslat)) then
    efieldpol=.true.
    tshift=.false.
    ! electric field vector in lattice coordinates
    call r3mv(ainv,efieldc,efieldl)
  end if
  afieldpol=.false.
  if ((abs(afieldc(1)).gt.epslat).or.(abs(afieldc(2)).gt.epslat).or.(abs(afieldc(3)).gt.epslat)) then
    afieldpol=.true.
    ! vector potential added in second-variational step
    tevecsv=.true.
  end if

  !---------------------------------!
  !     crystal symmetry set up     !  ! EXCITING has this part together with crystal structure set up
  !---------------------------------!

  ! find Bravais lattice symmetries
  call findsymlat
  ! use only the identity if required
  if (nosym) nsymlat=1
  ! find the crystal symmetries and shift atomic positions if required
  call findsymcrys
  ! find the site symmetries
  call findsymsite
  ! check if fixed spin moments are invariant under the symmetry group
  call checkfsm

  !-----------------------!
  !     radial meshes     !  
  !-----------------------!

  nrmtmax=1
  nrcmtmax=1
  do is=1,nspecies
    ! make the muffin-tin mesh commensurate with lradstp
    nrmt(is)=nrmt(is)-mod(nrmt(is)-1,lradstp)
    nrmtmax=max(nrmtmax,nrmt(is))
    ! number of coarse radial mesh points
    nrcmt(is)=(nrmt(is)-1)/lradstp+1
    nrcmtmax=max(nrcmtmax,nrcmt(is))
  end do
  ! set up atomic and muffin-tin radial meshes
  call genrmesh
  
  !--------------------------------------!
  !     charges and number of states     !  ! SIRIUS starts to enter here
  !--------------------------------------!

  chgzn=0.d0
  chgcr=0.d0
  chgval=0.d0
  spnstmax=0
  do is=1,nspecies
    !nuclear charge
    chgzn=chgzn+spzn(is)*dble(natoms(is))
    !find the maximum number of atomic states
    spnstmax=max(spnstmax,spnst(is))
    !compute the electronic charge for each species, as well as the total core and valence charge
    spze(is)=0.d0
    do ist=1,spnst(is)
      spze(is)=spze(is)+spocc(ist,is)
      if (spcore(ist,is)) then
        chgcr=chgcr+dble(natoms(is))*spocc(ist,is)
      else
        chgval=chgval+dble(natoms(is))*spocc(ist,is)
      end if
    end do
  end do
  ! add excess charge
  chgval=chgval+chgexs
  ! total charge
  chgtot=chgcr+chgval
  if (chgtot.lt.1.d-8) then
    write(*,*)
    write(*,'("Error(init0): zero total charge")')
    write(*,*)
    stop
  end if
  ! effective Wigner radius
  rwigner=(3.d0/(fourpi*(chgtot/omega)))**(1.d0/3.d0)
  
  
  
  ! ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !
  !                                                     initialize SIRIUS                                                         !
  ! ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !                   
  
  ! Global on/off switch that makes ALL other switches effective.
  ! TODO: make the value assignment done in elk.in.
  use_sirius_library=.true.
  
  ! other control switches declared as parameter constant in modsirius.f90, copied here for memo. 
  !
  !use_sirius_vha=.true.
  !use_sirius_vxc=.true.
  !use_sirius_eigen_states=.true. 
  !use_sirius_density=.true.  
  !use_sirius_cfun=.true. 
  !use_sirius_gvec=.true.
  !
  !use_sirius_radial_solver=.false.   ! have to set FALSE b/c "sirius_radial_solver" seems gone from the API
  !use_sirius_apwfr=.false.           ! set FALSE b/c we want to pass apwfr to SIRIUS 
  !use_sirius_lofr=.false.            ! set FALSE b/c we want to pass lofr to SIRIUS 
  !use_sirius_olprad=.false.          ! set FALSE b/c we want to pass Olp radial integrals to SIRIUS
  !use_sirius_hmlrad=.false.          ! set FALSE b/c we want to pass Ham radial integrals to SIRIUS
  !
  !use_sirius_autoenu=.false.         ! let sirius automatically determine linearization energy, for APW and LO
  !use_sirius_rhoinit=.false.         ! false at the moment b/c it needs to be used with sirius species files
  
  ! ==================================================== 
  if (use_sirius_library) then
#ifdef _SIRIUS_
        ! init sirius (set to .false. b/c we assume "MPI_Init" is already done, otherwise it should be set to .true.)
        call sirius_initialize(logical(.false.,kind=c_bool))
        ! set global MPI communicator for sirius
        sctx = sirius_create_context(MPI_COMM_WORLD)
        ! set json string to specify calculation type fp-lapw or pp-pw
        call sirius_import_parameters(sctx, string('{"parameters" : {"electronic_structure_method" : "full_potential_lapwlo"}}'))
        ! set important parameters
        call sirius_set_parameters(sctx,&
                                  &lmax_apw=lmaxapw,&
                                  &lmax_rho=lmaxvr,&
                                  &lmax_pot=lmaxvr,&
                                  &pw_cutoff=gmaxvr,&
                                  &num_mag_dims=ndmag,&
                                  &auto_rmt=0,&
                                  &core_rel=string('none'),&
                                  &valence_rel=string('none')    )
        ! set lattice vectors
        call sirius_set_lattice_vectors( sctx, avec(1,1), avec(1,2), avec(1,3) )

        ! ----------------------------------------------------------------------- loop species
        ! In main, readinput (which called readspecies) has been called BEFORE 
        ! executing gndstate (which has init0/init1 at the beginning).
        ! So, at this point, all species-related variables: spfname,spsymb,spname,
        ! spzn,spmass,sprmin,rmt,sprmax,nrmt,spnst,spn,spl,spk,spocc,spcore  
        ! (check species file) have been declared/allocated/assigned values. 
        
        do is = 1, nspecies
          
          ! set the label of the atom type. 
          ! NOTE: species file name, e.g. "Ni.in", is spfname(is) in EP.
          !       species chemical names, e.g. "nickel" or "Ni", are spname(is) or spsymb(is) in EP.
          !       Because spsymb(is) has been used as an (optional)input argument in 
          !       "sirius_add_atom_type" below, here I set "label" to be spfname(is). 
          label = adjustl(spfname(is))
          
          ! add atom type. 
          ! NOTE: the optional argument "fname" is not specified, b/c it is for passing the .json 
          ! species file generated by sirius, which is currently experimental. 
          call sirius_add_atom_type(sctx, string(trim(label)), zn=nint(-spzn(is)), symbol=string(trim(spsymb(is))), mass=spmass(is))
          
          ! From Exciting comments: radial grid can be set from the EP side or generated by SIRIUS, 
          !                         in the second case it must be brought back to the Fortran side and all 
          !                         the relevant arrays must be reallocated. 
          ! So, at the moment, set it in EP side. 
          call sirius_set_atom_type_radial_grid(sctx, string(trim(label)), nrmt(is), spr(1, is))
          
          ! set elecronic configuration for each atom specie,
          ! (spnst - number of states, spn/spl/spk - quantum numbers)
          do ist = 1, spnst(is)          
            call sirius_set_atom_type_configuration( sctx, &
                 &string(trim(label)), spn(ist,is), spl(ist,is), spk(ist,is), spocc(ist,is), logical(spcore(ist,is),kind=c_bool) )
          enddo
          
          ! set sirius apw descriptor.
          !
          ! 1. EXCITING sets "-1" to the 3rd argument, meaning the library will figure out principle quantum number. 
          !    Here we use apwpqn(l,is), where l=0,lmaxapw. This variable is declared and read in in readspecies.f90. 
          !    EXCITING does not have this quantity.
          !
          ! 2. memo: 
          !    apwve is .true. if the linearisation energies are allowed to vary
          !    apwe0 is value of initial linearisation energy
          !    apwe  is value of linearisation energy
          !    apwdm is order of the energy derivative
          !    (these are all read in from species file)
          !
          ! 3. we can set linearisation energy (enu) to be apwe0 AND set autoenu to be true, 
          !    OR
          !    set linearisation energy (enu) to be some assigned value AND autoenu to be true or false.
          
          do l = 0, lmaxapw
            do io = 1, apword(l, is) 
              autoenu = .false.
              if (use_sirius_autoenu.and.apwve(io,l,is)) autoenu = .true.
              call sirius_add_atom_type_aw_descriptor(sctx, string(trim(label)), apwpqn(l,is), l, &
                                                      &apwe0(io, l, is), apwdm(io, l, is), autoenu)
            enddo
          enddo

          ! set sirius lo descriptor: 
          !
          ! 1. n = -1 means that the library will try to figure out principal quantum number itself.
          !    we can also use: lopqn(ilo,is), which is also species file. 
          !    Exciting has the similar variable: lorbpqn(io, ilo, is), in mod_APW_LO.F90, but I am not sure where its value comes?
          !
          ! 2. lorbve, lorbe0, lorbdm, autoenu are all same as related apw quantities (note there is no "lorbe").
          !    lorbl is lo angular momentum number.
          
          do ilo = 1, nlorb(is)
            do io = 1, lorbord(ilo, is)
              autoenu = .false.
              if (use_sirius_autoenu.and.lorbve(io, ilo, is)) autoenu = .true.
              call sirius_add_atom_type_lo_descriptor(sctx, string(trim(label)), ilo, lopqn(ilo,is), lorbl(ilo, is)&
                                                      &lorbe0(io, ilo, is), lorbdm(io, ilo, is), autoenu)
            enddo
          enddo

        enddo 
        ! ----------------------------------------------------------------------- loop species


        ! add atoms to the unit cell, 
        ! (atposl = atom position in lattice coord; bfcmt0 = starting magnetization;)
        do is = 1, nspecies
          label = adjustl(spfname(is))
          do ia = 1, natoms(is)
            call sirius_add_atom( sctx, string(trim(label)), atposl(1,ia,is), bfcmt0(1,ia,is) )
          enddo
        enddo
        
        ! set equivalent atoms
        ncls = 0
        allocate(icls(natmtot))
        icls = 0
        do is = 1, nspecies
          do ia1 = 1 ,natoms(is)
            if (icls(idxas(ia1,is)).eq.0) then
              ncls = ncls + 1
              do ia2 = 1, natoms(is) 
                if (eqatoms(ia1, ia2, is)) then
                  icls(idxas(ia2, is)) = ncls
                endif
              enddo
            endif
          enddo
        enddo
        call sirius_set_equivalent_atoms(sctx, icls(1))
        deallocate(icls)
        
        ! LDA+U is not activated for the moment.
        ! if (ldapu.ne.0) call sirius_set_uj_correction(1)
        
        ! add the XC functionals to sirius.
        ! I assume the 2nd argument in the call is LibXC names of functionals, from https://tddft.org/programs/libxc/functionals/
        call sirius_add_xc_functional(sctx, string('GGA_X_PBE'))
        call sirius_add_xc_functional(sctx, string('GGA_C_PBE'))       
        
        
        
        ! use square distribution over k-points
        ! note: this part seems for band structure plot, commented out for now
        !kpt_groups=input%groundstate%kptgroups
        !if (kpt_groups.lt.1) kpt_groups=procs
        !ranks_per_kpt=procs/kpt_groups
        !if (kpt_groups*ranks_per_kpt.ne.procs) then
        !  write(*,*) 'number of MPI ranks is incompatible with the number of kpt groups'
        !  stop
        !endif
        !cols_per_kpt=sqrt(dble(ranks_per_kpt) + 1d-8)
        !do while (mod(ranks_per_kpt,cols_per_kpt).ne.0)
        !  cols_per_kpt=cols_per_kpt-1
        !enddo
        !rows_per_kpt=ranks_per_kpt/cols_per_kpt
        
        
        
        ! as initial test, use sequential diagonalisation in lapack, i.e. sirius not running parallel.
        !mpi_grid(1)=rows_per_kpt 
        !mpi_grid(2)=cols_per_kpt
        mpi_grid(1) = 1 
        mpi_grid(2) = 1                
        call sirius_set_mpi_grid_dims(sctx, 2, mpi_grid(1))


        ! set number of fv states and aw cutoff, the same is also set in init1.f90
        ! (EXCITING has: nstfv = Int (chgval/2.d0) + input%groundstate%nempty + 1. I don't understant why +1 ?)
        nstfv = int(chgval/2.d0) + nempty
        call sirius_set_parameters(sctx, num_fv_states=nstfv, aw_cutoff=rgkmax)

        ! set the solver type, the string can be "exact" or "davidson"
        ! as initial test, use direct diagonalisation i.e. "exact"
        call sirius_set_parameters(sctx, iter_solver_type=string('exact'))

        ! initialize global variables
        call sirius_initialize_context(sctx)
        
         
         
        ! ------------------ FFT MPI related variables 
        
        ! get fft mpi communicator from sirius
        call sirius_get_fft_comm(sctx, sirius_fft_comm)                       
        ! assign value to ngr_local (local # of fft grid points)
        ngr_loc = sirius_get_num_fft_grid_points(sctx)                        
        ! assign value to sirius_fft_comm_size/rank
        call mpi_comm_size(sirius_fft_comm, sirius_fft_comm_size, ierr)
        call mpi_comm_rank(sirius_fft_comm, sirius_fft_comm_rank, ierr)    
        ! assign values to the local fraction of the array ngr_loc_all. Why shift 1??   
        allocate(ngr_loc_all(0:sirius_fft_comm_size-1))                       
        ngr_loc_all = 0                                                       
        ngr_loc_all(sirius_fft_comm_rank) = ngr_loc
        ! combine local fractions to get complete ngr_loc_all
        call mpi_allreduce(MPI_IN_PLACE, ngr_loc_all, sirius_fft_comm_size,& 
                          &MPI_INT, MPI_SUM, sirius_fft_comm, ierr)

        ! ------------------ k-points MPI related variables
        
        ! memo from Sirius documents:
        ! sirius_get_kpoint_inter_comm --- Get communicator which is used to split k-points. 
        ! sirius_get_kpoint_inner_comm --- Get communicator which is used to parallise band problem. 
        
        ! get inner/inter-k mpi communicators from sirius
        call sirius_get_kpoint_inter_comm(sctx, all_kpts_comm)
        call sirius_get_kpoint_inner_comm(sctx, kpt_inner_comm)
        ! assign value to kpt_inner_comm_size/rank
        call mpi_comm_size(kpt_inner_comm, kpt_inner_comm_size, ierr)
        call mpi_comm_rank(kpt_inner_comm, kpt_inner_comm_rank, ierr) 
        ! assign value all_kpts_comm_size/rank        
        call mpi_comm_size(all_kpts_comm, all_kpts_comm_size, ierr)
        call mpi_comm_rank(all_kpts_comm, all_kpts_comm_rank, ierr)

#else
    stop sirius_error
#endif
  endif  
  ! ==================================================== 
  
  if (.not.use_sirius_library) then
    ! EXCITING sets square distribution of k-points, 
    ! at the moment I do nothing here. 
  endif
  
  ! ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !
  !                                                     initialize SIRIUS                                                         !
  ! ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !



!-------------------------!
!     G-vector arrays     !
!-------------------------!
! determine gkmax from rgkmax and the muffin-tin radius
if (nspecies.gt.0) then
  if ((isgkmax.ge.1).and.(isgkmax.le.nspecies)) then
! use user-specified muffin-tin radius
    gkmax=rgkmax/rmt(isgkmax)
  else if (isgkmax.eq.-1) then
! use average muffin-tin radius
    sum=0.d0
    do is=1,nspecies
      sum=sum+dble(natoms(is))*rmt(is)
    end do
    sum=sum/dble(natmtot)
    gkmax=rgkmax/sum
  else
! use minimum muffin-tin radius
    gkmax=rgkmax/minval(rmt(1:nspecies))
  end if
else
  gkmax=rgkmax/2.d0
end if
gmaxvr=max(gmaxvr,2.d0*gkmax+epslat)
! find the G-vector grid sizes
call gridsize
! generate the G-vectors
call gengvec
! generate the spherical harmonics of the G-vectors
call genylmg
! allocate structure factor array for G-vectors
if (allocated(sfacg)) deallocate(sfacg)
allocate(sfacg(ngvec,natmtot))
! generate structure factors for G-vectors
call gensfacgp(ngvec,vgc,ngvec,sfacg)
! generate the characteristic function
call gencfun

!-------------------------!
!     atoms and cores     !
!-------------------------!
! allocate global species charge density and potential arrays
if (allocated(sprho)) deallocate(sprho)
allocate(sprho(spnrmax,nspecies))
if (allocated(spvr)) deallocate(spvr)
allocate(spvr(spnrmax,nspecies))
! allocate core state eigenvalue array and set to default
if (allocated(evalcr)) deallocate(evalcr)
allocate(evalcr(spnstmax,natmtot))
evalcr=-1.d0
! allocate core state radial wavefunction array
if (allocated(rwfcr)) deallocate(rwfcr)
allocate(rwfcr(spnrmax,2,spnstmax,natmtot))
! number of core spin channels
if (spincore) then
  nspncr=2
else
  nspncr=1
end if
! allocate core state charge density array
if (allocated(rhocr)) deallocate(rhocr)
allocate(rhocr(spnrmax,natmtot,nspncr))
! partial core density array
if (pt_core) then
 if (allocated(ptrhocr)) deallocate(ptrhocr)
 allocate(ptrhocr(spnrmax,natmtot,nspncr))
endif

!---------------------------------------!
!     charge density and potentials     !
!---------------------------------------!

! allocate charge density arrays
if (allocated(rhomt)) deallocate(rhomt)
allocate(rhomt(lmmaxvr,nrmtmax,natmtot))
if (allocated(rhoir)) deallocate(rhoir)
allocate(rhoir(ngrtot))
if (allocated(vmad)) deallocate (vmad)  ! added for sirius interface
allocate (vmad(natmtot))                ! added for sirius interface
! allocate charge density arrays
if (rho_val.or.pt_core) then
 if (allocated(rhomt_val)) deallocate(rhomt_val)
 allocate(rhomt_val(lmmaxvr,nrmtmax,natmtot))
 if (allocated(rhoir_val)) deallocate(rhoir_val)
 allocate(rhoir_val(ngrtot))
endif
! allocate magnetisation arrays
if (allocated(magmt)) deallocate(magmt)
if (allocated(magir)) deallocate(magir)
if (spinpol) then
  allocate(magmt(lmmaxvr,nrmtmax,natmtot,ndmag))
  allocate(magir(ngrtot,ndmag))
end if
! Coulomb potential
if (allocated(vclmt)) deallocate(vclmt)
allocate(vclmt(lmmaxvr,nrmtmax,natmtot))
if (allocated(vclir)) deallocate(vclir)
allocate(vclir(ngrtot))
! exchange-correlation potential
if (allocated(vxcmt)) deallocate(vxcmt)
allocate(vxcmt(lmmaxvr,nrmtmax,natmtot))
if (allocated(vxcir)) deallocate(vxcir)
allocate(vxcir(ngrtot))
! exchange-correlation potential from valence density
if (rho_val.or.pt_core) then
 if (allocated(vxcmt_val)) deallocate(vxcmt_val)
 allocate(vxcmt_val(lmmaxvr,nrmtmax,natmtot))
 if (allocated(vxcir_val)) deallocate(vxcir_val)
 allocate(vxcir_val(ngrtot))
endif
! exchange-correlation magnetic and effective fields
if (allocated(bxcmt)) deallocate(bxcmt)
if (allocated(bxcir)) deallocate(bxcir)
if (allocated(beffmt)) deallocate(beffmt)
if (spinpol) then
  allocate(bxcmt(lmmaxvr,nrmtmax,natmtot,ndmag))
  allocate(bxcir(ngrtot,ndmag))
  allocate(beffmt(lmmaxvr,nrcmtmax,natmtot,ndmag))
end if
! spin-orbit coupling radial function
if (allocated(socfr)) deallocate(socfr)
if (spinorb) then
  allocate(socfr(nrcmtmax,natmtot))
end if
! exchange energy density
if (allocated(exmt)) deallocate(exmt)
allocate(exmt(lmmaxvr,nrmtmax,natmtot))
if (allocated(exir)) deallocate(exir)
allocate(exir(ngrtot))
! correlation energy density
if (allocated(ecmt)) deallocate(ecmt)
allocate(ecmt(lmmaxvr,nrmtmax,natmtot))
if (allocated(ecir)) deallocate(ecir)
allocate(ecir(ngrtot))
! effective potential
if (allocated(veffmt)) deallocate(veffmt)
allocate(veffmt(lmmaxvr,nrmtmax,natmtot))
if (allocated(veffir)) deallocate(veffir)
allocate(veffir(ngrtot))
if (allocated(veffig)) deallocate(veffig)
allocate(veffig(ngvec))
! allocate muffin-tin charge and moment arrays
if (allocated(chgmt)) deallocate(chgmt)
allocate(chgmt(natmtot))
if (allocated(mommt)) deallocate(mommt)
allocate(mommt(3,natmtot))

!--------------------------------------------!
!     forces and structural optimisation     !
!--------------------------------------------!
if (allocated(forcehf)) deallocate(forcehf)
allocate(forcehf(3,natmtot))
if (allocated(forcecr)) deallocate(forcecr)
allocate(forcecr(3,natmtot))
if (allocated(forceibs)) deallocate(forceibs)
allocate(forceibs(3,natmtot))
if (allocated(forcetot)) deallocate(forcetot)
allocate(forcetot(3,natmtot))
if (allocated(forcetp)) deallocate(forcetp)
allocate(forcetp(3,natmtot))
if (allocated(tauatm)) deallocate(tauatm)
allocate(tauatm(natmtot))
! initialise the previous force
forcetp(:,:)=0.d0
! initial step sizes
tauatm(:)=tau0atm

!-------------------------!
!     LDA+U variables     !
!-------------------------!
if ((ldapu.ne.0).or.(task.eq.17)) then
! LDA+U requires second-variational eigenvectors
  tevecsv=.true.
! density matrices
  if (allocated(dmatlu)) deallocate(dmatlu)
  allocate(dmatlu(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
! potential matrix elements
  if (allocated(vmatlu)) deallocate(vmatlu)
  allocate(vmatlu(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
! zero the potential
  vmatlu(:,:,:,:,:)=0.d0
! energy for each atom
  if (allocated(engyalu)) deallocate(engyalu)
  allocate(engyalu(natmtot))
! interpolation constants (alpha)
  if (allocated(alphalu)) deallocate(alphalu)
  allocate(alphalu(natmtot))
end if

!-----------------------!
!     miscellaneous     !
!-----------------------!
! determine the nuclear-nuclear energy
call energynn
! get smearing function description
call getsdata(stype,sdescr)
! get mixing type description
call getmixdata(mixtype,mixdescr)
! generate the spherical harmonic transform (SHT) matrices
call genshtmat
! allocate 1D plotting arrays
if (allocated(dvp1d)) deallocate(dvp1d)
allocate(dvp1d(nvp1d))
if (allocated(vplp1d)) deallocate(vplp1d)
allocate(vplp1d(3,npp1d))
if (allocated(dpp1d)) deallocate(dpp1d)
allocate(dpp1d(npp1d))
! zero self-consistent loop number
iscl=0
tlast=.false.
! if reducebf < 1 then reduce the external magnetic fields immediately for
! non-self-consistent calculations
if ((reducebf.lt.1.d0).and.(task.gt.3).and.(task.ne.700)) then
  bfieldc(:)=0.d0
  bfcmt(:,:,:)=0.d0
end if
! set the Fermi energy to zero
efermi=0.d0

call timesec(ts1)
timeinit=timeinit+ts1-ts0
return
end subroutine
!EOC

