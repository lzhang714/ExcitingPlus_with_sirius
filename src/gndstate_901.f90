! 
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

subroutine gndstate_901
  
  use modmain
  use mod_wannier
  use mod_sic
  use mod_seceqn
  use mod_addons
#ifdef _MPI_
  use mpi
#endif
  
  implicit none

  ! ------------------------------------------------------------------------------------------------- 
  ! L.Z., 2019 Nov. 02 
  ! 1. This task will only do init0 ---> init1 ---> SIRIUS Full SCF.
  ! 2. This task does not involve parsing of eigen vectors, all related arrays are not used.  
  ! 3. "evalfv" is re-named to "evalfvloc". 
  ! ------------------------------------------------------------------------------------------------- 

  ! local variables
  logical exist
  integer ik,is,ia,ias,ist,i,j,k,ikloc,idm,inmat
  integer n,nwork
  real(8) dv,etp,de,timetot
  real(16) x1,y1,length1,x2,y2,length2
  ! allocatable arrays
  real(8), allocatable :: v(:)
  real(8), allocatable :: work(:)
  character(100) :: label
  ! 
  ! new, two temporary holders to keep the fv eigen values/vectors from SIRIUS,  
  !      they are supposed to have an additional dimension for "nspnfv", which is 
  !      missing at the moment, so this only works for spinsprl=.false. !!  
  !
  real(8),    allocatable :: evalfv_tmp(:,:)   
  complex(8), allocatable :: evecfv_tmp(:,:,:) 
  
  call timer_start(t_init,reset=.true.)  
  call init0                          
  call init1_plainmpi                    
  call timer_stop(t_init)
  
  !!! -----------------------------------------------------
  !!! set nkptloc=nkpt again (already done in init1.f90)
  nkptloc = nkpt 
  nkptnrloc = nkptnr  
  !!! if (.not.mpi_grid_in()) return
  !!! wproc = mpi_grid_root()
  if(iproc.eq.0) wproc=.true.
  !!! -----------------------------------------------------
  
  ! allocate eigen value/vector arrays
  ! the 2 tmp arrays are for all k-points
  allocate(evalfv_tmp(nstfv,nkpt))                 
  allocate(evecfv_tmp(nmatmax,nstfv,nkpt))
  ! "tsveqn" assigned .true. in mod_addons.f90, we haven't met a case of tsveqn=.false.
  if (tsveqn) then 
    allocate(evalfvloc(nstfv,nspnfv,nkptloc))           ! declared addons/mod_addons.f90
    allocate(evecfvloc(nmatmax,nstfv,nspnfv,nkptloc))   ! declared addons/mod_addons.f90
    allocate(evecsvloc(nstsv,nstsv,nkptloc))            ! declared addons/mod_addons.f90
  else
    allocate(evecfdloc(nspinor*nmatmax,nstsv,nkptloc))  ! declared addons/mod_addons.f90
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
    !! LZ took it out b/c it's taking too long time for large molecules
    !!call writenn         ! defined in addons/writenn.f90         
    !!call writegshells    ! defined in addons/writegshells.f90
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
      call allatoms                             !! solve for isolated atom density sprho
      do ias=1,natmtot
        is=ias2is(ias)
        do ist=1,spnst(is)
          evalcr(ist,ias)=speval(ist,is)
        enddo
      enddo

      if (use_sirius_library.and.sirius_run_full_scf) then
        ! set the free atomic desnity for each atom type.
        ! we can do it only here after we called allatoms() to solve the isolated atom
        ! and generate free-atom density sprho.
        do is=1,nspecies
          label = adjustl(spfname(is))
          call sirius_add_atom_type_radial_function(sctx, string(trim(label)), string("ae_rho"),sprho(1, is), spnr(is))
        enddo
        ! ------------------------ ! 
        ! SIRIUS full scf 
        write(*,*) ' SIRIUS full scf started! '
        call sirius_find_ground_state(gs_handler, potential_tol=epspot, energy_tol=epsengy, niter=maxscl,&
                                     &save_state=bool(.false.))
        write(*,*) ' SIRIUS full scf done! '
      else
        write(*,*) ' Error : task 901 must call SIRIUS to run full scf. '
      endif
      
      call poteff      ! has use_sirius_vha and use_sirius_vxc
      call timer_stop(t_rhoinit)
  
  end if
  ! ---------------------- task selection 
  
  if (wproc) call flushifc(ioinfo)
  
  tstop=.false.   
  tlast=.false.   
  etp=0.d0

  
  ! -------------------------------------------------------------------------------------- !
  ! before EP SCF loop, get fv/sv eigen values and vectors from SIRIUS and save them in: 
  ! evalfv_sirius(nstfv,nspnfv,nkpt)
  ! evalsv_sirius(nstsv,nkpt)
  ! evecfv_sirius(nmatmax,nstfv,nspnfv,nkptloc)
  ! 
  ! option 1:  master rank loop all k-points, save it only in master copy of evalfv_sirius. 
  ! option 2:  loop local fraction of k-points, index back to global "ik" and save, then  
  !            do the mpi_grid_reduce. every mpi rank will have a copy of evalfv_sirius. 
  ! ---------------------------------------------------------------------------------------!

!!  ! option 2
!!  if (mpi_grid_side(dims=(/dim_k/))) then
!!    do i=0,mpi_grid_dim_size(dim_k)-1
!!      if (mpi_grid_dim_pos(dim_k).eq.i) then
!!        do ikloc=1,nkptloc


! evalfv_allk_sirius (nstfv,nspnfv,nkpt)
! evalfv_allk        (nstfv,nspnfv,nkpt)
! evecfv_allk_sirius (nmatmax,nstfv,nspnfv,nkpt)
! evecfv_allk        (nmatmax,nstfv,nspnfv,nkpt)
! 
! evalsv_allk_sirius (nstsv,nkpt)
! evalsv             (nstsv,nkpt)
! evecsv_allk_sirius (nstsv,nstsv,nkpt)
! evecsv_allk        (nstsv,nstsv,nkpt)




!        ! get local fraction of eigen-vectors
!        do ikloc=1,nkptloc
!          ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
!          call sirius_get_fv_eigen_vectors(ks_handler, ik, evecfvloc(1, 1, 1, ikloc), nmatmax, nstfv)
!          call sirius_get_sv_eigen_vectors(ks_handler, ik, evecsvloc(1, 1, ikloc), nstsv)
!          write(*,*)'evecfv=',evecfvloc(:,:,:,ikloc)
!          write(*,*)'evecsv=',evecsvloc(:,:,ikloc)
!        enddo !ikloc
!        ! get all eigen-values and band occupancies
!        do ik = 1, nkpt
!          if (ndmag.eq.0.or.ndmag.eq.3) then
!            call sirius_get_band_energies(ks_handler, ik, 0, evalsv(1, ik))
!            call sirius_get_band_occupancies(ks_handler, ik, 0, occsv(1, ik))
!          else
!            call sirius_get_band_energies(ks_handler, ik, 0, evalsv(1, ik))
!            call sirius_get_band_energies(ks_handler, ik, 1, evalsv(nstfv+1, ik))
!            call sirius_get_band_occupancies(ks_handler, ik, 0, occsv(1, ik))
!            call sirius_get_band_occupancies(ks_handler, ik, 1, occsv(nstfv+1, ik))
!          endif
!          write(*,*)'evalsv=',evalsv(:,ik)
!          write(*,*)'occsv=',occsv(:,ik)
!        enddo




!!! 2019 Oct 28 
!!!        ! zero the tmp
!!!        evalfv_tmp(:,:)=0.d0
!!!        evecfv_tmp(:,:,:)=zzero
!!!        do ikloc=1,nkptloc
!!!        ! zero the tmp
!!!          ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
!!!          ! SIRIUS API function does not include the dimension for "nspnfv".
!!!          call sirius_get_fv_eigen_values(  ks_handler, ik, evalfv_tmp(1,ik), nstfv )                         ! eval fv
!!!          call sirius_get_fv_eigen_vectors( ks_handler, ik, evecfv_allk_sirius(1,1,1,ik), nmatmax, nstfv )    ! evec fv 
!!!          call sirius_get_sv_eigen_vectors(ks_handler, ik, evecsvloc(1, 1, ikloc), nstsv)                     ! evec sv loc
!!!          !!call sirius_get_sv_eigen_vectors( ks_handler, ik, evecsv_allk_sirius(1,1,ik), nstsv)                ! evec sv allk
!!!          ! retrieve sv eigen values and vectors from sirius, by global index ik
!!!          if (ndmag.eq.0.or.ndmag.eq.3) then
!!!            call sirius_get_band_energies   (ks_handler, ik, 0, evalsv_allk_sirius(1,       ik) )             ! eval sv
!!!            call sirius_get_band_occupancies(ks_handler, ik, 0,              occsv(1,       ik) )             ! occ
!!!          else
!!!            call sirius_get_band_energies   (ks_handler, ik, 0, evalsv_allk_sirius(1,       ik) )
!!!            call sirius_get_band_energies   (ks_handler, ik, 1, evalsv_allk_sirius(nstfv+1, ik) )
!!!            call sirius_get_band_occupancies(ks_handler, ik, 0,              occsv(1,       ik) )
!!!            call sirius_get_band_occupancies(ks_handler, ik, 1,              occsv(nstfv+1, ik) )
!!!          endif
!!!        enddo
        
        
        
!!        end do ! ikloc
!!      end if
!!      call mpi_grid_barrier(dims=(/dim_k/))
!!    end do ! i
!!  endif



! evalfv_allk_sirius (nstfv,nspnfv,nkpt)
! evalfv_allk        (nstfv,nspnfv,nkpt)
! evecfv_allk_sirius (nmatmax,nstfv,nspnfv,nkpt)
! evecfv_allk        (nmatmax,nstfv,nspnfv,nkpt)
! 
! evalsv_allk_sirius (nstsv,nkpt)
! evalsv             (nstsv,nkpt)
! evecsv_allk_sirius (nstsv,nstsv,nkpt)
! evecsv_allk        (nstsv,nstsv,nkpt)
! 
! evalfv_tmp(nstfv,nkpt)                 
! evecfv_tmp(nmatmax,nstfv,nkpt)   


!!! 2019 Oct 28 
!!!  ! combine fractions of k-points
!!!  call mpi_grid_reduce( evalfv_tmp(1,1),             nstfv*nkpt,                 dims=(/dim_k/), all=.true. )
!!!  call mpi_grid_reduce( evecfv_allk_sirius(1,1,1,1), nmatmax*nstfv*nspnfv*nkpt,  dims=(/dim_k/), all=.true. )
!!!  call mpi_grid_reduce( evalsv_allk_sirius(1,1),     nstsv*nkpt,                 dims=(/dim_k/), all=.true. )
!!!  call mpi_grid_reduce( evecsv_allk_sirius(1,1,1),   nstsv*nstsv*nkpt,           dims=(/dim_k/), all=.true. )
!!!
!!!  ! pass the tmp arrays, then free the tmp arrays, only considering non-spinsprl for now
!!!  if (.not.spinsprl) then 
!!!      do ik=1,nkpt
!!!        do i=1,nstfv
!!!          evalfv_allk_sirius(i,1,ik) = evalfv_tmp(i,ik)  
!!!        enddo
!!!      enddo
!!!  else 
!!!      write(*,*) ' EP-SIRIUS interface error, it is not yet working with spinsprl=.true. '  
!!!  endif
  deallocate(evalfv_tmp)
  deallocate(evecfv_tmp)
  
  

  ! ================================================================== ! 
  !                            EP SCF loop starts                      ! 
  ! ================================================================== ! 
  ! 
  ! EP SCF iteration is deleted.  
  ! 
  ! ================================================================== ! 
  !                            EP SCF loop ends                        ! 
  ! ================================================================== ! 





  write(*,*) ' the whole gndstate scf loop done! '

  ! ----------------------------------------------------------------------------------------------------------
  ! collect EP fv eigen values for all k-points, for comparing with evalfv_sirius. 
  ! 
  ! NOTE: evalsv for all k-points is already there, done in the k-loop inside the scf loop, it can be 
  !       directly used to compare with evalsv_allk_sirius. 
  !       
  !       This is because there is local-index-to-global-index mapping in seceqnsv.f90 (see "mpi_grid_map"
  !       in the generated mod_mpi_grid.f90), we can use "mpi_grid_reduce" (in scf loop k loop) to get 
  !       evalsv of all k-points. 
  !       
  !       The fv eigen values are always handled with local fraction of k-points, so we can't use 
  !       "mpi_grid_reduce" directly. Again we index back to global "ik" and put it in the local 
  !       copy of the arrays for all k-points, then use "mpi_grid_reduce" on evalfv. 
  ! 
  !       The same idea for fv eigen vectors, index back from "ikloc" to global "ik" (from the existing 
  !       "evecfvloc" to the array evecfv_allk. 
  ! ----------------------------------------------------------------------------------------------------------
  ! 
  !   evalfv_allk_sirius (nstfv,nspnfv,nkpt)
  !   evalfv_allk        (nstfv,nspnfv,nkpt)
  !   evecfv_allk_sirius (nmatmax,nstfv,nspnfv,nkpt)
  !   evecfv_allk        (nmatmax,nstfv,nspnfv,nkpt)
  ! 
  !   evalsv_allk_sirius (nstsv,nkpt)
  !   evalsv             (nstsv,nkpt)
  !   evecsv_allk_sirius (nstsv,nstsv,nkpt)
  !   evecsv_allk        (nstsv,nstsv,nkpt)
  ! 
  !   evalfv_tmp(nstfv,nkpt)                 
  !   evecfv_tmp(nmatmax,nstfv,nkpt)
  ! 
  !   evalfvloc(nstfv,nspnfv,nkptloc)
  !   evecfvloc(nmatmax,nstfv,nspnfv,nkptloc)
  !   evecsvloc(nstsv,nstsv,nkptloc)
  ! 
  ! ----------------------------------------------------------------------------------------------------------

!!!  evalfv_allk(:,:,:)=0.d0          ! all nkpt 
!!!  evecfv_allk(:,:,:,:)=zzero       ! all nkpt 
!!!  evecsv_allk(:,:,:)=zzero         ! all nkpt

!!!  !!  if (mpi_grid_side(dims=(/dim_k/))) then
!!!  !!    do i=0,mpi_grid_dim_size(dim_k)-1
!!!  !!      if (mpi_grid_dim_pos(dim_k).eq.i) then
!!!          do ikloc=1,nkptloc
!!!              ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
!!!              !do i=1,nstfv
!!!              !  do j=1,nspnfv
!!!              !    evalfv_allk(i,j,ik) = evalfvloc(i,j,ikloc)
!!!              !    do k=1,nmatmax
!!!              !      evecfv_allk(k,i,j,ik) = evecfvloc(k,i,j,ikloc)
!!!              !    enddo
!!!              !  enddo
!!!              !enddo            
!!!              evalfv_allk(:,:,ik)   = evalfvloc(:,:,ikloc)
!!!              evecfv_allk(:,:,:,ik) = evecfvloc(:,:,:,ikloc)
!!!              evecsv_allk(:,:,ik)   = evecsvloc(:,:,ikloc)
!!!          end do
!!!  !!      end if
!!!  !!      call mpi_grid_barrier(dims=(/dim_k/))
!!!  !!    end do
!!!  !!  endif 
!!!  
!!!    call mpi_grid_reduce( evalfv_allk(1,1,1),    nstfv*nspnfv*nkpt,         dims=(/dim_k/), all=.true. )
!!!    call mpi_grid_reduce( evecfv_allk(1,1,1,1),  nmatmax*nstfv*nspnfv*nkpt, dims=(/dim_k/), all=.true. )  
!!!    call mpi_grid_reduce( evecsv_allk(1,1,1),    nstsv*nstsv*nkpt,          dims=(/dim_k/), all=.true. )    


  
  
!!!  ! write eigen values/vectors and occupancies to file
!!!  if (mpi_grid_side(dims=(/dim_k/))) then
!!!    do i=0,mpi_grid_dim_size(dim_k)-1
!!!      if (mpi_grid_dim_pos(dim_k).eq.i) then
!!!        do ikloc=1,nkptloc
!!!          ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
!!!          if (tsveqn) then
!!!            call putevalfv(ik,evalfvloc(1,1,ikloc))      ! had renamed evalfv to evalfvloc
!!!            call putevecfv(ik,evecfvloc(1,1,1,ikloc))
!!!            call putevecsv(ik,evecsvloc(1,1,ikloc))
!!!          else
!!!            call putevecfd(ikloc,evecfdloc(1,1,ikloc))
!!!          endif
!!!          call putevalsv(ik,evalsv(1,ik))
!!!          call putoccsv(ik,occsv(1,ik))
!!!        end do
!!!      end if
!!!      call mpi_grid_barrier(dims=(/dim_k/))
!!!    end do
!!!  endif
!!!  call mpi_grid_bcast(tstop)
  
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
    deallocate(evalfvloc)      ! had renamed evalfv to evalfvloc
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


