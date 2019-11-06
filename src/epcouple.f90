
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine epcouple
use modmain
use modqpt
implicit none
! local variables
integer is,ia,ias,js,ja,jas
integer i,j,n,isym,ip,iv(3)
integer iq,ik,jk,ikq
integer ist,jst,ir,irc
real(8) vkql(3),x
real(8) t1,t2,t3,t4,t5
complex(8) zt1
! allocatable arrays
real(8), allocatable :: wq(:,:)
real(8), allocatable :: gq(:,:)
real(8), allocatable :: evalfv(:,:)
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: dynq(:,:,:)
complex(8), allocatable :: ev(:,:)
complex(8), allocatable :: dveffmt(:,:,:,:)
complex(8), allocatable :: dveffir(:,:)
complex(8), allocatable :: zfmt(:,:,:)
complex(8), allocatable :: gzfmt(:,:,:,:)
complex(8), allocatable :: zfir(:)
complex(8), allocatable :: epmat(:,:,:)
! external functions
real(8) sdelta,stheta
external sdelta,stheta
! initialise universal variables
call init0
call init1
call init2
! check k-point grid is commensurate with q-point grid
iv(:)=mod(ngridk(:),ngridq(:))
if ((iv(1).ne.0).or.(iv(2).ne.0).or.(iv(3).ne.0)) then
  write(*,*)
  write(*,'("Error(epcouple): k-point grid incommensurate with q-point grid")')
  write(*,'(" ngridk : ",3I6)') ngridk
  write(*,'(" ngridq : ",3I6)') ngridq
  write(*,*)
  stop
end if
n=3*natmtot
! allocate local arrays
allocate(wq(n,nqpt))
allocate(gq(n,nqpt))
allocate(dynq(n,n,nqpt))
allocate(ev(n,n))
allocate(dveffmt(lmmaxvr,nrcmtmax,natmtot,n))
allocate(dveffir(ngrtot,n))
allocate(zfmt(lmmaxvr,nrcmtmax,natmtot))
allocate(gzfmt(lmmaxvr,nrcmtmax,3,natmtot))
allocate(zfir(ngrtot))
! read in the density and potentials from file
call readstate
! read in the Fermi energy
call readfermi
! find the linearisation energies
call linengy
! set the speed of light >> 1
solsc=sol*20.d0
! new file extension for eigenvector files with c >> 1
filext='_EP.OUT'
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! compute the overlap radial integrals
call olprad
! compute the Hamiltonian radial integrals
call hmlrad
! generate muffin-tin effective magnetic fields and s.o. coupling functions
call genbeffmt
! begin parallel loop over k-points
  do ik=1,nkpt
! every thread should allocate its own arrays
    allocate(evalfv(nstfv,nspnfv))
    allocate(evecfv(nmatmax,nstfv,nspnfv))
    allocate(evecsv(nstsv,nstsv))
! solve the first- and second-variational secular equations
    call seceqn(ik,evalfv,evecfv,evecsv)
! write the eigenvectors to file
    call putevecfv(ik,evecfv)
    call putevecsv(ik,evecsv)
    deallocate(evalfv,evecfv,evecsv)
  end do
! reset the speed of light
solsc=sol
! compute the occupancies and density of states at the Fermi energy
call occupy
! read in the dynamical matrices
call readdyn(dynq)
! apply the acoustic sum rule
call sumrule(dynq)
! loop over species
do is=1,nspecies
! loop over atoms
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! convert potential to complex spherical harmonic expansion
    irc=0
    do ir=1,nrmt(is),lradstp
      irc=irc+1
      call rtozflm(lmaxvr,veffmt(:,ir,ias),zfmt(:,irc,ias))
    end do
! compute the gradients of the effective potential for the rigid-ion term
    call gradzfmt(lmaxvr,nrcmt(is),rcmt(:,is),lmmaxvr,nrcmtmax,zfmt(:,:,ias), &
     gzfmt(:,:,:,ias))
  end do
end do
! loop over phonon q-points
do iq=1,nqpt
  write(*,'("Info(epcouple): ",I6," of ",I6," q-points")') iq,nqpt
! diagonalise the dynamical matrix
  call dyndiag(dynq(:,:,iq),wq(:,iq),ev)
! loop over phonon branches
  do j=1,n
! zero any negative frequencies
    if (wq(j,iq).lt.0.d0) wq(j,iq)=0.d0
! find change effective potential for mode j
    dveffmt(:,:,:,j)=0.d0
    dveffir(:,j)=0.d0
    i=0
    do is=1,nspecies
! prefactor
      t1=2.d0*spmass(is)*wq(j,iq)
      if (t1.gt.1.d-8) then
        t1=1.d0/sqrt(t1)
      else
        t1=0.d0
      end if
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        do ip=1,3
          i=i+1
! read in the Cartesian change in effective potential
          call readdveff(iq,is,ia,ip,zfmt,zfir)
! add the rigid-ion term
          do irc=1,nrcmt(is)
            zfmt(:,irc,ias)=zfmt(:,irc,ias)-gzfmt(:,irc,ip,ias)
          end do
! multiply with eigenvector component and add to total
          zt1=t1*ev(i,j)
          do js=1,nspecies
            do ja=1,natoms(js)
              jas=idxas(ja,js)
              do irc=1,nrcmt(js)
                dveffmt(:,irc,jas,j)=dveffmt(:,irc,jas,j)+zt1*zfmt(:,irc,jas)
              end do
            end do
          end do
          dveffir(:,j)=dveffir(:,j)+zt1*zfir(:)
        end do
      end do
    end do
! multiply the interstitial potential with the characteristic function
    dveffir(:,j)=dveffir(:,j)*cfunir(:)
  end do
! zero the phonon linewidths array
  gq(:,iq)=0.d0
! begin parallel loop over non-reduced k-points
  do ik=1,nkptnr
    allocate(epmat(nstsv,nstsv,n))
! equivalent reduced k-point
    jk=ikmap(ivknr(1,ik),ivknr(2,ik),ivknr(3,ik))
! compute the electron-phonon coupling matrix elements
    call genepmat(iq,vklnr(:,ik),dveffmt,dveffir,epmat)
! k+q-vector in lattice coordinates
    vkql(:)=vklnr(:,ik)+vql(:,iq)
! index to k+q-vector
    call findkpt(vkql,isym,ikq)
    t1=twopi*wkptnr(ik)*(occmax/2.d0)
! loop over second-variational states
    do ist=1,nstsv
      x=(evalsv(ist,ikq)-efermi)/swidth
      t2=1.d0-stheta(stype,x)
! loop over phonon branches
      do i=1,n
        x=(evalsv(ist,ikq)-efermi+wq(i,iq))/swidth
        t3=1.d0-stheta(stype,x)
        do jst=1,nstsv
          if (wq(i,iq).gt.1.d-8) then
            t4=(t2-t3)/wq(i,iq)
          else
            t4=0.d0
          end if
          x=(evalsv(jst,jk)-evalsv(ist,ikq)-wq(i,iq))/swidth
          t4=t4*sdelta(stype,x)/swidth
          t5=dble(epmat(ist,jst,i))**2+aimag(epmat(ist,jst,i))**2
          gq(i,iq)=gq(i,iq)+wq(i,iq)*t1*t4*t5
        end do
      end do
    end do
    deallocate(epmat)
! end loop over k-points
  end do
! end loop over phonon q-points
end do
filext='.OUT'
! write the phonon linewidths to file
call writegamma(gq)
! write electron-phonon coupling constants to file
call writelambda(wq,gq)
deallocate(wq,gq,dynq,ev,dveffmt,dveffir)
deallocate(zfmt,gzfmt,zfir)
return
end subroutine

