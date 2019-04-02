
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine elnes
use modmain
use modtest
implicit none
! local variables
integer ik,jk,ikq,ist,jst,ikloc
integer isym,n,nsk(3),iw
real(8) vecqc(3),qc,vkql(3)
real(8) v(3),wd,dw,w,t1
! allocatable arrays
real(8), allocatable :: e(:,:,:)
real(8), allocatable :: f(:,:,:)
real(8), allocatable :: ddcs(:)
complex(8), allocatable :: zqrmt(:,:,:)
complex(8), allocatable :: emat(:,:)
! initialise universal variables
call init0
call init1
! check q-vector is commensurate with k-point grid
v(:)=dble(ngridk(:))*vecql(:)
v(:)=abs(v(:)-nint(v(:)))
if ((v(1).gt.epslat).or.(v(2).gt.epslat).or.(v(3).gt.epslat)) then
  write(*,*)
  write(*,'("Error(elnes): q-vector incommensurate with k-point grid")')
  write(*,'(" ngridk : ",3I6)') ngridk
  write(*,'(" vecql : ",3G18.10)') vecql
  write(*,*)
  stop
end if
! allocate local arrays
allocate(e(nstsv,nstsv,nkptnr))
allocate(f(nstsv,nstsv,nkptnr))
allocate(ddcs(nwdos))
allocate(zqrmt(lmmaxvr,nrcmtmax,natmtot))
! read in the density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! get the second-variational eigenvalues and occupancies from file
do ik=1,nkpt
  call getevalsv(vkl(:,ik),evalsv(:,ik))
  call getoccsv(vkl(:,ik),occsv(:,ik))
end do
! generate the phase factor function exp(iq.r) in the muffin-tins
call genzqrmt(zqrmt)
e(:,:,:)=0.d0
f(:,:,:)=0.d0
allocate(emat(nstsv,nstsv))
! begin parallel loop over non-reduced k-points
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  write(*,'("Info(elnes): ",I6," of ",I6," k-points")') ik,nkptnr
! equivalent reduced k-point
  jk=ikmap(ivknr(1,ik),ivknr(2,ik),ivknr(3,ik))
! k+q-vector in lattice coordinates
  vkql(:)=vklnr(:,ik)+vecql(:)
! index to k+q-vector
  call findkpt(vkql,isym,ikq)
! compute < i,k+G+q | exp(iq.r) | j,k > matrix elements
  call genexpiqr(vklnr(:,ik),zqrmt,emat)
! add to the double differential scattering cross-section
  do jst=1,nstsv
    if (evalsv(jst,jk).lt.emaxelnes) then
      do ist=1,nstsv
        e(ist,jst,ik)=evalsv(ist,ikq)-evalsv(jst,jk)
        t1=dble(emat(ist,jst))**2+aimag(emat(ist,jst))**2
        f(ist,jst,ik)=t1*occsv(jst,jk)*(occmax-occsv(ist,ikq))
      end do
    end if
  end do
end do
deallocate(emat)
call mpi_grid_reduce(e(1,1,1),nstsv*nstsv*nkptnr,dims=(/dim_k/))
call mpi_grid_reduce(f(1,1,1),nstsv*nstsv*nkptnr,dims=(/dim_k/))
! number of subdivisions used for interpolation
nsk(:)=max(ngrdos/ngridk(:),1)
n=nstsv*nstsv
if (mpi_grid_root()) then
! integrate over the Brillouin zone
  call brzint(nsmdos,ngridk,nsk,ikmapnr,nwdos,wdos,n,n,e,f,ddcs)
! q-vector in Cartesian coordinates
  call r3mv(bvec,vecql,vecqc)
  qc=sqrt(vecqc(1)**2+vecqc(2)**2+vecqc(3)**2)
  t1=2.d0/(omega*occmax)
  if (qc.gt.epslat) t1=t1/qc**4
  ddcs(:)=t1*ddcs(:)
  open(50,file='ELNES.OUT',action='WRITE',form='FORMATTED')
  wd=wdos(2)-wdos(1)
  dw=wd/dble(nwdos)
  do iw=1,nwdos
    w=dw*dble(iw-1)+wdos(1)
    write(50,'(2G18.10)') w,ddcs(iw)
  end do
  close(50)
  write(*,*)
  write(*,'("Info(elnes):")')
  write(*,'(" ELNES double differential cross-section written to ELNES.OUT")')
  write(*,*)
! write ELNES distribution to test file
  call writetest(140,'ELNES cross-section',nv=nwdos,tol=1.d-2,rva=ddcs)
end if
deallocate(e,f,ddcs,zqrmt)
return
end subroutine

