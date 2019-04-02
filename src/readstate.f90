
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: readstate
! !INTERFACE:
subroutine readstate
! !USES:
use modmain
use modldapu
! !DESCRIPTION:
!   Reads in the charge density and other relevant variables from the file
!   {\tt STATE.OUT}. Checks for version and parameter compatibility.
!
! !REVISION HISTORY:
!   Created May 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
logical spinpol_
integer iostat
integer is,ia,ias,lmmax,lm,ir,jr
integer idm,ngm,i1,i2,i3,j1,j2,j3
integer version_(3),nspecies_,lmmaxvr_
integer natoms_,nrmt_(maxspecies),nrmtmax_
integer ngrid_(3),ngrtot_,ngvec_,ndmag_
integer nspinor_,ldapu_,lmmaxlu_
integer natmtot_,spnstmax_,spnrmax_
real(8) t1
! allocatable arrays
integer, allocatable :: mapir(:)
real(8), allocatable :: spr_(:,:)
real(8), allocatable :: rhomt_(:,:,:)
real(8), allocatable :: rhoir_(:)
real(8), allocatable :: rhomt_val_(:,:,:)
real(8), allocatable :: rhoir_val_(:)
real(8), allocatable :: vclmt_(:,:,:)
real(8), allocatable :: vclir_(:)
real(8), allocatable :: vxcmt_(:,:,:)
real(8), allocatable :: vxcir_(:)
real(8), allocatable :: vxcmt_val_(:,:,:)
real(8), allocatable :: vxcir_val_(:)
real(8), allocatable :: veffmt_(:,:,:)
real(8), allocatable :: veffir_(:)
real(8), allocatable :: magmt_(:,:,:,:)
real(8), allocatable :: magir_(:,:)
real(8), allocatable :: bxcmt_(:,:,:,:)
real(8), allocatable :: bxcir_(:,:)
complex(8), allocatable :: veffig_(:)
complex(8), allocatable :: vmatlu_(:,:,:,:,:)
open(50,file='STATE'//trim(filext),action='READ',form='UNFORMATTED', &
 status='OLD',iostat=iostat)
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(readstate): error opening ",A)') 'STATE'//trim(filext)
  write(*,*)
  stop
end if
read(50) version_
if ((version(1).ne.version_(1)).or.(version(2).ne.version_(2)) &
 .or.(version(3).ne.version_(3))) then
  write(*,*)
  write(*,'("Warning(readstate): different versions")')
  write(*,'(" current   : ",I3.3,".",I3.3,".",I3.3)') version
  write(*,'(" STATE.OUT : ",I3.3,".",I3.3,".",I3.3)') version_
end if
read(50) spinpol_
read(50) nspecies_
if (nspecies.ne.nspecies_) then
  write(*,*)
  write(*,'("Error(readstate): differing nspecies")')
  write(*,'(" current   : ",I4)') nspecies
  write(*,'(" STATE.OUT : ",I4)') nspecies_
  write(*,*)
  stop
end if
read(50) lmmaxvr_
read(50) nrmtmax_
allocate(spr_(nrmtmax_,nspecies))
do is=1,nspecies
  read(50) natoms_
  if (natoms(is).ne.natoms_) then
    write(*,*)
    write(*,'("Error(readstate): differing natoms for species ",I4)') is
    write(*,'(" current   : ",I4)') natoms(is)
    write(*,'(" STATE.OUT : ",I4)') natoms_
    write(*,*)
    stop
  end if
  read(50) nrmt_(is)
  read(50) spr_(1:nrmt_(is),is)
end do
read(50) ngrid_
read(50) ngvec_
read(50) ndmag_
if ((spinpol_).and.(ndmag_.ne.1).and.(ndmag_.ne.3)) then
  write(*,*)
  write(*,'("Error(readstate): invalid ndmag in STATE.OUT : ",I8)') ndmag_
  write(*,*)
  stop
end if
read(50) nspinor_
read(50) ldapu_
read(50) lmmaxlu_
ngrtot_=ngrid_(1)*ngrid_(2)*ngrid_(3)
allocate(mapir(ngrtot))
allocate(rhomt_(lmmaxvr_,nrmtmax_,natmtot))
allocate(rhoir_(ngrtot_))
allocate(vclmt_(lmmaxvr_,nrmtmax_,natmtot))
allocate(vclir_(ngrtot_))
allocate(vxcmt_(lmmaxvr_,nrmtmax_,natmtot))
allocate(vxcir_(ngrtot_))
allocate(veffmt_(lmmaxvr_,nrmtmax_,natmtot))
allocate(veffir_(ngrtot_))
allocate(veffig_(ngvec_))
if (rho_val.or.pt_core) then !new
 allocate(rhomt_val_(lmmaxvr_,nrmtmax_,natmtot))
 allocate(rhoir_val_(ngrtot_))
 allocate(vxcmt_val_(lmmaxvr_,nrmtmax_,natmtot))
 allocate(vxcir_val_(ngrtot_))
endif
! read muffin-tin density
read(50) rhomt_,rhoir_
! read Coulomb potential (spin independent)
read(50) vclmt_,vclir_
! read exchange-correlation potential
read(50) vxcmt_,vxcir_
! read effective potential
read(50) veffmt_,veffir_,veffig_
! read magnetisation and effective field
if (spinpol_) then
  allocate(magmt_(lmmaxvr_,nrmtmax_,natmtot,ndmag_))
  allocate(magir_(ngrtot_,ndmag_))
  allocate(bxcmt_(lmmaxvr_,nrmtmax_,natmtot,ndmag_))
  allocate(bxcir_(ngrtot_,ndmag_))
  read(50) magmt_,magir_
  read(50) bxcmt_,bxcir_
end if
! read LDA+U potential matrix elements
if ((ldapu.ne.0).and.(ldapu_.ne.0)) then
  allocate(vmatlu_(lmmaxlu_,lmmaxlu_,nspinor_,nspinor_,natmtot))
  read(50) vmatlu_
  lmmax=min(lmmaxlu,lmmaxlu_)
  vmatlu(:,:,:,:,:)=0.d0
  if (nspinor.eq.nspinor_) then
    vmatlu(1:lmmax,1:lmmax,:,:,:)=vmatlu_(1:lmmax,1:lmmax,:,:,:)
  else if ((nspinor.eq.1).and.(nspinor_.eq.2)) then
    vmatlu(1:lmmax,1:lmmax,1,1,:)=0.5d0*(vmatlu_(1:lmmax,1:lmmax,1,1,:) &
     +vmatlu_(1:lmmax,1:lmmax,2,2,:))
  else
    vmatlu(1:lmmax,1:lmmax,1,1,:)=vmatlu_(1:lmmax,1:lmmax,1,1,:)
    vmatlu(1:lmmax,1:lmmax,2,2,:)=vmatlu_(1:lmmax,1:lmmax,1,1,:)
  end if
  deallocate(vmatlu_)
end if
read(50) natmtot_
if (natmtot_.ne.natmtot) then
  write(*,*)
  write(*,'("Error(readstate): differing natmtot")')
  write(*,'(" current   : ",I4)') natmtot
  write(*,'(" STATE.OUT : ",I4)') natmtot_
  write(*,*)
  stop
end if
read(50) spnstmax_
if (spnstmax_.ne.spnstmax) then
  write(*,*)
  write(*,'("Error(readstate): differing spnstmax")')
  write(*,'(" current   : ",I4)') spnstmax
  write(*,'(" STATE.OUT : ",I4)') spnstmax_
  write(*,*)
  stop
end if
read(50) spnrmax_
if (spnrmax_.ne.spnrmax) then
  write(*,*)
  write(*,'("Error(readstate): differing spnrmax")')
  write(*,'(" current   : ",I4)') spnrmax
  write(*,'(" STATE.OUT : ",I4)') spnrmax_
  write(*,*)
  stop
end if
read(50) evalcr
read(50) spvr
read(50) bfcmt
if (rho_val.or.pt_core) then
 read(50) rhomt_val_,rhoir_val_
 read(50) vxcmt_val_,vxcir_val_
endif
close(50)
!---------------------------!
!     muffin-tin arrays     !
!---------------------------!
rhomt(:,:,:)=0.d0
vclmt(:,:,:)=0.d0
vxcmt(:,:,:)=0.d0
veffmt(:,:,:)=0.d0
if (rho_val.or.pt_core) then
 rhomt_val(:,:,:)=0.d0
 vxcmt_val(:,:,:)=0.d0
endif
if (spinpol) then
  magmt(:,:,:,:)=0.d0
  bxcmt(:,:,:,:)=0.d0
end if
lmmax=min(lmmaxvr,lmmaxvr_)
! interpolate the old arrays on the new radial mesh
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do lm=1,lmmax
      call rfinterp(nrmt_(is),spr_(:,is),lmmaxvr_,rhomt_(lm,1,ias),nrmt(is), &
       spr(:,is),lmmaxvr,rhomt(lm,1,ias))
      call rfinterp(nrmt_(is),spr_(:,is),lmmaxvr_,vclmt_(lm,1,ias),nrmt(is), &
       spr(:,is),lmmaxvr,vclmt(lm,1,ias))
      call rfinterp(nrmt_(is),spr_(:,is),lmmaxvr_,vxcmt_(lm,1,ias),nrmt(is), &
       spr(:,is),lmmaxvr,vxcmt(lm,1,ias))
      if (rho_val.or.pt_core) then !new
       call rfinterp(nrmt_(is),spr_(:,is),lmmaxvr_,rhomt_val_(lm,1,ias), &
       nrmt(is), spr(:,is),lmmaxvr,rhomt_val(lm,1,ias))
       call rfinterp(nrmt_(is),spr_(:,is),lmmaxvr_,vxcmt_val_(lm,1,ias), &
        nrmt(is),spr(:,is),lmmaxvr,vxcmt_val(lm,1,ias))
      endif
      call rfinterp(nrmt_(is),spr_(:,is),lmmaxvr_,veffmt_(lm,1,ias),nrmt(is), &
       spr(:,is),lmmaxvr,veffmt(lm,1,ias))
    end do
    if ((spinpol).and.(spinpol_)) then
      if (ndmag.eq.ndmag_) then
        do idm=1,ndmag
          do lm=1,lmmax
            call rfinterp(nrmt_(is),spr_(:,is),lmmaxvr_,magmt_(lm,1,ias,idm), &
             nrmt(is),spr(:,is),lmmaxvr,magmt(lm,1,ias,idm))
            call rfinterp(nrmt_(is),spr_(:,is),lmmaxvr_,bxcmt_(lm,1,ias,idm), &
             nrmt(is),spr(:,is),lmmaxvr,bxcmt(lm,1,ias,idm))
          end do
        end do
      else
        do lm=1,lmmax
          call rfinterp(nrmt_(is),spr_(:,is),lmmaxvr_,magmt_(lm,1,ias,ndmag_), &
           nrmt(is),spr(:,is),lmmaxvr,magmt(lm,1,ias,ndmag))
          call rfinterp(nrmt_(is),spr_(:,is),lmmaxvr_,bxcmt_(lm,1,ias,ndmag_), &
           nrmt(is),spr(:,is),lmmaxvr,bxcmt(lm,1,ias,ndmag))
        end do
      end if
    end if
  end do
end do
!-----------------------------!
!     interstitial arrays     !
!-----------------------------!
rhoir(:)=0.d0
vclir(:)=0.d0
vxcir(:)=0.d0
if (rho_val.or.pt_core) then !new
 rhoir_val(:)=0.d0
 vxcir_val(:)=0.d0
endif
veffir(:)=0.d0
veffig(:)=0.d0
if (spinpol) then
  magir(:,:)=0.d0
  bxcir(:,:)=0.d0
end if
! map from new grid to old
do i3=0,ngrid(3)-1
  t1=dble(i3*ngrid_(3))/dble(ngrid(3))
  j3=modulo(nint(t1),ngrid_(3))
  do i2=0,ngrid(2)-1
    t1=dble(i2*ngrid_(2))/dble(ngrid(2))
    j2=modulo(nint(t1),ngrid_(2))
    do i1=0,ngrid(1)-1
      t1=dble(i1*ngrid_(1))/dble(ngrid(1))
      j1=modulo(nint(t1),ngrid_(1))
      ir=i3*ngrid(2)*ngrid(1)+i2*ngrid(1)+i1+1
      jr=j3*ngrid_(2)*ngrid_(1)+j2*ngrid_(1)+j1+1
      mapir(ir)=jr
    end do
  end do
end do
do ir=1,ngrtot
  jr=mapir(ir)
  rhoir(ir)=rhoir_(jr)
  vclir(ir)=vclir_(jr)
  vxcir(ir)=vxcir_(jr)
  if (rho_val.or.pt_core) then
   rhoir_val(ir)=rhoir_val_(ir)
   vxcir_val(ir)=vxcir_val_(ir)
  endif
  veffir(ir)=veffir_(jr)
end do
ngm=min(ngvec,ngvec_)
veffig(1:ngm)=veffig_(1:ngm)
if ((spinpol).and.(spinpol_)) then
  do ir=1,ngrtot
    jr=mapir(ir)
    if (ndmag.eq.ndmag_) then
      magir(ir,:)=magir_(jr,:)
      bxcir(ir,:)=bxcir_(jr,:)
    else
      magir(ir,ndmag)=magir_(jr,ndmag_)
      bxcir(ir,ndmag)=bxcir_(jr,ndmag_)
    end if
  end do
end if
deallocate(mapir,spr_,rhomt_,rhoir_,vclmt_,vclir_)
deallocate(vxcmt_,vxcir_,veffmt_,veffir_,veffig_)
if (rho_val.or.pt_core) deallocate(vxcmt_val_,vxcir_val_,rhomt_val_,rhoir_val_)
if (spinpol_) deallocate(magmt_,magir_,bxcmt_,bxcir_)
return
end subroutine
!EOC
