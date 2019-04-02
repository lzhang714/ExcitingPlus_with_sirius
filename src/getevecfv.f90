
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: getevecfv
! !INTERFACE:
subroutine getevecfv(vpl,vgpl,evecfv)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   vpl    : p-point vector in lattice coordinates (in,real(3))
!   vgpl   : G+p-vectors in lattice coordinates (out,real(3,ngkmax))
!   evecfv : first-variational eigenvectors (out,complex(nmatmax,nstfv,nspnfv))
! !DESCRIPTION:
!   Reads in a first-variational eigenvector from file. If the input $k$-point,
!   ${\bf p}$, is not in the reduced set, then the eigenvector of the equivalent
!   point is read in and the required rotation/translation operations applied.
!
! !REVISION HISTORY:
!   Created Feburary 2007 (JKD)
!   Fixed transformation error, October 2007 (JKD, Anton Kozhevnikov)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: vpl(3)
real(8), intent(in) :: vgpl(3,ngkmax,nspnfv)
complex(8), intent(out) :: evecfv(nmatmax,nstfv,nspnfv)
! local variables
integer isym,lspl,ilspl
integer ilo,l,lm,i,j,ispn
integer ik,igp,igk,ig
integer is,ia,ja,ias,jas
integer recl,nmatmax_,nstfv_,nspnfv_
real(8) vkl_(3),v(3),t1
real(8) si(3,3)
complex(8) zt1
! allocatable arrays
complex(8), allocatable :: evecfvt(:,:)
real(8), allocatable :: vgkl_(:,:,:)
integer, allocatable :: igkig_(:,:)
allocate(vgkl_(3,ngkmax,nspnfv))
allocate(igkig_(ngkmax,nspnfv))
evecfv=zzero
! find the equivalent k-point number and crystal symmetry element
call findkpt(vpl,isym,ik)
! find the record length
inquire(iolength=recl) vkl_,nmatmax_,nstfv_,nspnfv_,evecfv,vgkl_,igkig_
open(70,file=trim(scrpath)//'EVECFV'//trim(filext),action='READ', &
 form='UNFORMATTED',access='DIRECT',recl=recl)
read(70,rec=ik) vkl_,nmatmax_,nstfv_,nspnfv_,evecfv,vgkl_,igkig_
close(70)
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getevecfv): differing vectors for k-point ",I8)') ik
  write(*,'(" current    : ",3G18.10)') vkl(:,ik)
  write(*,'(" EVECFV.OUT : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (nmatmax.ne.nmatmax_) then
  write(*,*)
  write(*,'("Error(getevecfv): differing nmatmax for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nmatmax
  write(*,'(" EVECFV.OUT : ",I8)') nmatmax_
  write(*,*)
  stop
end if
if (nstfv.ne.nstfv_) then
  write(*,*)
  write(*,'("Error(getevecfv): differing nstfv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nstfv
  write(*,'(" EVECFV.OUT : ",I8)') nstfv_
  write(*,*)
  stop
end if
if (nspnfv.ne.nspnfv_) then
  write(*,*)
  write(*,'("Error(getevecfv): differing nspnfv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nspnfv
  write(*,'(" EVECFV.OUT : ",I8)') nspnfv_
  write(*,*)
  stop
end if
! if p = k then return
t1=abs(vpl(1)-vkl(1,ik))+abs(vpl(2)-vkl(2,ik))+abs(vpl(3)-vkl(3,ik))
if (t1.lt.epslat) then
  deallocate(vgkl_,igkig_)
  return
end if
allocate(evecfvt(nmatmax,nstfv))
! index to spatial rotation in lattice point group
lspl=lsplsymc(isym)
! the inverse of the spatial symmetry rotates k into p
ilspl=isymlat(lspl)
si(:,:)=dble(symlat(:,:,ilspl))
!-----------------------------------------------!
!     translate and rotate APW coefficients     !
!-----------------------------------------------!
do ispn=1,nspnfv
  do igk=1,ngk(ispn,ik)
    ig=igkig_(igk,ispn)
    v(:)=dble(ivg(:,ig))
    t1=-twopi*dot_product(v(:),vtlsymc(:,isym))
    zt1=cmplx(cos(t1),sin(t1),8)
    evecfvt(igk,:)=zt1*evecfv(igk,:,ispn)
  end do
! inverse rotation used because transformation is passive
  do igk=1,ngk(ispn,ik)
    call r3mtv(si,vgkl_(:,igk,ispn),v)
    do igp=1,ngk(ispn,ik)
      t1=abs(v(1)-vgpl(1,igp,ispn)) &
        +abs(v(2)-vgpl(2,igp,ispn)) &
        +abs(v(3)-vgpl(3,igp,ispn))
      if (t1.lt.epslat) then
        evecfv(igp,:,ispn)=evecfvt(igk,:)
        goto 10
      end if
    end do
10 continue
  end do
end do
!---------------------------------------------------------!
!     translate and rotate local-orbital coefficients     !
!---------------------------------------------------------!
if (nlotot.gt.0) then
! rotate k-point by inverse symmetry matrix
  call r3mtv(si,vkl(:,ik),v)
! loop over the first-variational spins
  do ispn=1,nspnfv
! make a copy of the local-orbital coefficients
    do i=ngk(ispn,ik)+1,nmat(ispn,ik)
      evecfvt(i,:)=evecfv(i,:,ispn)
    end do
    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)
! equivalent atom for this symmetry
        ja=ieqatom(ia,is,isym)
        jas=idxas(ja,is)
! phase factor from translation
        t1=-twopi*dot_product(vkl(:,ik),atposl(:,ja,is))
        zt1=cmplx(cos(t1),sin(t1),8)
        t1=twopi*dot_product(v(:),atposl(:,ia,is))
        zt1=zt1*cmplx(cos(t1),sin(t1),8)
! rotate local orbitals (active transformation)
        do ilo=1,nlorb(is)
          l=lorbl(ilo,is)
          lm=idxlm(l,-l)
          i=ngk(ispn,ik)+idxlo(lm,ilo,ias)
          j=ngk(ispn,ik)+idxlo(lm,ilo,jas)
          call rotzflm(symlatc(:,:,lspl),l,l,nstfv,nmatmax,evecfvt(j,1), &
           evecfv(i,1,ispn))
          evecfv(i:i+2*l,:,ispn)=zt1*evecfv(i:i+2*l,:,ispn)
        end do
      end do
    end do
  end do
end if
deallocate(evecfvt)
deallocate(vgkl_,igkig_)
return
end subroutine
!EOC

