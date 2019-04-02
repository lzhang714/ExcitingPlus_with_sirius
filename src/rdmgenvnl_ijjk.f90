
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rdmgenvnl_ijjk
! !INTERFACE:
subroutine rdmgenvnl_ijjk(ikp,vnlijji,vnlijjk)
! !USES:
use modrdm
use modmain
use modqpt
! !INPUT/OUTPUT PARAMETERS:
!   ikp     : k-point from non-reduced k-point set (in,integer)
!   vnlijji : non-local Coulomb matrix elements (out,real(nstsv,nstsv,nkpt))
!   vnlijjk : non-local Coulomb matrix elements
!             (out,complex(nstsv,nstsv,nstsv,nkpt))
! !DESCRIPTION:
!   Calculates non-local Coulomb matrix elements of the type $(i-jj-k)$ and
!   $(i-jj-i)$.
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ikp
real(8), intent(out) :: vnlijji(nstsv,nstsv,nkpt)
complex(8), intent(out) :: vnlijjk(nstsv,nstsv,nstsv,nkpt)
! local variables
integer ngknr,ik,ist1,ist2,ist3
integer lmax,ig,iq,igq0,iv(3)
real(8) cfq,v(3),t1
complex(8) zrho01,zrho02,zt1,zt2
! automatic arrays
real(8) zn(nspecies)
complex(8) sfacgq0(natmtot)
! allocatable arrays
integer, allocatable :: igkignr(:)
real(8), allocatable :: vgklnr(:,:)
real(8), allocatable :: vgkcnr(:,:)
real(8), allocatable :: gkcnr(:)
real(8), allocatable :: tpgkcnr(:,:)
real(8), allocatable :: vgqc(:,:)
real(8), allocatable :: tpgqc(:,:)
real(8), allocatable :: gqc(:)
real(8), allocatable :: jlgqr(:,:,:)
real(8), allocatable :: jlgq0r(:,:,:)
real(8), allocatable :: evalsvp(:)
real(8), allocatable :: evalsvnr(:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: sfacgknr(:,:)
complex(8), allocatable :: ylmgq(:,:)
complex(8), allocatable :: sfacgq(:,:)
complex(8), allocatable :: wfmt1(:,:,:,:,:)
complex(8), allocatable :: wfmt2(:,:,:,:,:)
complex(8), allocatable :: wfir1(:,:,:)
complex(8), allocatable :: wfir2(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:)
complex(8), allocatable :: zrhoir(:)
complex(8), allocatable :: zvclmt(:,:,:)
complex(8), allocatable :: zvclir(:)
! external functions
complex(8) zfinp
external zfinp
! allocate local arrays
allocate(igkignr(ngkmax))
allocate(vgklnr(3,ngkmax))
allocate(vgkcnr(3,ngkmax))
allocate(gkcnr(ngkmax))
allocate(tpgkcnr(2,ngkmax))
allocate(vgqc(3,ngvec))
allocate(tpgqc(2,ngvec))
allocate(gqc(ngvec))
allocate(jlgqr(0:lmaxvr+npsden+1,ngvec,nspecies))
allocate(jlgq0r(0:lmaxvr,nrcmtmax,nspecies))
allocate(evalsvp(nstsv))
allocate(evalsvnr(nstsv))
allocate(sfacgknr(ngkmax,natmtot))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv))
allocate(evecsv(nstsv,nstsv))
allocate(ylmgq(lmmaxvr,ngvec))
allocate(sfacgq(ngvec,natmtot))
allocate(wfmt1(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfmt2(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir1(ngrtot,nspinor,nstsv))
allocate(wfir2(ngrtot,nspinor,nstsv))
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot))
allocate(zrhoir(ngrtot))
allocate(zvclmt(lmmaxvr,nrcmtmax,natmtot))
allocate(zvclir(ngrtot))
! factor for long-range term
cfq=0.5d0*(omega/pi)**2
! set the point charges to zero
zn(:)=0.d0
! start loop over reduced k-point set
do ik=1,nkpt
! get the eigenvectors and values from file
  call getevalsv(vkl(:,ik),evalsvp)
  call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv(vkl(:,ik),evecsv)
! find the matching coefficients
  call match(ngk(1,ik),gkc(:,1,ik),tpgkc(:,:,1,ik),sfacgk(:,:,1,ik),apwalm)
! calculate the wavefunctions for all states for the reduced k-point
  call genwfsv(.false.,.false.,ngk(1,ik),igkig(:,1,ik),evalsvp,apwalm,evecfv, &
   evecsv,wfmt1,wfir1)
! generate G+k vectors for non-reduced k-point ikp
  call gengpvec(vklnr(:,ikp),vkcnr(:,ikp),ngknr,igkignr,vgklnr,vgkcnr,gkcnr, &
   tpgkcnr)
! get the eigenvalues/vectors from file for non-reduced k-point ikp
  call getevalsv(vklnr(:,ikp),evalsvnr)
  call getevecfv(vklnr(:,ikp),vgklnr,evecfv)
  call getevecsv(vklnr(:,ikp),evecsv)
! generate the structure factors
  call gensfacgp(ngknr,vgkcnr,ngkmax,sfacgknr)
! find the matching coefficients
  call match(ngknr,gkcnr,tpgkcnr,sfacgknr,apwalm)
! determine q-vector
  iv(:)=ivk(:,ik)-ivknr(:,ikp)
  iv(:)=modulo(iv(:),ngridk(:))
  iq=iqmap(iv(1),iv(2),iv(3))
  v(:)=vkc(:,ik)-vkcnr(:,ikp)
  do ig=1,ngvec
! determine G+q vectors
    vgqc(:,ig)=vgc(:,ig)+v(:)
! G+q-vector length and (theta, phi) coordinates
    call sphcrd(vgqc(:,ig),gqc(ig),tpgqc(:,ig))
! spherical harmonics for G+q-vectors
    call genylm(lmaxvr,tpgqc(:,ig),ylmgq(:,ig))
  end do
! structure factors for G+q
  call gensfacgp(ngvec,vgqc,ngvec,sfacgq)
! find the shortest G+q-vector
  call findigp0(ngvec,gqc,igq0)
  sfacgq0(:)=sfacgq(igq0,:)
! compute the required spherical Bessel functions
  lmax=lmaxvr+npsden+1
  call genjlgpr(lmax,gqc,jlgqr)
  call genjlgq0r(gqc(igq0),jlgq0r)
! calculate the wavefunctions for all states for passed non-reduced k-point ikp
  call genwfsv(.false.,.false.,ngknr,igkignr,evalsvnr,apwalm,evecfv,evecsv, &
   wfmt2,wfir2)
!----------------------------------------------!
!     valence-valence-valence contribution     !
!----------------------------------------------!
  do ist1=1,nstsv
    do ist2=1,nstsv
! calculate the complex overlap density
      call vnlrho(.true.,wfmt2(:,:,:,:,ist2),wfmt1(:,:,:,:,ist1), &
       wfir2(:,:,ist2),wfir1(:,:,ist1),zrhomt,zrhoir)
! compute the potential and G=0 coefficient of the density
      call zpotcoul(nrcmt,nrcmtmax,nrcmtmax,rcmt,igq0,gqc,jlgqr,ylmgq,sfacgq, &
       zn,zrhomt,zrhoir,zvclmt,zvclir,zrho02)
      zt1=zfinp(.true.,zrhomt,zvclmt,zrhoir,zvclir)
      t1=cfq*wiq2(iq)*(dble(zrho02)**2+aimag(zrho02)**2)
      vnlijjk(ist1,ist1,ist2,ik)=wkptnr(ikp)*dble(zt1)+t1
      do ist3=1,nstsv
        if (ist1.lt.ist3) then
! calculate the complex overlap density
          call vnlrho(.true.,wfmt2(:,:,:,:,ist2),wfmt1(:,:,:,:,ist3), &
           wfir2(:,:,ist2),wfir1(:,:,ist3),zrhomt,zrhoir)
          zt1=zfinp(.true.,zrhomt,zvclmt,zrhoir,zvclir)
! compute the density coefficient of the smallest G+q-vector
          call zrhogp(jlgq0r,ylmgq(:,igq0),sfacgq0,zrhomt,zrhoir,zrho01)
          zt2=cfq*wiq2(iq)*(conjg(zrho01)*zrho02)
          vnlijjk(ist3,ist1,ist2,ik)=wkptnr(ikp)*zt1+zt2
! end loop over ist3
        end if
      end do
! end loop over ist2
    end do
! end loop over ist1
  end do
! calculate the lower diagonal
  do ist1=1,nstsv
    do ist3=1,nstsv
      if (ist1.gt.ist3) then
        vnlijjk(ist3,ist1,:,ik)=conjg(vnlijjk(ist1,ist3,:,ik))
      end if
    end do
  end do
! copy the matrix elements of the type i-jj-i to vnlijji
  do ist1=1,nstsv
    do ist2=1,nstsv
      vnlijji(ist1,ist2,ik)=dble(vnlijjk(ist1,ist1,ist2,ik))
    end do
  end do
! end loop over reduced k-point set
end do
deallocate(igkignr,vgklnr,vgkcnr,gkcnr,tpgkcnr)
deallocate(vgqc,tpgqc,gqc,jlgqr,jlgq0r)
deallocate(evalsvp,evalsvnr)
deallocate(apwalm,evecfv,evecsv,sfacgknr,ylmgq,sfacgq)
deallocate(wfmt1,wfmt2,wfir1,wfir2)
deallocate(zrhomt,zrhoir,zvclmt,zvclir)
return
end subroutine
!EOC

