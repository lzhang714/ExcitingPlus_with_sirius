
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine exxengy
use modmain
use modqpt
implicit none
! local variables
integer is,ia,nrc,m1,m2
integer ik,ist,jst
real(8) evv,ecv,ecc
complex(8) zt1
! allocatable arrays
complex(8), allocatable :: wfcr1(:,:,:)
complex(8), allocatable :: wfcr2(:,:,:)
complex(8), allocatable :: zrhomt(:,:)
complex(8), allocatable :: zvclmt(:,:)
complex(8), allocatable :: zfmt(:,:)
! external functions
complex(8) zfmtinp
external zfmtinp
allocate(wfcr1(lmmaxvr,nrcmtmax,2))
allocate(wfcr2(lmmaxvr,nrcmtmax,2))
allocate(zrhomt(lmmaxvr,nrcmtmax))
allocate(zvclmt(lmmaxvr,nrcmtmax))
allocate(zfmt(lmmaxvr,nrcmtmax))
evv=0.d0
ecv=0.d0
ecc=0.d0
do ik=1,nkpt
  write(*,'("Info(exxengy): ",I6," of ",I6," k-points")') ik,nkpt
  call exxengyk(ik,evv,ecv)
end do
!-----------------------------------!
!    core-core-core contribution    !
!-----------------------------------!
! begin loops over atoms and species
do is=1,nspecies
  nrc=nrcmt(is)
  do ia=1,natoms(is)
    do jst=1,spnst(is)
      if (spcore(jst,is)) then
        do m2=-spk(jst,is),spk(jst,is)-1
          call wavefcr(lradstp,is,ia,jst,m2,nrcmtmax,wfcr2)
          do ist=1,spnst(is)
            if (spcore(ist,is)) then
              do m1=-spk(ist,is),spk(ist,is)-1
                call wavefcr(lradstp,is,ia,ist,m1,nrcmtmax,wfcr1)
! calculate the complex overlap density
                call vnlrhomt(.true.,is,wfcr1(:,:,1),wfcr2(:,:,1),zrhomt)
                call vnlrhomt(.true.,is,wfcr1(:,:,2),wfcr2(:,:,2),zfmt)
                zrhomt(:,1:nrc)=zrhomt(:,1:nrc)+zfmt(:,1:nrc)
! calculate the Coulomb potential
                call zpotclmt(ptnucl,lmaxvr,nrc,rcmt(:,is),0.d0,lmmaxvr, &
                 zrhomt,zvclmt)
                zt1=zfmtinp(.true.,lmaxvr,nrc,rcmt(:,is),lmmaxvr,zrhomt,zvclmt)
                ecc=ecc-0.5d0*dble(zt1)
              end do
! end loop over ist
            end if
          end do
        end do
! end loop over jst
      end if
    end do
! end loops over atoms and species
  end do
end do
! total exchange energy
engyx=evv+ecv+ecc
deallocate(wfcr1,wfcr2,zrhomt,zvclmt,zfmt)
return
end subroutine

