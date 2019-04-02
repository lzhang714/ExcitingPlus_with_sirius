
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.


!BOP
! !ROUTINE: rdmwritededn
! !INTERFACE:
subroutine rdmwritededn(dedn)
! !USES:
use modrdm
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   dedn : derivative of energy (in,real(nstsv,nkpt))
! !DESCRIPTION:
!   Writes the derivative of total energy w.r.t. occupation numbers to file
!   {\tt RDM\_DEDN.OUT}.
!
! !REVISION HISTORY:
!   Created 2008 (Sharma)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: dedn(nstsv,nkpt)
! local variables
integer ik,ist
open(93,file='RDM_DEDN.OUT',action='WRITE',form='FORMATTED')
write(93,'(I6," : nkpt")') nkpt
write(93,'(I6," : nstsv")') nstsv
do ik=1,nkpt
  write(93,*)
  write(93,'(I6,3G18.10," : k-point, vkl")') ik,vkl(:,ik)
  write(93,'("     (state, occupancy and derivative below)")')
  do ist=1,nstsv
    write(93,'(I6,4G18.10)') ist,occsv(ist,ik),-dedn(ist,ik)
  end do
end do
close(93)
return
end subroutine
!EOC
