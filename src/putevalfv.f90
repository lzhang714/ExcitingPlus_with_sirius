
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putevalfv(ik,evalfv)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
real(8), intent(in) :: evalfv(nstfv,nspnfv)
! local variables
integer recl
! find the record length
inquire(iolength=recl) vkl(:,ik),nstfv,nspnfv,evalfv
open(70,file=trim(scrpath)//'EVALFV'//trim(filext),action='WRITE', &
 form='UNFORMATTED',access='DIRECT',recl=recl)
write(70,rec=ik) vkl(:,ik),nstfv,nspnfv,evalfv
close(70)
return
end subroutine


!!subroutine putevalfv_EP_vs_SIRIUS(ik,evalfv_EP,evalfv_SR)
!!use modmain
!!implicit none
!!! arguments
!!integer, intent(in) :: ik
!!real(8), intent(in) :: evalfv_EP(nstfv,nspnfv)
!!real(8), intent(in) :: evalfv_SR(nstfv,nspnfv)
!!! local variables
!!integer recl
!!! find the record length
!!inquire(iolength=recl) vkl(:,ik),nstfv,nspnfv,evalfv_EP,evalfv_SR
!!open(70,file=trim(scrpath)//'EVALFV_EP_vs_SIRIUS'//trim(filext),action='WRITE', &
!! form='UNFORMATTED',access='DIRECT',recl=recl)
!!write(70,rec=ik) vkl(:,ik),nstfv,nspnfv,evalfv_EP,evalfv_SR
!!close(70)
!!return
!!end subroutine
