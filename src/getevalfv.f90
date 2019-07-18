
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getevalfv(vpl,evalfv_org)
use modmain
implicit none
! arguments
real(8), intent(in) :: vpl(3)
real(8), intent(out) :: evalfv_org(nstfv,nspnfv)
! local variables
integer isym,ik
integer recl,nstfv_,nspnfv_
real(8) vkl_(3),t1
! find the k-point number
call findkpt(vpl,isym,ik)
! find the record length
inquire(iolength=recl) vkl_,nstfv_,nspnfv_,evalfv_org
open(70,file=trim(scrpath)//'EVALFV'//trim(filext),action='READ', &
 form='UNFORMATTED',access='DIRECT',recl=recl)
read(70,rec=ik) vkl_,nstfv_,nspnfv_,evalfv_org
close(70)
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getevalfv): differing vectors for k-point ",I8)') ik
  write(*,'(" current    : ",3G18.10)') vkl(:,ik)
  write(*,'(" EVALFV.OUT : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (nstfv.ne.nstfv_) then
  write(*,*)
  write(*,'("Error(getevalfv): differing nstfv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nstfv
  write(*,'(" EVALFV.OUT : ",I8)') nstfv_
  write(*,*)
  stop
end if
if (nspnfv.ne.nspnfv_) then
  write(*,*)
  write(*,'("Error(getevalfv): differing nspnfv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nspnfv
  write(*,'(" EVALFV.OUT : ",I8)') nspnfv_
  write(*,*)
  stop
end if
return
end subroutine

