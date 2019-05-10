
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.


subroutine genffacg(is,ngv,ffacg)
use modmain
implicit none
! arguments
integer, intent(in) :: is
integer, intent(in) :: ngv
real(8), intent(out) :: ffacg(ngv)
! local variables
integer ig
real(8) t1,t2
t1=fourpi/omega
ffacg(1)=(t1/3.d0)*rmt(is)**3
do ig=2,ngv                     ! Long: btw, why the loop "ig=2,ngv" rather than "ig=1,ngv"?
  t2=gc(ig)*rmt(is)
          if (ig.eq.2)  write(*,*)' debug flag, genffacg, 1 '
  ffacg(ig)=t1*(sin(t2)-t2*cos(t2))/(gc(ig)**3)
          if (ig.eq.2)  write(*,*)' debug flag, genffacg, 2 '
end do
return
end subroutine


! --------------------------------------------------------------
! copy below the complete original EXCITING-PLUS genffacg.f90
! --------------------------------------------------------------
!subroutine genffacg(is,ngv,ffacg)
!use modmain
!implicit none
!! arguments
!integer, intent(in) :: is
!integer, intent(in) :: ngv
!real(8), intent(out) :: ffacg(ngv)
!! local variables
!integer ig
!real(8) t1,t2
!t1=fourpi/omega
!ffacg(1)=(t1/3.d0)*rmt(is)**3
!do ig=2,ngv                     ! Long: btw, why the loop "ig=2,ngv" rather than "ig=1,ngv"?
!  t2=gc(ig)*rmt(is)
!          if (ig .eq. 2)  write(*,*)' debug flag, genffacg, 1 '
!  ffacg(ig)=t1*(sin(t2)-t2*cos(t2))/(gc(ig)**3)
!          if (ig .eq. 2)  write(*,*)' debug flag, genffacg, 2 '
!end do
!return
!end subroutine


! --------------------------------------------------------------
! copy below the complete original EXCITING genffacg.f90
! --------------------------------------------------------------
!Subroutine genffacg (is, ffacg)
!      Use modmain
!      Use modinput
!      Implicit None
!! arguments
!      Integer, Intent (In) :: is
!      Real (8), Intent (Out) :: ffacg (ngvec)
!! local variables
!      Integer :: ig
!      Real (8) :: t1, t2, t3, t4
!      t1 = fourpi / omega
!      t2 = input%groundstate%cfdamp / input%groundstate%gmaxvr
!      Do ig = 1, ngvec
!         If (gc(ig) .Gt. input%structure%epslat) Then
!            If (input%groundstate%cfdamp .Ne. 0.d0) Then
!! use damping if required
!               t3 = Exp (-(t2*gc(ig))**2)
!            Else
!               t3 = 1.d0
!            End If
!            t4 = gc (ig) * rmt (is)
!            ffacg (ig) = t1 * t3 * (Sin(t4)-t4*Cos(t4)) / (gc(ig)**3)
!         Else
!            ffacg (ig) = (t1/3.d0) * rmt (is) ** 3
!         End If
!      End Do
!      Return
!End Subroutine


! --------------------------------------------------------------
! copy below the complete original ELK 5.2.14  genffacg.f90
! --------------------------------------------------------------
!subroutine genffacgp(is,gpc,ffacgp)
!use modmain
!implicit none
!! arguments
!integer, intent(in) :: is
!real(8), intent(in) :: gpc(ngtot)
!real(8), intent(out) :: ffacgp(ngtot)
!! local variables
!integer ig
!real(8) t1,t2
!t1=fourpi/omega
!do ig=1,ngtot
!  if (gpc(ig).gt.epslat) then
!    t2=gpc(ig)*rmt(is)
!    ffacgp(ig)=t1*(sin(t2)-t2*cos(t2))/(gpc(ig)**3)
!  else
!    ffacgp(ig)=(t1/3.d0)*rmt(is)**3
!  end if
!end do
!return
!end subroutine


