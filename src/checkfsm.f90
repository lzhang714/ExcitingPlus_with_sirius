
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine checkfsm
use modmain
implicit none
! local variables
integer is,ia,ja
integer isym,lspn
real(8) sc(3,3),v(3),t1
! external functions
real(8) r3taxi
external r3taxi
if (fixspin.eq.0) return
do isym=1,nsymcrys
  lspn=lspnsymc(isym)
  sc(:,:)=dble(symlatd(lspn))*symlatc(:,:,lspn)
! check invariance of global moment
  if ((abs(fixspin).eq.1).or.(abs(fixspin).eq.3)) then
    v(:)=sc(:,1)*momfix(1)+sc(:,2)*momfix(2)+sc(:,3)*momfix(3)
    t1=r3taxi(momfix,v)
    if (t1.gt.epslat) then
      write(*,*)
      write(*,'("Error(checkfsm): momfix not invariant under symmetry group")')
      write(*,*)
      stop
    end if
  end if
! check invariance of muffin-tin moments
  if ((abs(fixspin).eq.2).or.(abs(fixspin).eq.3)) then
    do is=1,nspecies
      do ia=1,natoms(is)
! equivalent atom
        ja=ieqatom(ia,is,isym)
        v(:)=sc(:,1)*mommtfix(1,ja,is)+sc(:,2)*mommtfix(2,ja,is) &
         +sc(:,3)*mommtfix(3,ja,is)
        t1=r3taxi(mommtfix(:,ia,is),v)
        if (t1.gt.epslat) then
          write(*,*)
          write(*,'("Error(checkfsm): mommtfix not invariant under symmetry&
           & group")')
          write(*,'(" for species ",I4)') is
          write(*,'(" and equivalent atoms ",2I4)') ia,ja
          write(*,*)
          stop
        end if
      end do
    end do
  end if
end do
return
end subroutine

