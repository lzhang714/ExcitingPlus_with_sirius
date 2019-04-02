subroutine initmpigrid
use modmain
use mod_addons_q
use mod_wannier
implicit none
integer, allocatable :: d(:)
integer i1,i2,nd,i
select case(task)
  case(0,1,2,3,20,21,805,822,700,701,702)
    nd=2
    allocate(d(nd)); d(:)=1
    if (nproc.le.nkpt) then
      d(dim1)=nproc
    else
      d(dim1)=nkpt
      d(dim2)=nproc/nkpt
    endif
! response code also runs on a 2D grid but the total number of k-points is distributed 
  case(800,801,802,8022)
    nd=2
    allocate(d(nd)); d(:)=1
    if (nproc.le.nkptnr) then
      d(dim1)=nproc
    else  
      d(dim1)=nkptnr
      d(dim2)=nproc/nkptnr
    endif
  case(863)
    nd=2
    allocate(d(nd)); d(:)=1
    if (nproc.le.nrxyz(1)) then
       d(dim2)=nproc
    else
       d(dim2)=nrxyz(1)
       d(dim1)=nproc/nrxyz(1)
    endif
  case default
    nd=1
    allocate(d(nd))
    d=nproc
end select  
! overwrite default grid layout
if (lmpigrid) then
  deallocate(d)
  allocate(d(mpigrid_ndim))
  d(1:mpigrid_ndim)=mpigrid(1:mpigrid_ndim)
endif
call mpi_grid_initialize(d)
if (mpi_grid_root()) then
  write(*,*)
  write(*,'("[initmpigrid] mpi grid size : ",10I8)')(mpi_grid_dim_size(i),i=1,mpi_grid_nd)
endif
deallocate(d)
return
end
