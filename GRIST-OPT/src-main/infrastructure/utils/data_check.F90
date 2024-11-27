subroutine datacheckscalari4(dat, ref, name, file, line)
  use grist_mpi
  implicit none
  integer(kind=4)     :: dat, ref
  character(len=*) :: name, file
  integer          :: line
  integer          :: i
  if (dat /= ref) then
    print *, "proc=", mpi_rank(), "check of ", name, " failed in ", file, ":", line, dat, ref, "diff:", ABS(dat - ref)
    stop
  else
    if (mpi_rank() == 0) &
      print *, "Check of ", name, " success in ", file, ":", line
  endif
end subroutine datacheckscalari4


subroutine datacheckscalar(dat, ref, name, file, line)
  use grist_mpi
  implicit none
  real(kind=8)     :: dat, ref
  character(len=*) :: name, file
  integer          :: line
  integer          :: i
  if (dat /= ref) then
    print *, "proc=", mpi_rank(), "check of ", name, " failed in ", file, ":", line, dat, ref, "diff:", ABS(dat - ref)
    stop
  else
    if (mpi_rank() == 0) &
      print *, "Check of ", name, " success in ", file, ":", line
  endif
end subroutine datacheckscalar

subroutine datacheck1d(dat, ref, name, file, line)
  use grist_mpi
  implicit none
  real(kind=8)     :: dat(:), ref(:)
  character(len=*) :: name, file
  integer          :: line
  integer          :: i
  if (ANY((dat .ne. ref) .and. .not. (ref .ne. ref))) then
    print *, "proc=", mpi_rank(), "check of ", name, " failed in ", file, ":", line, "max diff:", MAXVAL(ABS(dat - ref))
    do i = lbound(dat, 1), ubound(dat,1)
      if (dat(i) /= ref(i)) then
        print *, "proc=", mpi_rank(), "first err at:", i, "of (", lbound(dat, 1), ":", ubound(dat,1), ")", "dat=", dat(i), "ref=", ref(i)
        stop
      end if
    end do
  else
    if (mpi_rank() == 0) &
      print *, "Check of ", name, " success in ", file, ":", line
  endif
end subroutine datacheck1d

subroutine datacheck2d(dat, ref, name, file, line)
  use grist_mpi
  implicit none
  real(kind=8)     :: dat(:,:), ref(:,:)
  character(len=*) :: name, file
  integer          :: line
  integer          :: i, j
  if (ANY((dat .ne. ref) .and. .not. (ref .ne. ref))) then
    print *, "proc=", mpi_rank(), "check of ", name, " failed in ", file, ":", line, "max diff:", MAXVAL(ABS(dat - ref))
    do i = lbound(dat, 2), ubound(dat,2)
      do j = lbound(dat, 1), ubound(dat, 1)
        if (dat(j,i) /= ref(j, i)) then
          print *, "proc=", mpi_rank(), "first err at:", j, i, "of (", lbound(dat, 1), ":", ubound(dat,1), ",", lbound(dat, 2), ":", ubound(dat,2), ")", "dat=", dat(j,i), "ref=", ref(j,i)
          stop
        end if
      end do
    end do
  else
    if (mpi_rank() == 0) &
      print *, "Check of ", name, " success in ", file, ":", line
  endif
end subroutine datacheck2d

subroutine datacheck3d(dat, ref, name, file, line)
  use grist_mpi
  implicit none
  real(kind=8)     :: dat(:,:,:), ref(:,:,:)
  character(len=*) :: name, file
  integer          :: line
  integer          :: i, j, k
  if (ANY((dat .ne. ref) .and. .not. (ref .ne. ref))) then
    print *, "proc=", mpi_rank(), "check of ", name, " failed in ", file, ":", line, "max diff:", MAXVAL(ABS(dat - ref))
    do i = lbound(dat, 3), ubound(dat, 3)
      do j = lbound(dat, 2), ubound(dat,2)
        do k = lbound(dat, 1), ubound(dat, 1)
          if (dat(k,j,i) /= ref(k,j,i)) then
            print *, "first err at:", k, j, i, "dat=", dat(k,j,i), "ref=", ref(k,j,i)
            stop
          end if
        end do
      end do
    end do
    print *, "Check of ", name, " failed in ", file, ":", line, "max diff:", MAXVAL(ABS(dat - ref)), 1e-15
  else
    if (mpi_rank() == 0) &
      print *, "Check of ", name, " success in ", file, ":", line
  endif
end subroutine datacheck3d
