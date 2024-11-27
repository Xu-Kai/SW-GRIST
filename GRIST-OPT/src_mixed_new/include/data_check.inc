interface data_check_intf
  subroutine datacheckscalari4(dat, ref, name, file, line)
    implicit none
    integer(kind=4)     :: dat, ref
    character(len=*) :: name, file
    integer          :: line
  end subroutine datacheckscalari4


  subroutine datacheckscalar(dat, ref, name, file, line)
    implicit none
    real(kind=8)     :: dat, ref
    character(len=*) :: name, file
    integer          :: line
  end subroutine datacheckscalar

  subroutine datacheck1d(dat, ref, name, file, line)
    implicit none
    real(kind=8)     :: dat(:), ref(:)
    character(len=*) :: name, file
    integer          :: line
  end subroutine datacheck1d

  subroutine datacheck2d(dat, ref, name, file, line)
    implicit none
    real(kind=8)     :: dat(:,:), ref(:,:)
    character(len=*) :: name, file
    integer          :: line
  end subroutine datacheck2d
  subroutine datacheck3d(dat, ref, name, file, line)
    implicit none
    real(kind=8)     :: dat(:,:,:), ref(:,:,:)
    character(len=*) :: name, file
    integer          :: line
  end subroutine datacheck3d

end interface data_check_intf
#ifndef data_check
#define data_check(dat, ref) data_check_intf(dat, ref, "dat", __FILE__, __LINE__)
#endif