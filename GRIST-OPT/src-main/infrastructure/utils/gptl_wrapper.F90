#ifdef WITH_GPTL
  subroutine t_startf(name)
    implicit none
#include "gptl.inc"
    character(len=*) :: name
    integer :: ierr
    ierr = gptlstart(name)
  end subroutine t_startf
  subroutine t_stopf(name)
    implicit none
#include "gptl.inc"
    character(len=*) :: name
    integer :: ierr
    ierr = gptlstop(name)
  end subroutine t_stopf
  subroutine t_initf()
    implicit none
#include "gptl.inc"
    integer :: ierr
    ierr = gptlinitialize()
  end subroutine t_initf
  subroutine t_prf()
    use mpi
    implicit none
#include "gptl.inc"
    integer :: ierr, rank
    ierr = gptlpr_summary_file(MPI_COMM_WORLD, "grist_timing_stats")
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
    if (rank .eq. 0) then
      ierr = gptlpr_file("grist_timing.0")
    endif
  end subroutine
#else
  subroutine t_startf(name)
    implicit none
    character(len=*) :: name
    integer :: ierr
  end subroutine t_startf
  subroutine t_stopf(name)
    implicit none
    character(len=*) :: name
    integer :: ierr
  end subroutine t_stopf
  subroutine t_initf()
    implicit none
    integer :: ierr
  end subroutine t_initf
  subroutine t_prf()
    implicit none
    integer :: ierr
  end subroutine
#endif
