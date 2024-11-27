!================================================
!  Created by zhangyi on 16/8/5.
!================================================

 module grist_handle_error

   implicit none

   private

   public :: endrun, &
             showinfo
 
   character(len=200), public :: grist_message

   contains

   subroutine endrun(info, number)
!
! io
!
   character*(*), optional, intent(in)   :: info
   integer, optional, intent(in)         :: number

     if(present(info).and.present(number))then
        print*,"--------------------------------------------"
        print*,"                                            "
        print*,"          ENDRUN: ",trim(info), number
        print*,"                                            "
        print*,"--------------------------------------------"
     end if

     stop
   end subroutine endrun

   subroutine showinfo(info, number)
!
! io
!
   character*(*), optional, intent(in)   :: info
   integer, optional, intent(in)         :: number

     if(present(info).and.present(number))then
        print*,"--------------------------------------------"
        print*,"                                            "
        print*,"          SHOWINFO: ",trim(info), number
        print*,"                                            "
        print*,"--------------------------------------------"
     end if

   end subroutine showinfo

 end module grist_handle_error
