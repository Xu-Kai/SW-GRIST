
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description: Procedures related with hydrostatic equation, pgf, and terms 
!              contain hpressure like exner, delhp, etc..
!              They only handle vertical 1D data structure, should not be aware
!              of the horizontal data structure; except calc_hpe_get_full_mass
! Revision history:
!----------------------------------------------------------------------------

  module grist_hpe_hydro_pgf

    use grist_constants,      only: i4, r8, zero, rdry,cp, p00
    use grist_hpe_constants,  only: p0, nlev, nlevp, eta_face_a, eta_face_b, eta_full_a,eta_full_b

    implicit none

    private

    public :: calc_hpe_hydro               , &
              calc_hpe_hpressure_face_level, &
              calc_hpe_hpressure_full_level, &
              calc_hpe_delhp               , &
              calc_hpe_exner               , &
              calc_hpe_get_full_mass

   contains
!
!  vertical index
!-(1)--------------------------------- 1/2
!  ..............[1]..................
!-(2)--------------------------------- 1+1/2
!  ..............[2]..................
!-(3)--------------------------------- 2+1/2
!  ..............[3]..................
!-(4)--------------------------------- 3+1/2
!  ..............[4]..................
!-(5)--------------------------------- 4+1/2
!  ..............[5]..................
!-(6)------------------ilevel-1------- 5+1/2
!  ..............[6]...ilevel-1.......
!-(7)------------------ilevel--------- 6+1/2
!  ..............[7]...ilevel.........
!-(8)------------------ilevel+1------- 7+1/2
!  ..............[8]...ilevel+1.......
!-(9)------------------ilevel+2------- 8+1/2
!
!
    subroutine calc_hpe_hydro(scalar_temperature_full_level  , &
                              scalar_hpressure_face_level    , &
                              scalar_delhp_full_level        , &
                              scalar_geopotential_surface    , &
                              scalar_geopotential_face_level , &
                              scalar_geopotential_full_level )
!
! io
!
      real(r8), dimension(nlev),    intent(in)    :: scalar_temperature_full_level
      real(r8), dimension(nlevp),   intent(in)    :: scalar_hpressure_face_level
      real(r8), dimension(nlev),    intent(in)    :: scalar_delhp_full_level
      real(r8),                     intent(in)    :: scalar_geopotential_surface
      real(r8), dimension(nlevp),   intent(inout) :: scalar_geopotential_face_level
      real(r8), dimension(nlev),    intent(inout) :: scalar_geopotential_full_level
!
! local
!
      real(r8), dimension(nlev)            :: alpha
      real(r8)                             :: vert_sum 
      integer(i4)                          :: ilev
      integer(i4)                          :: ilevel
!$omp declare target

!      if(.not.allocated(alpha)) allocate(alpha(nlev))

!      call calc_hpe_alpha_sb(scalar_hpressure_face_level, scalar_delhp_full_level, alpha)
!
! 3.5 of SS81
!
      do ilevel = 1, nlev     ! for each face ilevel except surface
         vert_sum = 0._r8        
         do ilev = ilevel, nlev ! for each full level ilev
            vert_sum = vert_sum+rdry*scalar_temperature_full_level(ilev)*&
                                 log(scalar_hpressure_face_level(ilev+1)/&
                                     scalar_hpressure_face_level(ilev))
         end do
         scalar_geopotential_face_level(ilevel) = scalar_geopotential_surface+vert_sum
      end do

         scalar_geopotential_face_level(nlev+1) = scalar_geopotential_surface
!
! 3.6 of SS81
!
!      do ilev = 1, nlev   ! for each full level
!         scalar_geopotential_full_level(ilev) = scalar_geopotential_face_level(ilev+1)+&
!                                                  alpha(ilev)*rdry*&
!                                                  scalar_temperature_full_level(ilev)
!      end do
!
! original hdc use above, comment below for bit reproduce
! 
      do ilev = 1, nlev   ! for each full level
         scalar_geopotential_full_level(ilev) = 0.5_r8*(scalar_geopotential_face_level(ilev+1)+&
                                                        scalar_geopotential_face_level(ilev))
      end do

!      deallocate(alpha)
      
      return
    end subroutine calc_hpe_hydro

    subroutine calc_hpe_hpressure_face_level(ps,scalar_hpressure_face_level)
!$omp declare target
      real(r8),                   intent(in)     :: ps
      real(r8), dimension(nlevp), intent(inout)  :: scalar_hpressure_face_level
      integer(i4)                        :: ilev

      do ilev = 1, nlev+1
         scalar_hpressure_face_level(ilev) = eta_face_a(ilev)*p0+eta_face_b(ilev)*ps
      end do

      return
    end subroutine calc_hpe_hpressure_face_level

    subroutine calc_hpe_delhp(scalar_hpressure_face_level,scalar_delhp_full_level)
!$omp declare target
      real(r8), dimension(nlevp),  intent(in)    :: scalar_hpressure_face_level
      real(r8), dimension(nlev),  intent(inout)  :: scalar_delhp_full_level
      integer(i4)   :: ilev

        do ilev = 1, nlev     ! for each full level
           scalar_delhp_full_level(ilev) = scalar_hpressure_face_level(ilev+1)-&
                                             scalar_hpressure_face_level(ilev) 
        end do
      return
    end subroutine calc_hpe_delhp

    subroutine calc_hpe_hpressure_full_level(ps, &
                                            scalar_hpressure_face_level,&
                                            scalar_delhp_full_level    ,&
                                            scalar_hpressure_full_level)
!$omp declare target
      real(r8),                  intent(in)     ::  ps
      real(r8), dimension(nlevp),intent(in)     ::  scalar_hpressure_face_level
      real(r8), dimension(nlev), intent(in)     ::  scalar_delhp_full_level
      real(r8), dimension(nlev), intent(inout)  ::  scalar_hpressure_full_level
!
! local
!
      integer(i4)                        ::  ilev
!
! 3.17 of SS81
!
       scalar_hpressure_full_level(1) = 0.5_r8*scalar_delhp_full_level(1)

       do ilev = 2, nlev
          scalar_hpressure_full_level(ilev) = scalar_delhp_full_level(ilev)/&
                                           log(scalar_hpressure_face_level(ilev+1)/&
                                               scalar_hpressure_face_level(ilev))
       end do
!
! original hdc version use above, comment below for bit reproduce of hdc
!
      do ilev = 1, nlev
         scalar_hpressure_full_level(ilev) = eta_full_a(ilev)*p0+eta_full_b(ilev)*ps
      end do

      return
    end subroutine calc_hpe_hpressure_full_level

    subroutine calc_hpe_exner(scalar_pressure_full_level,&
                              scalar_exner_full_level)

      real(r8), dimension(nlev),  intent(in)     :: scalar_pressure_full_level
      real(r8), dimension(nlev),  intent(inout)  ::    scalar_exner_full_level
      integer  :: ilev

       do ilev = 1, nlev
          scalar_exner_full_level(ilev) = (scalar_pressure_full_level(ilev)/p00)**(rdry/cp)
       end do

       return
    end subroutine calc_hpe_exner

    subroutine calc_hpe_alpha_sb(scalar_hpressure_face_level, &
                                 scalar_delhp_full_level    , &
                                 alpha)

      real(r8), dimension(nlevp), intent(in)     :: scalar_hpressure_face_level
      real(r8), dimension(nlev),  intent(in)     :: scalar_delhp_full_level
      real(r8), dimension(nlev),  intent(inout)  :: alpha
      integer(i4)                         :: ilev
!
! 3.9 of SS81
!
      do ilev = 1, nlev
         alpha(ilev) = 1._r8-(scalar_hpressure_face_level(ilev)/scalar_delhp_full_level(ilev))*&
                           log(scalar_hpressure_face_level(ilev+1)/scalar_hpressure_face_level(ilev))
      end do

      return
    end subroutine calc_hpe_alpha_sb

    subroutine calc_hpe_get_full_mass(nlev, ncell, working_mode          , & ! in
                                          scalar_hpressure_at_pc_face_level  , & ! in
                                          scalar_hpressure_at_pc_full_level  , & ! in
                                          scalar_delhp_at_pc_full_level      , & ! in
                                          scalar_mif_at_pc_full_level        , & ! in
                                          scalar_pressure_at_pc_full_level   , & ! out
                                          scalar_pressure_at_pc_face_level)      ! out
! io
     integer(i4)        , intent(in)    ::  nlev, ncell
     character(len=*),    intent(in)    ::  working_mode
     real(r8),            intent(in)    ::  scalar_hpressure_at_pc_face_level(:,:) 
     real(r8),            intent(in)    ::  scalar_hpressure_at_pc_full_level(:,:) 
     real(r8),            intent(in)    ::  scalar_delhp_at_pc_full_level(:,:)
     real(r8),            intent(in)    ::  scalar_mif_at_pc_full_level(:,:)
     real(r8),            intent(out)   ::  scalar_pressure_at_pc_full_level(:,:)
     real(r8),            intent(out)   ::  scalar_pressure_at_pc_face_level(:,:)
! local
     integer(i4)  :: iv, ilev, ilevel

#if defined(CHK_HPE_GET_FULL_MASS) || defined(CHK_ALL)
#include "data_check.inc"
      real(r8), allocatable :: ref_scalar_pressure_at_pc_face_level(:,:)
      
#endif
      integer(i4) itcell
     select case(trim(working_mode))

     case('dycore') ! for reference
call t_startf("hpe_get_full_mass_copy")
!$omp target
!$omp parallel workshare
       scalar_pressure_at_pc_face_level(:,:)  =  scalar_hpressure_at_pc_face_level(:,:)
!$omp end parallel workshare
!$omp end target
!$omp target
!$omp parallel workshare
       scalar_pressure_at_pc_full_level(:, :)  =  scalar_hpressure_at_pc_full_level(:,:)
!$omp end parallel workshare
!$omp end target 
call t_stopf("hpe_get_full_mass_copy")
       return
     case('tracer') ! for reference
       return
     case default ! this should produce very similar  pressure as hpressure in a dry mode;
                 ! but not bit reproducity due to machine round off;
                 ! at the first top, assume p=hp

      ! scalar_pressure_at_pc_face_level(1,1:ncell) = scalar_hpressure_at_pc_face_level(1,1:ncell)
#if defined(CHK_HPE_GET_FULL_MASS) || defined(CHK_ALL)
      allocate(ref_scalar_pressure_at_pc_face_level, source=scalar_pressure_at_pc_face_level)
#endif
call t_startf("hpe_get_full_mass_target")

!$omp target
!$omp parallel private(itcell, ilev)
!$omp do
      do itcell = 1, ncell
            scalar_pressure_at_pc_face_level(1,itcell) = scalar_hpressure_at_pc_face_level(1,itcell)
       do ilev = 2, nlev+1 ! for each full level ilev
          scalar_pressure_at_pc_face_level(ilev,itcell) = scalar_pressure_at_pc_face_level(ilev-1,itcell)+&
                                                     scalar_delhp_at_pc_full_level(ilev-1,itcell)/scalar_mif_at_pc_full_level(ilev-1,itcell)
       end do
      enddo
!$omp end do nowait
!$omp end parallel 
!$omp end target
call t_stopf("hpe_get_full_mass_target")

#if defined(CHK_HPE_GET_FULL_MASS) || defined(CHK_ALL)

call t_startf("hpe_get_full_mass_target_1")
      do ilev = 2, nlev+1 ! for each full level ilev
            ref_scalar_pressure_at_pc_face_level(ilev,1:ncell) = ref_scalar_pressure_at_pc_face_level(ilev-1,1:ncell)+&
                                                     scalar_delhp_at_pc_full_level(ilev-1,1:ncell)/scalar_mif_at_pc_full_level(ilev-1,1:ncell)
       end do
call t_stopf("hpe_get_full_mass_target_1")
call data_check(ref_scalar_pressure_at_pc_face_level, scalar_pressure_at_pc_face_level)
deallocate(ref_scalar_pressure_at_pc_face_level)
#endif
call t_startf("hpe_get_full_mass_copy_2")
!$omp target
!$omp parallel workshare
       scalar_pressure_at_pc_full_level(1:nlev,1:ncell)  = 0.5_r8*(scalar_pressure_at_pc_face_level(1:nlev,1:ncell)+scalar_pressure_at_pc_face_level(2:nlev+1,1:ncell))
!$omp end parallel workshare
!$omp end target
call t_stopf("hpe_get_full_mass_copy_2")
     end select
    
     return
   end subroutine calc_hpe_get_full_mass

 end module grist_hpe_hydro_pgf
