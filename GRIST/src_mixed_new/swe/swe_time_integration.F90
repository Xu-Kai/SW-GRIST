
!----------------------------------------------------------------------------
! Created on 2016
! Version 1.0
! Description: Time integration schemes, only retain rk23 and fit now
! Revision history: 
!----------------------------------------------------------------------------

 module swe_time_integration

  use grist_constants,       only: i4, r8
  use grist_domain_types,    only: global_domain
  use grist_data_types,      only: scalar_1d_field,exchange_field_list_1d
  use swe_spatial_module,    only: calc_spatial_tend
  use swe_spatial_module_tr, only: calc_spatial_tend_tr
  use grist_nml_module,      only: time_scheme,use_limiter_rk3, &
                                   isolate_advection_test
  use grist_config_partition,only: exchange_data_1d_add,exchange_data_1d

  implicit none

  private

  public :: time_integration_rk       ,&
            time_integration_rk3_tr   ,& ! only integrate advection
            time_integration_fit_tr   ,& ! only integrate advection
            time_integration_fit         ! only for testing advection in isolation

  contains

!================================================
!                 Unified RK3&2
!================================================

  subroutine time_integration_rk(mesh,scalar_normal_velocity_at_edge     ,&
                                      scalar_height_at_prime_cell        ,&
                                      scalar_topo_at_prime_cell          ,&
                                      scalar_normal_velocity_at_edge_next,&
                                      scalar_height_at_prime_cell_next   ,&
                                      tend_normal_velocity_at_edge       ,&
                                      tend_height_at_prime_cell          ,&
                                      itimestep,dtime)

! io
  type(global_domain), intent(inout) :: mesh
  type(scalar_1d_field),   intent(in)    :: scalar_normal_velocity_at_edge
  type(scalar_1d_field),   intent(in)    :: scalar_height_at_prime_cell
  type(scalar_1d_field),   intent(in)    :: scalar_topo_at_prime_cell
  type(scalar_1d_field),   intent(inout) :: scalar_normal_velocity_at_edge_next
  type(scalar_1d_field),   intent(inout) :: scalar_height_at_prime_cell_next
  type(scalar_1d_field),   intent(inout) :: tend_normal_velocity_at_edge
  type(scalar_1d_field),   intent(inout) :: tend_height_at_prime_cell
  integer(i4)       ,   intent(in)    :: itimestep
  real(r8)          ,   intent(in)    :: dtime
! local
  type(scalar_1d_field)                  :: scalar_normal_velocity_tmpnew
  type(scalar_1d_field)                  :: scalar_normal_velocity_tend
  type(scalar_1d_field)                  :: scalar_height_tmpnew
  type(scalar_1d_field)                  :: scalar_height_tend
  integer(i4)                         :: irk_step
  integer(i4)                         :: nrk
  logical                             :: flag(3)
  type(exchange_field_list_1d),pointer:: field_head

  field_head => null()

     if(trim(time_scheme).eq.'rk3')then
        nrk     = 3
        flag(1) = .false.
        flag(2) = .false.
        flag(3) = use_limiter_rk3
     end if
     if(trim(time_scheme).eq.'rk2')then
        nrk     = 2
        flag(1) = .false.
        flag(2) = use_limiter_rk3 
        flag(3) = .false.
     end if
     if(trim(time_scheme).eq.'fit')then
        nrk     = 1
        flag(1) = use_limiter_rk3 
        flag(2) = .false.
        flag(3) = .false.
     end if

!================================================
!    Just initialization, will be overwritten
!================================================

    scalar_normal_velocity_tmpnew = scalar_normal_velocity_at_edge
    scalar_normal_velocity_tend   = scalar_normal_velocity_at_edge
    scalar_height_tmpnew          = scalar_height_at_prime_cell
    scalar_height_tend            = scalar_height_at_prime_cell

!================================================
!                 1st call
!================================================
     do irk_step = 1, nrk
  
          call calc_spatial_tend(mesh,scalar_normal_velocity_tmpnew,&
                                      scalar_height_tmpnew         ,&
                                      scalar_topo_at_prime_cell    ,&
                                      scalar_normal_velocity_tend  ,&
                                      scalar_height_tend           ,&
                                      itimestep, dtime, flag(irk_step),irk_step)

        scalar_normal_velocity_tmpnew%f(:) = scalar_normal_velocity_at_edge%f(:)+&
                                             dtime*scalar_normal_velocity_tend%f(:)/(nrk-irk_step+1)

        scalar_height_tmpnew%f(:)          = scalar_height_at_prime_cell%f(:)+&
                                             dtime*scalar_height_tend%f(:)/(nrk-irk_step+1)

        call exchange_data_1d_add(mesh,field_head,scalar_height_tmpnew)
        call exchange_data_1d_add(mesh,field_head,scalar_normal_velocity_tmpnew)
        call exchange_data_1d(mesh%local_block,field_head)

     end do

!================================================
!               Final Update
!================================================

     scalar_normal_velocity_at_edge_next%f(:) = scalar_normal_velocity_tmpnew%f(:)
     scalar_height_at_prime_cell_next%f(:)    = scalar_height_tmpnew%f(:)
     tend_normal_velocity_at_edge%f(:)        = scalar_normal_velocity_tend%f(:)
     tend_height_at_prime_cell%f(:)           = scalar_height_tend%f(:)

     if(isolate_advection_test) scalar_normal_velocity_at_edge_next%f(:) = 0._r8

    return
  end subroutine time_integration_rk


!================================================
!          RK3 at dual cell (WS2002)
!================================================

  subroutine time_integration_rk3_tr(mesh,scalar_normal_velocity_at_tr_edge     ,&
                                          scalar_height_at_dual_cell            ,&
                                          scalar_normal_velocity_at_tr_edge_next,&
                                          scalar_height_at_dual_cell_next       ,&
                                          itimestep, dtime)
! io
  type(global_domain), intent(inout) :: mesh
  type(scalar_1d_field),   intent(in)    :: scalar_normal_velocity_at_tr_edge
  type(scalar_1d_field),   intent(in)    :: scalar_height_at_dual_cell
  type(scalar_1d_field),   intent(inout) :: scalar_normal_velocity_at_tr_edge_next
  type(scalar_1d_field),   intent(inout) :: scalar_height_at_dual_cell_next
  integer(i4)       ,   intent(in)    :: itimestep
  real(r8)          ,   intent(in)    :: dtime
! local
  type(scalar_1d_field)      :: scalar_normal_velocity_tmpnew
  type(scalar_1d_field)      :: scalar_normal_velocity_tend1
  type(scalar_1d_field)      :: scalar_normal_velocity_tend2
  type(scalar_1d_field)      :: scalar_normal_velocity_tend3
  type(scalar_1d_field)      :: scalar_height_tmpnew
  type(scalar_1d_field)      :: scalar_height_tend1
  type(scalar_1d_field)      :: scalar_height_tend2
  type(scalar_1d_field)      :: scalar_height_tend3

!================================================
!    Just initialization, will be overwritten
!================================================

  scalar_normal_velocity_tmpnew = scalar_normal_velocity_at_tr_edge
  scalar_normal_velocity_tend1  = scalar_normal_velocity_at_tr_edge
  scalar_normal_velocity_tend2  = scalar_normal_velocity_at_tr_edge
  scalar_normal_velocity_tend3  = scalar_normal_velocity_at_tr_edge

  scalar_height_tmpnew          = scalar_height_at_dual_cell
  scalar_height_tend1           = scalar_height_at_dual_cell
  scalar_height_tend2           = scalar_height_at_dual_cell
  scalar_height_tend3           = scalar_height_at_dual_cell

!================================================
!                 1st call
!================================================

     call calc_spatial_tend_tr(mesh                             ,&
                               scalar_normal_velocity_at_tr_edge,&
                               scalar_height_at_dual_cell       ,&
                               scalar_height_tend1              ,&
                               itimestep                        ,&
                               dtime, .false., .true., 1)

     scalar_height_tmpnew%f(:)          = scalar_height_at_dual_cell%f(:)+&
                                          dtime*scalar_height_tend1%f(:)/3._r8

!================================================
!                 2nd call
!================================================
  
     call calc_spatial_tend_tr(mesh                         ,&
                               scalar_normal_velocity_tmpnew,&
                               scalar_height_tmpnew         ,&
                               scalar_height_tend2          ,&
                               itimestep                    ,& 
                               dtime, .false., .true., 2)

     scalar_height_tmpnew%f(:)          = scalar_height_at_dual_cell%f(:)+&
                                          dtime*scalar_height_tend2%f(:)/2._r8

!================================================
!                 3rd call
!================================================

     call calc_spatial_tend_tr(mesh                         ,&
                               scalar_normal_velocity_tmpnew,&
                               scalar_height_tmpnew         ,&
                               scalar_height_tend3          ,&
                               itimestep                    ,&
                               dtime, use_limiter_rk3, .true., 3)


!================================================
!               Final Update
!================================================

     scalar_height_at_dual_cell_next%f(:)   = scalar_height_at_dual_cell%f(:)+&
                                              dtime*scalar_height_tend3%f(:)

     return
  end subroutine time_integration_rk3_tr

  subroutine time_integration_fit_tr(mesh,scalar_normal_velocity_at_tr_edge     ,&
                                          scalar_height_at_dual_cell            ,&
                                          scalar_normal_velocity_at_tr_edge_next,&
                                          scalar_height_at_dual_cell_next       ,&
                                          itimestep, dtime)
! io
  type(global_domain), intent(inout) :: mesh
  type(scalar_1d_field),   intent(in)    :: scalar_normal_velocity_at_tr_edge
  type(scalar_1d_field),   intent(in)    :: scalar_height_at_dual_cell
  type(scalar_1d_field),   intent(inout) :: scalar_normal_velocity_at_tr_edge_next
  type(scalar_1d_field),   intent(inout) :: scalar_height_at_dual_cell_next
  integer(i4)       ,   intent(in)    :: itimestep
  real(r8)          ,   intent(in)    :: dtime
! local
  type(scalar_1d_field)      :: scalar_normal_velocity_tmpnew
  type(scalar_1d_field)      :: scalar_normal_velocity_tend1
  type(scalar_1d_field)      :: scalar_normal_velocity_tend2
  type(scalar_1d_field)      :: scalar_normal_velocity_tend3
  type(scalar_1d_field)      :: scalar_height_tmpnew
  type(scalar_1d_field)      :: scalar_height_tend1
  type(scalar_1d_field)      :: scalar_height_tend2
  type(scalar_1d_field)      :: scalar_height_tend3

!================================================
!    Just initialization, will be overwritten
!================================================

  scalar_normal_velocity_tmpnew = scalar_normal_velocity_at_tr_edge
  scalar_normal_velocity_tend1  = scalar_normal_velocity_at_tr_edge
  scalar_normal_velocity_tend2  = scalar_normal_velocity_at_tr_edge
  scalar_normal_velocity_tend3  = scalar_normal_velocity_at_tr_edge

  scalar_height_tmpnew          = scalar_height_at_dual_cell
  scalar_height_tend1           = scalar_height_at_dual_cell
  scalar_height_tend2           = scalar_height_at_dual_cell
  scalar_height_tend3           = scalar_height_at_dual_cell

!================================================
!                 1st call
!================================================

     call calc_spatial_tend_tr(mesh                             ,&
                               scalar_normal_velocity_at_tr_edge,&
                               scalar_height_at_dual_cell       ,&
                               scalar_height_tend1              ,&
                               itimestep                        ,&
                               dtime, use_limiter_rk3, .true.,1)

!================================================
!               Final Update
!================================================

     scalar_height_at_dual_cell_next%f(:)   = scalar_height_at_dual_cell%f(:)+&
                                              dtime*scalar_height_tend1%f(:)

     return
  end subroutine time_integration_fit_tr


!================================================
!                Forward-In-Time
!================================================

   subroutine time_integration_fit(mesh,scalar_normal_velocity_at_edge     ,&
                                        scalar_height_at_prime_cell        ,&
                                        scalar_topo_at_prime_cell          ,&
                                        scalar_normal_velocity_at_edge_next,&
                                        scalar_height_at_prime_cell_next   ,&
                                        tend_normal_velocity_at_edge       ,&
                                        tend_height_at_prime_cell          ,&
                                        itimestep,dtime)
! io
   type(global_domain), intent(inout) :: mesh
   type(scalar_1d_field),   intent(in)    :: scalar_normal_velocity_at_edge
   type(scalar_1d_field),   intent(in)    :: scalar_height_at_prime_cell
   type(scalar_1d_field),   intent(in)    :: scalar_topo_at_prime_cell
   type(scalar_1d_field),   intent(inout) :: scalar_normal_velocity_at_edge_next
   type(scalar_1d_field),   intent(inout) :: scalar_height_at_prime_cell_next
   type(scalar_1d_field),   intent(inout) :: tend_normal_velocity_at_edge
   type(scalar_1d_field),   intent(inout) :: tend_height_at_prime_cell
   integer(i4)       ,   intent(in)    :: itimestep
   real(r8)          ,   intent(in)    :: dtime

    call calc_spatial_tend(mesh,scalar_normal_velocity_at_edge,& ! in
                                scalar_height_at_prime_cell   ,& ! in
                                scalar_topo_at_prime_cell     ,& ! in
                                tend_normal_velocity_at_edge  ,& ! out
                                tend_height_at_prime_cell     ,& ! out
                                itimestep, dtime, use_limiter_rk3,1)       ! in

    scalar_normal_velocity_at_edge_next%f(:) = scalar_normal_velocity_at_edge%f(:)+&
                                               dtime*tend_normal_velocity_at_edge%f(:)

    scalar_height_at_prime_cell_next%f(:)    = scalar_height_at_prime_cell%f(:)+&
                                               dtime*tend_height_at_prime_cell%f(:)

     return
   end subroutine time_integration_fit

  end module swe_time_integration
