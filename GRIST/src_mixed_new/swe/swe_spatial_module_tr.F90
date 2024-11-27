
!----------------------------------------------------------------------------
! Created on 2016
! Version 1.0
! Description: Divergence operator at dual cell
! Revision history: 
!----------------------------------------------------------------------------

 module swe_spatial_module_tr
!
! grist
!
  use grist_constants,      only: i4, r8
  use grist_domain_types,   only: global_domain
  use grist_data_types,     only: scalar_1d_field
  use grist_flux_operators, only: calc_normal_flux_tr_edge
  use grist_gcd_module    , only: divergence_operator_tr,divergence_operator_tr_smth
  use grist_limiter_module, only: flux_limiter_fct_tr 
  use grist_nml_module,     only: pv_order
!
! swe
!
  use swe_init_module,      only: overwrite_wind 
  use swe_vars_module,      only: scalar_n=>scalar_height_at_dual_cell, &
                                  scalar_tangent_velocity_at_edge
 
  implicit none

  private
  public :: calc_spatial_tend_tr

  contains

  subroutine calc_spatial_tend_tr(mesh,scalar_normal_velocity_at_tr_edge ,& ! input of V
                                       scalar_height_at_dual_cell        ,& ! input of H
                                       tend_scalar_height_at_dual_cell   ,&
                                       itimestep                         ,&
                                       timestep                          ,&
                                       use_limiter                       ,&
                                       called_by_tr, irk_step    ) ! called by tr for isolated advection test
! io
   type(global_domain),  intent(inout) :: mesh
   type(scalar_1d_field),intent(in)    :: scalar_normal_velocity_at_tr_edge
   type(scalar_1d_field),intent(in)    :: scalar_height_at_dual_cell
   type(scalar_1d_field),intent(inout) :: tend_scalar_height_at_dual_cell
   integer(i4)       ,   intent(in)    :: itimestep
   real(r8)          ,   intent(in)    :: timestep
   logical           ,   intent(in)    :: use_limiter
   logical           ,   intent(in)    :: called_by_tr
   integer(i4)       ,   intent(in)    :: irk_step
! local
   integer(i4)             :: it
   integer(i4)             :: ie
   type(scalar_1d_field)   :: scalar_normal_velocity_at_tr_edge_local
   type(scalar_1d_field)   :: scalar_at_dual_edge
   type(scalar_1d_field)   :: flux_edge
   type(scalar_1d_field)   :: flux_up
   type(scalar_1d_field)   :: flux_lm


    allocate(scalar_normal_velocity_at_tr_edge_local%f(1:mesh%ne))
    allocate(scalar_at_dual_edge%f(1:mesh%ne))
    allocate(flux_edge%f(1:mesh%ne))
    allocate(flux_up%f(1:mesh%ne))
    allocate(flux_lm%f(1:mesh%ne))

    scalar_normal_velocity_at_tr_edge_local%pos = 2

    call overwrite_wind(mesh, scalar_normal_velocity_at_tr_edge_local,&
                              scalar_tangent_velocity_at_edge, &
                              itimestep, timestep,.true.)

!================================================
!       current order of flux operator
!================================================

    call calc_normal_flux_tr_edge(mesh,scalar_height_at_dual_cell%f         ,&
                                       scalar_normal_velocity_at_tr_edge%f  ,&
                                       scalar_tangent_velocity_at_edge%f    ,&  ! for ffsl
                                       scalar_at_dual_edge%f                ,&
                                       called_by_tr                         ,&
                                       pv_order(irk_step), timestep)
    do ie = 1, mesh%ne
       flux_edge%f(ie) = scalar_normal_velocity_at_tr_edge%f(ie)*scalar_at_dual_edge%f(ie)
    end do

!================================================
!                   limiter 
!================================================

    if(use_limiter) then
       call calc_normal_flux_tr_edge(mesh,scalar_n%f                          , & ! in
                                          scalar_normal_velocity_at_tr_edge%f , & ! in
                                          scalar_tangent_velocity_at_edge%f   , & ! for ffsl
                                          scalar_at_dual_edge%f               , & ! out
                                          called_by_tr                        , & ! in
                                          1   , timestep)
       do ie = 1, mesh%ne
          flux_up%f(ie) = scalar_normal_velocity_at_tr_edge%f(ie)*scalar_at_dual_edge%f(ie)
       end do

       call flux_limiter_fct_tr(mesh,scalar_n%f  ,&
                                     flux_edge%f ,&
                                     flux_up%f   ,&
                                     flux_lm%f   ,&
                                     timestep)
       flux_edge = flux_lm
    end if

    call divergence_operator_tr(mesh,flux_edge%f, tend_scalar_height_at_dual_cell%f)
    !call divergence_operator_tr_smth(mesh,flux_edge%f, tend_scalar_height_at_dual_cell%f)

    do it = 1, mesh%nt
       tend_scalar_height_at_dual_cell%f(it) = -1._r8*tend_scalar_height_at_dual_cell%f(it)
    end do

    deallocate(scalar_normal_velocity_at_tr_edge_local%f)
    deallocate(scalar_at_dual_edge%f)
    deallocate(flux_edge%f)
    deallocate(flux_up%f)
    deallocate(flux_lm%f)

    return
  end subroutine calc_spatial_tend_tr

  end module swe_spatial_module_tr
