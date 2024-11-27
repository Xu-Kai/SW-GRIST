
!----------------------------------------------------------------------------
! Created on 2016
! Version 1.0
! Description: Module for calc spatial derivatives of SWE
! Revision history: 
!----------------------------------------------------------------------------

 module swe_spatial_module

  use grist_constants,        only: gravity, i4, r8, rearth
  use grist_domain_types,     only: global_domain
  use grist_data_types,       only: scalar_1d_field,exchange_field_list_1d
  use grist_flux_operators,   only: calc_normal_flux_hx_edge
  use grist_gcd_module,       only: gradient_operator, divergence_operator
  use grist_limiter_module,   only: flux_limiter_fct_hx!, flux_limiter_thuburn_hx
  use grist_advection_module, only: divergence_operator_tspas
  use grist_hori_swe_module,  only: calc_kinetic_energy   , &
                                    calc_tend_velocity_coriolis
  use grist_nml_module,       only: advection_scheme      , &
                                    pv_order              , &
                                    isolate_advection_test, &
                                    use_tspas, use_streamf, &
                                    limiter_type          , &
                                    stencil_width         , &
                                    stencil_exchange_flag
  use swe_init_module,        only: overwrite_wind, overwrite_wind_streamf
  use swe_vars_module,        only: scalar_n=>scalar_height_at_prime_cell, &
                                    scalar_tangent_velocity_at_edge ! only used for FFSL
  use grist_config_partition, only: debug_data_1d,exchange_data_1d_add,exchange_data_1d

  implicit none

   real(r8) :: global_max_cfl  ! max cfl of all time
  
   private
   public   :: calc_spatial_tend
 

  contains

!================================================
!     Spatial tendency for velocity and height
!================================================
  
   subroutine calc_spatial_tend(mesh,scalar_normal_velocity_at_edge   ,& ! input of V
                                     scalar_height_at_prime_cell      ,& ! input of H
                                     scalar_topo_at_prime_cell        ,& ! input of TOPO
                                     tend_normal_velocity_at_edge     ,& ! tend of V
                                     tend_scalar_height_at_prime_cell ,& ! tend of H
                                     itimestep,timestep,use_limiter,irk_step)     ! time index

   use swe_vars_module,  only: scalar_normal_flux_at_edge
! io
   type(global_domain),  intent(inout) :: mesh
   type(scalar_1d_field),intent(in)    :: scalar_normal_velocity_at_edge
   type(scalar_1d_field),intent(in)    :: scalar_height_at_prime_cell
   type(scalar_1d_field),intent(in)    :: scalar_topo_at_prime_cell
   type(scalar_1d_field),intent(inout) :: tend_normal_velocity_at_edge
   type(scalar_1d_field),intent(inout) :: tend_scalar_height_at_prime_cell
   integer(i4)       ,   intent(in)    :: itimestep
   real(r8)          ,   intent(in)    :: timestep
   logical           ,   intent(in)    :: use_limiter
   integer(i4)       ,   intent(in)    :: irk_step
! local
   type(scalar_1d_field)               :: scalar_normal_velocity_at_edge_local
   type(scalar_1d_field)               :: scalar_normal_flux_at_edge_lo
   type(scalar_1d_field)               :: scalar_normal_flux_at_edge_lm
   real(r8)                            :: cfl,max_cfl
   integer(i4)                         :: ie
   type(exchange_field_list_1d),pointer:: field_head

   field_head => null()
! 

     allocate(scalar_normal_velocity_at_edge_local%f(1:mesh%ne_full))
     allocate(scalar_normal_flux_at_edge_lo%f(1:mesh%ne_full))
     allocate(scalar_normal_flux_at_edge_lm%f(1:mesh%ne_full))

     if(isolate_advection_test)then
        if(use_limiter.or.advection_scheme.eq.99.or.use_tspas) then
          mesh%nv = mesh%nv_halo(2)
          mesh%nt = mesh%nt_halo(2)
          mesh%ne = mesh%ne_halo(2)
        end if
        if(use_streamf)then
          call overwrite_wind_streamf(mesh, scalar_normal_velocity_at_edge_local,&
                                            scalar_tangent_velocity_at_edge, &
                                            itimestep, timestep, .false.)
        else
           call overwrite_wind(mesh, scalar_normal_velocity_at_edge_local,&
                                     scalar_tangent_velocity_at_edge, &
                                     itimestep, timestep, .false.)
        end if
        mesh%nv = mesh%nv_compute
        mesh%nt = mesh%nt_compute
        mesh%ne = mesh%ne_compute
     else
        scalar_normal_velocity_at_edge_local%f = scalar_normal_velocity_at_edge%f
     end if

   !  max_cfl      = 0._r8
   !  do ie = 1, mesh%ne
   !     cfl       = scalar_normal_velocity_at_edge_local%f(ie)*timestep/(mesh%edt(ie)%leng*rearth)
   !     max_cfl   = max(cfl,max_cfl)
   !  end do
   !  global_max_cfl = max(max_cfl,global_max_cfl)
   !  print*,"max_cfl=",max_cfl, "global_maxcfl=",global_max_cfl
    
     IF(use_tspas)then ! tspas can only be used with FIT and isolated transport test
         mesh%ne = mesh%ne_halo(1)
         mesh%nv = mesh%nv_halo(1)
         call divergence_operator_tspas(mesh,scalar_normal_velocity_at_edge_local%f,&
                                             scalar_height_at_prime_cell%f,         &
                                             timestep,&
                                             tend_scalar_height_at_prime_cell%f) 
         mesh%ne = mesh%ne_compute
         mesh%nv = mesh%nv_compute
     ELSE

!================================================
!        calculate normal mass flux
!================================================
        if(.not. stencil_exchange_flag) mesh%ne = mesh%ne_halo(1)
        call calc_normal_flux_hx_edge(mesh,scalar_normal_velocity_at_edge_local%f,& ! in
                                           scalar_tangent_velocity_at_edge%f     ,&
                                           scalar_height_at_prime_cell%f         ,& ! in
                                           scalar_normal_flux_at_edge%f          ,& ! out
                                           advection_scheme, &
                                           timestep)
        scalar_normal_flux_at_edge%f = scalar_normal_flux_at_edge%f*scalar_normal_velocity_at_edge_local%f
        if(.not. stencil_exchange_flag) mesh%ne = mesh%ne_compute

        if(isolate_advection_test.and.use_limiter)then
          select case(trim(limiter_type))
          case('fct')

            mesh%ne = mesh%ne_halo(2)
            call calc_normal_flux_hx_edge(mesh,scalar_normal_velocity_at_edge_local%f,& ! in
                                               scalar_tangent_velocity_at_edge%f     ,& ! dum
                                               scalar_n%f                            ,& ! in
                                               scalar_normal_flux_at_edge_lo%f       ,& ! out
                                               1, &
                                               timestep  )
            scalar_normal_flux_at_edge_lo%f = scalar_normal_flux_at_edge_lo%f*scalar_normal_velocity_at_edge_local%f
            mesh%ne = mesh%ne_compute
            !exchange scheme
            if(stencil_exchange_flag) then
              call exchange_data_1d_add(mesh,field_head,scalar_normal_flux_at_edge)
              call exchange_data_1d(mesh%local_block,field_head)
            end if
  
            !out compute,in halo(2);halo(1);halo(2)
            call flux_limiter_fct_hx(mesh,scalar_n%f                     ,&
                                          scalar_normal_flux_at_edge%f   ,&
                                          scalar_normal_flux_at_edge_lo%f,&
                                          scalar_normal_flux_at_edge_lm%f,&
                                          timestep)

            scalar_normal_flux_at_edge%f = scalar_normal_flux_at_edge_lm%f
          !case('thuburn') ! input is q-edge
              ! note that this line may be ill-conditioned since velocity may be zero, but just for a loose test
              ! on mac this prblem comes out when using streamf as init
          !    scalar_normal_flux_at_edge_lo%f = scalar_normal_flux_at_edge%f/scalar_normal_velocity_at_edge_local%f
          !    call flux_limiter_thuburn_hx(mesh,scalar_normal_velocity_at_edge_local%f, & ! in
          !                                      scalar_n%f                            , &
          !                                      scalar_normal_flux_at_edge_lo%f       , &
          !                                      scalar_normal_flux_at_edge_lm%f       , &
          !                                      timestep)
          !    scalar_normal_flux_at_edge%f = scalar_normal_flux_at_edge_lm%f*scalar_normal_velocity_at_edge_local%f
          case default
               print*,"you must select a limiter between fct and thuburn"
          end select
        end if

        if(.not.isolate_advection_test)then
          !exchange scheme
           if(stencil_exchange_flag) then
             call exchange_data_1d_add(mesh,field_head,scalar_normal_flux_at_edge)
             call exchange_data_1d(mesh%local_block,field_head)
           end if
   
           call calc_tend_velocity(mesh,scalar_normal_velocity_at_edge_local,&
                                        scalar_height_at_prime_cell   ,&
                                        scalar_topo_at_prime_cell     ,&
                                        scalar_normal_flux_at_edge    ,&
                                        tend_normal_velocity_at_edge  ,&
                                        itimestep,timestep,irk_step)
        end if
!
! advection tendency
!
        call divergence_operator(mesh,scalar_normal_flux_at_edge%f,&
                                      tend_scalar_height_at_prime_cell%f)

        tend_scalar_height_at_prime_cell%f(:) = -1._r8*tend_scalar_height_at_prime_cell%f(:)

     END IF

     deallocate(scalar_normal_velocity_at_edge_local%f)
     deallocate(scalar_normal_flux_at_edge_lo%f)
     deallocate(scalar_normal_flux_at_edge_lm%f)

     return
  end subroutine calc_spatial_tend

!*************************
!         PRIVATE
!*************************

  subroutine calc_tend_velocity(mesh,scalar_normal_velocity_at_edge,&
                                     scalar_height_at_prime_cell   ,&
                                     scalar_topo_at_prime_cell     ,&
                                     scalar_normal_flux_at_edge    ,&
                                     tend_normal_velocity_at_edge  ,&
                                     itimestep, timestep, irk_step)
! io
   type(global_domain), intent(inout) :: mesh
   type(scalar_1d_field),   intent(in)    :: scalar_normal_velocity_at_edge
   type(scalar_1d_field),   intent(in)    :: scalar_height_at_prime_cell
   type(scalar_1d_field),   intent(in)    :: scalar_topo_at_prime_cell
   type(scalar_1d_field),   intent(in)    :: scalar_normal_flux_at_edge
   type(scalar_1d_field),   intent(inout) :: tend_normal_velocity_at_edge
   integer(i4)          ,   intent(in)    :: itimestep
   real(r8)             ,   intent(in)    :: timestep
   integer(i4)          ,   intent(in)    :: irk_step
! local
   type(scalar_1d_field)                  :: tend_velocity_energy
   type(scalar_1d_field)                  :: tend_velocity_coriolis
   type(scalar_1d_field)                  :: kinetic_energy
   type(scalar_1d_field)                  :: energy_for_grad
   integer(i4)                            :: ie , iv


   if(.not.allocated(tend_velocity_energy%f))   allocate(tend_velocity_energy%f(1:mesh%ne_full))
   if(.not.allocated(tend_velocity_coriolis%f)) allocate(tend_velocity_coriolis%f(1:mesh%ne_full))
   if(.not.allocated(kinetic_energy%f))         allocate(kinetic_energy%f(1:mesh%nv_full))
   if(.not.allocated(energy_for_grad%f))        allocate(energy_for_grad%f(1:mesh%nv_full))
 
   mesh%nv = mesh%nv_halo(1)
   mesh%nt = mesh%nt_halo(1)
    call calc_kinetic_energy(mesh,scalar_normal_velocity_at_edge%f,kinetic_energy%f)

!
! KE+PE
!
    do iv = 1, mesh%nv
       energy_for_grad%f(iv) = kinetic_energy%f(iv)+gravity*(scalar_height_at_prime_cell%f(iv)+&
                                                             scalar_topo_at_prime_cell%f(iv))
    end do

    mesh%nv = mesh%nv_compute
    mesh%nt = mesh%nt_compute

    call gradient_operator(mesh,energy_for_grad%f,tend_velocity_energy%f)
!
! GRAD(KE+PE)
!
    tend_velocity_energy%f(:) = -1._r8*tend_velocity_energy%f(:)


    call calc_tend_velocity_coriolis(mesh,scalar_height_at_prime_cell%f    ,&
                                          scalar_normal_velocity_at_edge%f ,&
                                          scalar_normal_flux_at_edge%f     ,&
                                          tend_velocity_coriolis%f,timestep,irk_step)

    tend_normal_velocity_at_edge%f(:) = tend_velocity_energy%f(:)+tend_velocity_coriolis%f(:)
!
! clean 
!
   deallocate(tend_velocity_energy%f)
   deallocate(tend_velocity_coriolis%f)
   deallocate(kinetic_energy%f)
   deallocate(energy_for_grad%f)
 
   return
  end subroutine calc_tend_velocity

  end module swe_spatial_module
