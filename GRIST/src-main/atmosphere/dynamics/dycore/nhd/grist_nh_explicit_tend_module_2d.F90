
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description: Explicit Tendency of w-face and phi-face prepared before nh 
!              integration. these tends inlcude horizontal and vertical
!              tends, and are computed as residules of flux
! Revision history:
!              1. testing hori/vert method2 based on direct evaluation
!----------------------------------------------------------------------------

  module grist_nh_explicit_tend_module_2d
! para
    use grist_constants,              only: i4, r8, fillvalue, half, rearth, gravity, zero
#ifndef SEQ_GRIST
    use grist_data_types,             only: scalar_1d_field, scalar_2d_field, exchange_field_list_2d, exchange_field_list_1d
#else
    use grist_data_types,             only: scalar_1d_field, scalar_2d_field
#endif
    use grist_domain_types,           only: global_domain
    use grist_nml_module,             only: ver_adv_flag, nlev, nlevp, stencil_exchange_flag, write_verbose, eqs_vert_diff
! hori
    use grist_dycore_primal_flux_operators_2d,only: calc_primal_normal_flux_at_edge, calc_primal_normal_flux_at_edge_var2
    use grist_dycore_gcd_recon_module_2d,     only: divergence_operator_2d, divergence_operator_2d_var2
    use grist_dycore_diffusion_module,        only: perot_weight_at_pc
    use grist_dycore_ref_atmos,               only: scalar_geop_at_pc_face_level_bar
! vert
    use grist_hpe_constants,          only: deta_full, deta_face, deta_face_b
    use grist_hpe_vertical_advection, only: calc_hpe_tend_vert_mass_flux_face
! mpi
#ifndef SEQ_GRIST
    use grist_config_partition,       only: exchange_data_2d_add , exchange_data_2d, exchange_data_1d_add , exchange_data_1d, debug_data_2d
#endif
    use grist_mpi

    implicit none

     private
     public ::  grist_nh_et_init , &
                grist_nh_et_final, &
                grist_nh_et_adv_face!, &
                !grist_nh_et_www_face

      real(r8), allocatable   :: scalar_template_1d_nlev_a(:)
      real(r8), allocatable   :: scalar_template_1d_nlevp_a(:)
      real(r8), allocatable   :: scalar_template_1d_nlevp_b(:)
      real(r8), allocatable   :: scalar_template_1d_nlevp_c(:)

   CONTAINS

!===================================================
! same as grist_nh_et_www_face but evaluate 2 vars
!===================================================

   subroutine grist_nh_et_adv_face(mesh, dtime,  &
                                   tend_hpressure_at_pc_surface_cnty          ,& ! part(hps)/part(t)
                                   scalar_normal_mass_flux_at_edge_face_level ,& ! DELHP*V
                                   tend_mass_at_pc_face_level_hori_rk         ,& ! div.(mass)
                                   scalar_eta_mass_flux_at_pc_full_level      ,& ! m*etadot
                                   scalar_eta_mass_flux_at_pc_face_level      ,& ! m*etadot
                                   scalar_delhp_at_pc_face_level_rk           ,&
                                   scalar_delhp_at_pc_face_level_n            ,&
                                   scalar_delhp_at_pc_face_level_np1          ,&
                                   scalar_www_at_pc_face_level                ,&
                                   scalar_phi_at_pc_face_level                ,&
                                   tend_et_www_face_level                     ,& 
                                   tend_et_phi_face_level                     ,&
                                   adv_face_flag)
! io
      use omp_lib
      type(global_domain)  , intent(inout) :: mesh
      real(r8)             , intent(in)    :: dtime
      type(scalar_1d_field), intent(in)    :: tend_hpressure_at_pc_surface_cnty
      type(scalar_2d_field), intent(in)    :: scalar_normal_mass_flux_at_edge_face_level
      type(scalar_2d_field), intent(in)    :: tend_mass_at_pc_face_level_hori_rk
      type(scalar_2d_field), intent(in)    :: scalar_eta_mass_flux_at_pc_full_level
      type(scalar_2d_field), intent(in)    :: scalar_eta_mass_flux_at_pc_face_level
      type(scalar_2d_field), intent(in)    :: scalar_delhp_at_pc_face_level_rk
      type(scalar_2d_field), intent(in)    :: scalar_delhp_at_pc_face_level_n
      type(scalar_2d_field), intent(in)    :: scalar_delhp_at_pc_face_level_np1
      type(scalar_2d_field), intent(in)    :: scalar_www_at_pc_face_level
      type(scalar_2d_field), intent(in)    :: scalar_phi_at_pc_face_level
      type(scalar_2d_field), intent(inout) :: tend_et_www_face_level
      type(scalar_2d_field), intent(inout) :: tend_et_phi_face_level
      integer(i4)          , intent(in)    :: adv_face_flag
! local
      type(scalar_2d_field)                :: scalar_phipie_at_pc_face_level
      type(scalar_2d_field)                :: tend1_hori_final, tend2_hori_final
      type(scalar_2d_field)                :: tend1_vert_tmp,   tend2_vert_tmp
      type(scalar_2d_field)                :: tend1_vert_final, tend2_vert_final
      type(scalar_2d_field)                :: scalar_template_2d_ne_a
      type(scalar_2d_field)                :: scalar_template_2d_ne_b
      type(scalar_1d_field)                :: scalar_template_1d_ne_a
      type(scalar_1d_field)                :: scalar_template_1d_ne_b
#ifndef SEQ_GRIST
      type(exchange_field_list_1d),pointer :: field_head_1d
      type(exchange_field_list_2d),pointer :: field_head_2d
#endif
      integer(i4)                          :: iv, ilev, icell1, icell2, ie, inb
      real(r8)                             :: vector_sum1(3), vector_sum2(3), vector_sum3(3)

#ifndef SEQ_GRIST 
        field_head_1d => null()
        field_head_2d => null()
#endif
! init
        if(.not.allocated(scalar_template_2d_ne_a%f)) allocate(scalar_template_2d_ne_a%f(nlevp,mesh%ne_full))
        if(.not.allocated(scalar_template_2d_ne_b%f)) allocate(scalar_template_2d_ne_b%f(nlevp,mesh%ne_full))
        if(.not.allocated(scalar_template_1d_ne_a%f)) allocate(scalar_template_1d_ne_a%f(mesh%ne_full))
        if(.not.allocated(scalar_template_1d_ne_b%f)) allocate(scalar_template_1d_ne_b%f(mesh%ne_full))

        scalar_template_2d_ne_a%pos  = 6; scalar_template_2d_ne_a%f    = zero
        scalar_template_2d_ne_b%pos  = 6; scalar_template_2d_ne_b%f    = zero
 
        tend1_hori_final              = scalar_www_at_pc_face_level   ! just init
        tend1_vert_final              = scalar_www_at_pc_face_level   ! just init
        tend1_vert_tmp                = scalar_www_at_pc_face_level   ! just init
        tend2_hori_final              = scalar_phi_at_pc_face_level   ! just init
        tend2_vert_final              = scalar_phi_at_pc_face_level   ! just init
        tend2_vert_tmp                = scalar_phi_at_pc_face_level   ! just init
        scalar_phipie_at_pc_face_level= scalar_phi_at_pc_face_level

! info
        if(mpi_rank() == 0.and.write_verbose) print*,"before compute horizontal face level mass flux diff"

!
! hori evaluation
!
        call calc_primal_normal_flux_at_edge_var2(mesh, scalar_normal_mass_flux_at_edge_face_level%f, &
                                                        scalar_normal_mass_flux_at_edge_face_level%f, &
                                                        scalar_www_at_pc_face_level%f               , &
                                                        scalar_phipie_at_pc_face_level%f            , &
                                                        scalar_template_2d_ne_a%f                   , &
                                                        scalar_template_2d_ne_b%f                   , &
                                                        adv_face_flag, dtime, nlev+1, .true.)

        scalar_template_2d_ne_a%f = scalar_template_2d_ne_a%f*scalar_normal_mass_flux_at_edge_face_level%f
        scalar_template_2d_ne_b%f = scalar_template_2d_ne_b%f*scalar_normal_mass_flux_at_edge_face_level%f

#ifndef SEQ_GRIST
        call exchange_data_2d_add(mesh,field_head_2d,scalar_template_2d_ne_a)
        call exchange_data_2d_add(mesh,field_head_2d,scalar_template_2d_ne_b)
        call exchange_data_2d(mesh%local_block,field_head_2d)
#endif

        mesh%ne = mesh%ne_halo(1)
        mesh%nv = mesh%nv_halo(1)
        call divergence_operator_2d_var2(mesh, scalar_template_2d_ne_a%f, &
                                               scalar_template_2d_ne_b%f, &
                                               tend1_hori_final%f       , &
                                               tend2_hori_final%f       , &
                                               nlev+1)
        tend1_hori_final%f = tend1_hori_final%f-scalar_www_at_pc_face_level%f   *tend_mass_at_pc_face_level_hori_rk%f
        tend2_hori_final%f = tend2_hori_final%f-scalar_phipie_at_pc_face_level%f*tend_mass_at_pc_face_level_hori_rk%f

        mesh%ne = mesh%ne_compute
        mesh%nv = mesh%nv_compute

! info
        if(mpi_rank() == 0.and.write_verbose) print*,"finish compute horizontal face level mass flux diff"
!
! vert evaluation
!
           call  calc_hpe_tend_vert_mass_flux_face_2d_var2(mesh,&
                                                           scalar_www_at_pc_face_level%f    , & ! face level scalar1
                                                           scalar_phipie_at_pc_face_level%f , & ! face level scalar2
                                                           scalar_eta_mass_flux_at_pc_full_level%f , & ! full level etadot
                                                           ver_adv_flag                  , &
                                                           tend1_vert_tmp%f, tend2_vert_tmp%f)    ! face level tend

           scalar_template_1d_nlevp_c(1)      = zero
           scalar_template_1d_nlevp_c(nlev+1) = zero
!$omp parallel  private(iv) firstprivate(scalar_template_1d_nlevp_c)
!$omp do schedule(static,5)
           do iv = 1, mesh%nv_halo(1)
              scalar_template_1d_nlevp_c(2:nlev) = scalar_eta_mass_flux_at_pc_full_level%f(2:nlev,iv)-scalar_eta_mass_flux_at_pc_full_level%f(1:nlev-1,iv)
              tend1_vert_final%f(:,iv)    = tend1_vert_tmp%f(:,iv)-scalar_www_at_pc_face_level%f(:,iv)*scalar_template_1d_nlevp_c(:)
              tend2_vert_final%f(:,iv)    = tend2_vert_tmp%f(:,iv)-scalar_phipie_at_pc_face_level%f(:,iv)*scalar_template_1d_nlevp_c(:)
           end do
!$omp end do nowait
!$omp end parallel

! ASSIGN; still divided by local delhp at face leve
        tend_et_www_face_level%f   = (tend1_hori_final%f+tend1_vert_final%f)/scalar_delhp_at_pc_face_level_rk%f
        tend_et_phi_face_level%f   = (tend2_hori_final%f+tend2_vert_final%f)/scalar_delhp_at_pc_face_level_rk%f
! info
        if(mpi_rank() == 0.and.write_verbose) print*,"finish compute vertical www&phi mass flux diff"

! clean
        deallocate(scalar_template_2d_ne_a%f,scalar_template_2d_ne_b%f,scalar_template_1d_ne_a%f,scalar_template_1d_ne_b%f)
        deallocate(tend1_hori_final%f,tend2_hori_final%f,tend1_vert_final%f,tend2_vert_final%f)

     return
   end subroutine grist_nh_et_adv_face

   subroutine grist_nh_et_init(mesh)
! io
     type(global_domain),  intent(in) :: mesh
! local
      if(.not.allocated(scalar_template_1d_nlev_a))  allocate(scalar_template_1d_nlev_a(nlev))
      if(.not.allocated(scalar_template_1d_nlevp_a)) allocate(scalar_template_1d_nlevp_a(nlev+1))
      if(.not.allocated(scalar_template_1d_nlevp_b)) allocate(scalar_template_1d_nlevp_b(nlev+1))
      if(.not.allocated(scalar_template_1d_nlevp_c)) allocate(scalar_template_1d_nlevp_c(nlev+1))

     return
   end subroutine grist_nh_et_init
   
   subroutine grist_nh_et_final

      deallocate(scalar_template_1d_nlev_a)
      deallocate(scalar_template_1d_nlevp_a)
      deallocate(scalar_template_1d_nlevp_b)
      deallocate(scalar_template_1d_nlevp_c)

   end subroutine grist_nh_et_final

!================================================================
!                 PRIVATE ROUTINES BELOW
!================================================================

   subroutine calc_hpe_tend_vert_mass_flux_face_2d(mesh,&
                                                   scalar_at_pc_face_level             , &
                                                   scalar_eta_mass_velocity_at_pc_full , &
                                                   order                               , &
                                                   tend_vert_mass_flux_at_pc_face)
! io
    use omp_lib
    type(global_domain),    intent(in)   :: mesh
    real(r8), allocatable,  intent(in)   :: scalar_at_pc_face_level(:,:)
    real(r8), allocatable,  intent(in)   :: scalar_eta_mass_velocity_at_pc_full(:,:)
    integer(i4)          ,  intent(in)   :: order
    real(r8), allocatable,  intent(inout):: tend_vert_mass_flux_at_pc_face(:,:)  ! with minus sign
! local
    real(r8),  allocatable               :: scalar_at_pc_full_level(:,:)
    integer(i4)                          :: iv,ilev
    integer(i4)                          :: ii

     if(.not.allocated(scalar_at_pc_full_level)) allocate(scalar_at_pc_full_level(nlev,mesh%nv_full))

     call calc_hpe_vert_flux_operator_face_2d(mesh,&
                                              scalar_at_pc_face_level     ,&
                                              scalar_eta_mass_velocity_at_pc_full ,&
                                              order                    ,&
                                              scalar_at_pc_full_level )

!$omp parallel  private(ii,iv,ilev) 
!$omp do schedule(dynamic,100) 
    do ii=1,mesh%nv_halo(1)*(nlev-1),1
        iv=ceiling(ii/real(nlev-1,r8))
        ilev=1+ii-(iv-1)*(nlev-1)
!     do iv = 1, mesh%nv_halo(1)
!     do ilev =  2, nlev
! no minus sign here
        tend_vert_mass_flux_at_pc_face(ilev,iv) = (scalar_eta_mass_velocity_at_pc_full(ilev,iv)*scalar_at_pc_full_level(ilev,iv)-&
                                                  scalar_eta_mass_velocity_at_pc_full(ilev-1,iv)*scalar_at_pc_full_level(ilev-1,iv))

!     end do
!     end do
    end do
!$omp end do nowait
!$omp end parallel

     tend_vert_mass_flux_at_pc_face(1,:)      = 0._r8
     tend_vert_mass_flux_at_pc_face(nlev+1,:) = 0._r8

     deallocate(scalar_at_pc_full_level)

     return
   end subroutine calc_hpe_tend_vert_mass_flux_face_2d

!
! private, from face level to full level
!

   subroutine calc_hpe_vert_flux_operator_face_2d(mesh,&
                                                  scalar_at_pc_face_level            ,&
                                                  scalar_eta_mass_velocity_at_pc_full,&
                                                  order                              ,&
                                                  scalar_at_pc_full_level)
!
    use omp_lib
    type(global_domain),     intent(in)     :: mesh
    real(r8), allocatable,   intent(in)     :: scalar_at_pc_face_level(:,:)
    real(r8), allocatable,   intent(in)     :: scalar_eta_mass_velocity_at_pc_full(:,:)
    integer               ,  intent(in)     :: order
    real(r8), allocatable,   intent(inout)  :: scalar_at_pc_full_level(:,:)
! local
    integer(i4)    :: ilev, iv
    real(r8)       :: part1, der_k, der_kp1
    integer(i4)                             :: ii
    

!$omp parallel  private(iv) 
!$omp do schedule(static,5) 
     do iv = 1, mesh%nv_halo(1)
        scalar_at_pc_full_level(1,iv)      = half*(scalar_at_pc_face_level(1,iv)   +scalar_at_pc_face_level(2,iv))
        scalar_at_pc_full_level(nlev,iv)   = half*(scalar_at_pc_face_level(nlev,iv)+scalar_at_pc_face_level(nlev+1,iv))
     end do
!$omp end do nowait
!$omp end parallel 

     IF(EQS_VERT_DIFF)THEN
     select case(order)
     case(2)
       do iv = 1, mesh%nv_halo(1)
         do ilev = 2, nlev-1
            scalar_at_pc_full_level(ilev,iv) = half*(scalar_at_pc_face_level(ilev,iv)+scalar_at_pc_face_level(ilev+1,iv))
         end do
       end do
     case(3) 
!$omp parallel  private(ii,iv,ilev)      
!$omp do schedule(dynamic,100) 
    do ii=1,mesh%nv_halo(1)*(nlev-2),1
        iv=ceiling(ii/real(nlev-2,r8))
        ilev=1+ii-(iv-1)*(nlev-2)
!       do iv = 1, mesh%nv_halo(1)
!         do ilev = 2, nlev-1
            scalar_at_pc_full_level(ilev,iv) = (scalar_at_pc_face_level(ilev,iv)+scalar_at_pc_face_level(ilev+1,iv))*(7._r8/12._r8)-&
                                           (scalar_at_pc_face_level(ilev+2,iv)+scalar_at_pc_face_level(ilev-1,iv))*(1._r8/12._r8)+&
                                 sign(1._r8,scalar_eta_mass_velocity_at_pc_full(ilev,iv))*&
                                          ((scalar_at_pc_face_level(ilev+2,iv)-scalar_at_pc_face_level(ilev-1,iv))-&
                                     3._r8*(scalar_at_pc_face_level(ilev+1,iv)-scalar_at_pc_face_level(ilev,iv)))/12._r8
!         end do
       !end do
    end do
!$omp end do nowait
!$omp end parallel 
     case default
        print*," you must select a vert order in calc_hpe_vert_flux_operator_face"
        stop
     end select

     ELSE

     select case(order)
     case(2) 
       do iv = 1, mesh%nv_halo(1)
         do ilev = 2, nlev-1
            scalar_at_pc_full_level(ilev,iv) = half*(scalar_at_pc_face_level(ilev,iv)+scalar_at_pc_face_level(ilev+1,iv))
         end do
       end do
     case(3)
      do iv = 1, mesh%nv_halo(1)
         do ilev = 2, nlev-1
            part1 = half*(scalar_at_pc_face_level(ilev,iv)+scalar_at_pc_face_level(ilev+1,iv))

            der_k = 2._r8*(deta_face(ilev)**2)*((scalar_at_pc_face_level(ilev+1,iv)-scalar_at_pc_face_level(ilev,iv))/deta_full(ilev)-&
                                                (scalar_at_pc_face_level(ilev,iv)-scalar_at_pc_face_level(ilev-1,iv))/deta_full(ilev-1))/&
                                                (deta_full(ilev)+deta_full(ilev-1))

            der_kp1 = 2._r8*(deta_face(ilev+1)**2)*((scalar_at_pc_face_level(ilev+2,iv)-scalar_at_pc_face_level(ilev+1,iv))/deta_full(ilev+1)-&
                                                    (scalar_at_pc_face_level(ilev+1,iv)-scalar_at_pc_face_level(ilev,iv))/deta_full(ilev))/&
                                                    (deta_full(ilev+1)+deta_full(ilev))

            scalar_at_pc_full_level(ilev,iv)= part1-(der_k+der_kp1)/12._r8+sign(1._r8,scalar_eta_mass_velocity_at_pc_full(ilev,iv))*(der_kp1-der_k)
         end do
      end do
     case default
        print*," you must select a vert order in calc_hpe_vert_flux_operator_face_2d"
        stop
     end select

     END IF
     
     return
   end subroutine calc_hpe_vert_flux_operator_face_2d

! same as above, but for two vars

   subroutine calc_hpe_tend_vert_mass_flux_face_2d_var2(mesh,&
                                                   scalar1_at_pc_face_level            , &
                                                   scalar2_at_pc_face_level            , &
                                                   scalar_eta_mass_velocity_at_pc_full , &
                                                   order                               , &
                                                   tend1_vert_mass_flux_at_pc_face     , &
                                                   tend2_vert_mass_flux_at_pc_face     )
! io
    use omp_lib
    type(global_domain),    intent(in)   :: mesh
    real(r8), allocatable,  intent(in)   :: scalar1_at_pc_face_level(:,:)
    real(r8), allocatable,  intent(in)   :: scalar2_at_pc_face_level(:,:)
    real(r8), allocatable,  intent(in)   :: scalar_eta_mass_velocity_at_pc_full(:,:)
    integer(i4)          ,  intent(in)   :: order
    real(r8), allocatable,  intent(inout):: tend1_vert_mass_flux_at_pc_face(:,:)  ! with minus sign
    real(r8), allocatable,  intent(inout):: tend2_vert_mass_flux_at_pc_face(:,:)  ! with minus sign
! local
    real(r8),  allocatable               :: scalar1_at_pc_full_level(:,:)
    real(r8),  allocatable               :: scalar2_at_pc_full_level(:,:)
    integer(i4)                          :: iv,ilev, ii

    if(.not.allocated(scalar1_at_pc_full_level)) allocate(scalar1_at_pc_full_level(nlev,mesh%nv_full))
    if(.not.allocated(scalar2_at_pc_full_level)) allocate(scalar2_at_pc_full_level(nlev,mesh%nv_full))

    call calc_hpe_vert_flux_operator_face_2d_var2(mesh,&
                                             scalar1_at_pc_face_level   ,&
                                             scalar2_at_pc_face_level   ,&
                                             scalar_eta_mass_velocity_at_pc_full ,&
                                             order                      ,&
                                             scalar1_at_pc_full_level   ,&
                                             scalar2_at_pc_full_level )

!$omp parallel  private(ii,iv,ilev) 
!$omp do schedule(dynamic,100) 
    do ii=1,mesh%nv_halo(1)*(nlev-1),1
        iv=ceiling(ii/real(nlev-1,r8))
        ilev=1+ii-(iv-1)*(nlev-1)
!     do iv = 1, mesh%nv_halo(1)
!     do ilev =  2, nlev
! no minus sign here
        tend1_vert_mass_flux_at_pc_face(ilev,iv) = (scalar_eta_mass_velocity_at_pc_full(ilev  ,iv)*scalar1_at_pc_full_level(ilev  ,iv)-&
                                                    scalar_eta_mass_velocity_at_pc_full(ilev-1,iv)*scalar1_at_pc_full_level(ilev-1,iv))
        tend2_vert_mass_flux_at_pc_face(ilev,iv) = (scalar_eta_mass_velocity_at_pc_full(ilev  ,iv)*scalar2_at_pc_full_level(ilev  ,iv)-&
                                                    scalar_eta_mass_velocity_at_pc_full(ilev-1,iv)*scalar2_at_pc_full_level(ilev-1,iv))

!     end do
!     end do
    end do
!$omp end do nowait
!$omp end parallel

     tend1_vert_mass_flux_at_pc_face(1,:)      = 0._r8; tend2_vert_mass_flux_at_pc_face(1,:)      = 0._r8
     tend1_vert_mass_flux_at_pc_face(nlev+1,:) = 0._r8; tend2_vert_mass_flux_at_pc_face(nlev+1,:) = 0._r8

     deallocate(scalar1_at_pc_full_level)
     deallocate(scalar2_at_pc_full_level)

     return
   end subroutine calc_hpe_tend_vert_mass_flux_face_2d_var2

!
! private, from face level to full level
!

   subroutine calc_hpe_vert_flux_operator_face_2d_var2(mesh,&
                                                  scalar1_at_pc_face_level    ,&
                                                  scalar2_at_pc_face_level    ,&
                                                  scalar_eta_mass_velocity_at_pc_full,&
                                                  order                       ,&
                                                  scalar1_at_pc_full_level    ,&
                                                  scalar2_at_pc_full_level)
!
    use omp_lib
    type(global_domain),     intent(in)     :: mesh
    real(r8), allocatable,   intent(in)     :: scalar1_at_pc_face_level(:,:)
    real(r8), allocatable,   intent(in)     :: scalar2_at_pc_face_level(:,:)
    real(r8), allocatable,   intent(in)     :: scalar_eta_mass_velocity_at_pc_full(:,:)
    integer               ,  intent(in)     :: order
    real(r8), allocatable,   intent(inout)  :: scalar1_at_pc_full_level(:,:)
    real(r8), allocatable,   intent(inout)  :: scalar2_at_pc_full_level(:,:)
! local
    integer(i4)    :: ilev, iv, ii
    real(r8)       :: part1, der_k, der_kp1
    

!$omp parallel  private(iv) 
!$omp do schedule(static,5) 
     do iv = 1, mesh%nv_halo(1)
        scalar1_at_pc_full_level(1,iv)      = half*(scalar1_at_pc_face_level(1,iv)   +scalar1_at_pc_face_level(2,iv))
        scalar1_at_pc_full_level(nlev,iv)   = half*(scalar1_at_pc_face_level(nlev,iv)+scalar1_at_pc_face_level(nlev+1,iv))
        scalar2_at_pc_full_level(1,iv)      = half*(scalar2_at_pc_face_level(1,iv)   +scalar2_at_pc_face_level(2,iv))
        scalar2_at_pc_full_level(nlev,iv)   = half*(scalar2_at_pc_face_level(nlev,iv)+scalar2_at_pc_face_level(nlev+1,iv))
     end do
!$omp end do nowait
!$omp end parallel

    IF(EQS_VERT_DIFF)THEN
    select case(order)
    case(2)
       do iv = 1, mesh%nv_halo(1)
         do ilev = 2, nlev-1
            scalar1_at_pc_full_level(ilev,iv) = half*(scalar1_at_pc_face_level(ilev,iv)+scalar1_at_pc_face_level(ilev+1,iv))
            scalar2_at_pc_full_level(ilev,iv) = half*(scalar2_at_pc_face_level(ilev,iv)+scalar2_at_pc_face_level(ilev+1,iv))
         end do
       end do

    case(3) 
!$omp parallel  private(ii,iv,ilev)      
!$omp do schedule(dynamic,100) 
    do ii=1,mesh%nv_halo(1)*(nlev-2),1
        iv=ceiling(ii/real(nlev-2,r8))
        ilev=1+ii-(iv-1)*(nlev-2)
!       do iv = 1, mesh%nv_halo(1)
!         do ilev = 2, nlev-1
            scalar1_at_pc_full_level(ilev,iv) = (scalar1_at_pc_face_level(ilev,iv)  +scalar1_at_pc_face_level(ilev+1,iv))*(7._r8/12._r8)-&
                                                (scalar1_at_pc_face_level(ilev+2,iv)+scalar1_at_pc_face_level(ilev-1,iv))*(1._r8/12._r8)+&
                                 sign(1._r8,scalar_eta_mass_velocity_at_pc_full(ilev,iv))*&
                                          ((scalar1_at_pc_face_level(ilev+2,iv)-scalar1_at_pc_face_level(ilev-1,iv))-&
                                     3._r8*(scalar1_at_pc_face_level(ilev+1,iv)-scalar1_at_pc_face_level(ilev,iv)))/12._r8

            scalar2_at_pc_full_level(ilev,iv) = (scalar2_at_pc_face_level(ilev,iv)  +scalar2_at_pc_face_level(ilev+1,iv))*(7._r8/12._r8)-&
                                                (scalar2_at_pc_face_level(ilev+2,iv)+scalar2_at_pc_face_level(ilev-1,iv))*(1._r8/12._r8)+&
                                 sign(1._r8,scalar_eta_mass_velocity_at_pc_full(ilev,iv))*&
                                          ((scalar2_at_pc_face_level(ilev+2,iv)-scalar2_at_pc_face_level(ilev-1,iv))-&
                                     3._r8*(scalar2_at_pc_face_level(ilev+1,iv)-scalar2_at_pc_face_level(ilev,iv)))/12._r8

!         end do
       !end do
    end do
!$omp end do nowait
!$omp end parallel 
     case default
        print*," you must select a vert order in calc_hpe_vert_flux_operator_face"
        stop
     end select

     ELSE

     select case(order)
     case(2) 
       do iv = 1, mesh%nv_halo(1)
         do ilev = 2, nlev-1
            scalar1_at_pc_full_level(ilev,iv) = half*(scalar1_at_pc_face_level(ilev,iv)+scalar1_at_pc_face_level(ilev+1,iv))
            scalar2_at_pc_full_level(ilev,iv) = half*(scalar2_at_pc_face_level(ilev,iv)+scalar2_at_pc_face_level(ilev+1,iv))
         end do
       end do
     case(3)
      do iv = 1, mesh%nv_halo(1)
         do ilev = 2, nlev-1
            ! var1
            part1 = half*(scalar1_at_pc_face_level(ilev,iv)+scalar1_at_pc_face_level(ilev+1,iv))

            der_k = 2._r8*(deta_face(ilev)**2)*((scalar1_at_pc_face_level(ilev+1,iv)-scalar1_at_pc_face_level(ilev,iv))/deta_full(ilev)-&
                                                (scalar1_at_pc_face_level(ilev,iv)-scalar1_at_pc_face_level(ilev-1,iv))/deta_full(ilev-1))/&
                                                (deta_full(ilev)+deta_full(ilev-1))

            der_kp1 = 2._r8*(deta_face(ilev+1)**2)*((scalar1_at_pc_face_level(ilev+2,iv)-scalar1_at_pc_face_level(ilev+1,iv))/deta_full(ilev+1)-&
                                                    (scalar1_at_pc_face_level(ilev+1,iv)-scalar1_at_pc_face_level(ilev,iv))/deta_full(ilev))/&
                                                    (deta_full(ilev+1)+deta_full(ilev))

            scalar1_at_pc_full_level(ilev,iv)= part1-(der_k+der_kp1)/12._r8+sign(1._r8,scalar_eta_mass_velocity_at_pc_full(ilev,iv))*(der_kp1-der_k)
         
            ! var2
            part1 = half*(scalar2_at_pc_face_level(ilev,iv)+scalar2_at_pc_face_level(ilev+1,iv))

            der_k = 2._r8*(deta_face(ilev)**2)*((scalar2_at_pc_face_level(ilev+1,iv)-scalar2_at_pc_face_level(ilev,iv))/deta_full(ilev)-&
                                                (scalar2_at_pc_face_level(ilev,iv)  -scalar2_at_pc_face_level(ilev-1,iv))/deta_full(ilev-1))/&
                                                (deta_full(ilev)+deta_full(ilev-1))

            der_kp1 = 2._r8*(deta_face(ilev+1)**2)*((scalar2_at_pc_face_level(ilev+2,iv)-scalar2_at_pc_face_level(ilev+1,iv))/deta_full(ilev+1)-&
                                                    (scalar2_at_pc_face_level(ilev+1,iv)-scalar2_at_pc_face_level(ilev,iv))  /deta_full(ilev))/&
                                                    (deta_full(ilev+1)+deta_full(ilev))

            scalar2_at_pc_full_level(ilev,iv)= part1-(der_k+der_kp1)/12._r8+sign(1._r8,scalar_eta_mass_velocity_at_pc_full(ilev,iv))*(der_kp1-der_k)

         end do
      end do
     case default
        print*," you must select a vert order in calc_hpe_vert_flux_operator_face_2d"
        stop
     end select

     END IF
     
     return
   end subroutine calc_hpe_vert_flux_operator_face_2d_var2

  end module grist_nh_explicit_tend_module_2d
