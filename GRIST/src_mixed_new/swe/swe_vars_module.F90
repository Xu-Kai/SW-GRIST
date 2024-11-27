
!----------------------------------------------------------------------------
! Created on 2016
! Version 1.0
! Description: All edge points are intersecting points
! Revision history: 
!----------------------------------------------------------------------------

   module swe_vars_module

   use grist_constants,      only: r8, pi, rearth, i4
   use grist_domain_types,   only: global_domain
   use grist_data_types,     only: scalar_1d_field

   implicit none

   public

!================================================
!                Prognostic
!================================================

! vector component
   type(scalar_1d_field)   :: tend_normal_velocity_at_edge_nm2     ! n-2
   type(scalar_1d_field)   :: tend_normal_velocity_at_edge_nm1     ! n-1
   type(scalar_1d_field)   :: tend_normal_velocity_at_edge         ! n
   type(scalar_1d_field)   :: scalar_normal_velocity_at_edge_nm2   ! n-2
   type(scalar_1d_field)   :: scalar_normal_velocity_at_edge_nm1   ! n-1
   type(scalar_1d_field)   :: scalar_normal_velocity_at_edge       ! n
   type(scalar_1d_field)   :: scalar_tangent_velocity_at_edge      ! n
   type(scalar_1d_field)   :: scalar_normal_flux_at_edge           ! n
   type(scalar_1d_field)   :: scalar_normal_velocity_at_edge_next  ! n+1
! scalar
   type(scalar_1d_field)   :: tend_height_at_prime_cell_nm2        ! n-2
   type(scalar_1d_field)   :: tend_height_at_prime_cell_nm1        ! n-1
   type(scalar_1d_field)   :: tend_height_at_prime_cell            ! n
   type(scalar_1d_field)   :: scalar_height_at_prime_cell_nm2      ! n-2
   type(scalar_1d_field)   :: scalar_height_at_prime_cell_nm1      ! n-1
   type(scalar_1d_field)   :: scalar_height_at_prime_cell          ! n
   type(scalar_1d_field)   :: scalar_height_at_prime_cell_next     ! n+1

!================================================
!                  Diagnostic (13)
!================================================

   type(scalar_1d_field)   :: ue_init                              ! Initial U at edge
   type(scalar_1d_field)   :: ve_init                              ! Initial V at edge
   type(scalar_1d_field)   :: up_init                              ! Initial U at prime cell
   type(scalar_1d_field)   :: vp_init                              ! Initial V at prime cell

   type(scalar_1d_field)   :: ue_reconst                           ! Reconstructed U at edge
   type(scalar_1d_field)   :: ve_reconst                           ! Reconstructed V at edge
   type(scalar_1d_field)   :: up_reconst                           ! Reconstructed U at prime cell
   type(scalar_1d_field)   :: vp_reconst                           ! Reconstructed V at prime cell

   type(scalar_1d_field)   :: scalar_hb_at_prime_cell                 ! diagnosed, fluid+topo
   type(scalar_1d_field)   :: scalar_divergence_at_prime_cell         ! diagnosed
   type(scalar_1d_field)   :: scalar_absolute_vorticity_at_dual_cell  ! diagnosed
   type(scalar_1d_field)   :: scalar_relative_vorticity_at_dual_cell  ! diagnosed
   type(scalar_1d_field)   :: scalar_potential_vorticity_at_dual_cell ! diagnosed

!================================================
!                  Constant (5)
!================================================

   type(scalar_1d_field)   :: scalar_topo_at_prime_cell               ! constant
   type(scalar_1d_field)   :: scalar_area_of_prime_cell               ! constant
   type(scalar_1d_field)   :: scalar_area_of_dual_cell                ! constant
   type(scalar_1d_field)   :: scalar_leng_of_hx_edge                  ! constant
   type(scalar_1d_field)   :: scalar_leng_of_tr_edge                  ! constant
   type(scalar_1d_field)   :: lon_nt
   type(scalar_1d_field)   :: lat_nt
   type(scalar_1d_field)   :: lon_ne
   type(scalar_1d_field)   :: lat_ne
   type(scalar_1d_field)   :: lon_nv
   type(scalar_1d_field)   :: lat_nv

!================================================
!    Vars for testing advection on a Dual Cell
!================================================

   type(scalar_1d_field)   :: scalar_normal_velocity_at_tr_edge
   type(scalar_1d_field)   :: scalar_height_at_dual_cell
   type(scalar_1d_field)   :: scalar_normal_velocity_at_tr_edge_next
   type(scalar_1d_field)   :: scalar_height_at_dual_cell_next

    contains

    subroutine swe_vars_init(mesh)
     type(global_domain), intent(in) :: mesh
! local
     integer(i4)   :: it
     integer(i4)   :: ie
     integer(i4)   :: iv

!================================================
!                   Prognostic
!================================================

      allocate(tend_normal_velocity_at_edge_nm2%f(mesh%ne))
      allocate(tend_normal_velocity_at_edge_nm1%f(mesh%ne))
      allocate(tend_normal_velocity_at_edge%f(mesh%ne))
      allocate(scalar_normal_velocity_at_edge_nm2%f(mesh%ne))
      allocate(scalar_normal_velocity_at_edge_nm1%f(mesh%ne))
      allocate(scalar_normal_velocity_at_edge%f(mesh%ne))
      allocate(scalar_tangent_velocity_at_edge%f(mesh%ne))
      allocate(scalar_normal_flux_at_edge%f(mesh%ne))
      allocate(scalar_normal_velocity_at_edge_next%f(mesh%ne))

      allocate(tend_height_at_prime_cell_nm2%f(mesh%nv))
      allocate(tend_height_at_prime_cell_nm1%f(mesh%nv))
      allocate(tend_height_at_prime_cell%f(mesh%nv))
      allocate(scalar_height_at_prime_cell_nm2%f(mesh%nv))
      allocate(scalar_height_at_prime_cell_nm1%f(mesh%nv))
      allocate(scalar_height_at_prime_cell%f(mesh%nv))
      allocate(scalar_height_at_prime_cell_next%f(mesh%nv))

      tend_normal_velocity_at_edge_nm2%pos    = 6
      tend_normal_velocity_at_edge_nm1%pos    = 6
      tend_normal_velocity_at_edge%pos        = 6
      scalar_normal_velocity_at_edge_nm2%pos  = 6
      scalar_normal_velocity_at_edge_nm1%pos  = 6
      scalar_normal_velocity_at_edge%pos      = 6
      scalar_tangent_velocity_at_edge%pos     = 6
      scalar_normal_flux_at_edge%pos          = 6
      scalar_normal_velocity_at_edge_next%pos = 6

      tend_height_at_prime_cell_nm2%pos       = 0
      tend_height_at_prime_cell_nm1%pos       = 0
      tend_height_at_prime_cell%pos           = 0
      scalar_height_at_prime_cell_nm2%pos     = 0
      scalar_height_at_prime_cell_nm1%pos     = 0
      scalar_height_at_prime_cell%pos         = 0
      scalar_height_at_prime_cell_next%pos    = 0

!================================================
!              Diagnostic
!================================================

      allocate(ue_init%f(mesh%ne))
      allocate(ve_init%f(mesh%ne))
      allocate(up_init%f(mesh%nv))
      allocate(vp_init%f(mesh%nv))

      allocate(ue_reconst%f(mesh%ne))
      allocate(ve_reconst%f(mesh%ne))
      allocate(up_reconst%f(mesh%nv))
      allocate(vp_reconst%f(mesh%nv))

      allocate(scalar_hb_at_prime_cell%f(mesh%nv))
      allocate(scalar_divergence_at_prime_cell%f(mesh%nv))
      allocate(scalar_absolute_vorticity_at_dual_cell%f(mesh%nt))
      allocate(scalar_relative_vorticity_at_dual_cell%f(mesh%nt))
      allocate(scalar_potential_vorticity_at_dual_cell%f(mesh%nt))

      ue_init%pos                         = 6
      ve_init%pos                         = 6
      ue_reconst%pos                      = 6
      ve_reconst%pos                      = 6

      up_init%pos                         = 0
      vp_init%pos                         = 0
      up_reconst%pos                      = 0
      vp_reconst%pos                      = 0

      scalar_hb_at_prime_cell%pos                 = 0
      scalar_divergence_at_prime_cell%pos         = 0
      scalar_absolute_vorticity_at_dual_cell%pos  = 1
      scalar_relative_vorticity_at_dual_cell%pos  = 1
      scalar_potential_vorticity_at_dual_cell%pos = 1

!================================================
!                    Constant
!================================================

     allocate(scalar_topo_at_prime_cell%f(mesh%nv))
     allocate(scalar_area_of_prime_cell%f(mesh%nv))
     allocate(scalar_area_of_dual_cell%f(mesh%nt))
     allocate(scalar_leng_of_hx_edge%f(mesh%ne))
     allocate(scalar_leng_of_tr_edge%f(mesh%ne))
     allocate(lon_nt%f(mesh%nt))
     allocate(lat_nt%f(mesh%nt))
     allocate(lon_ne%f(mesh%ne))
     allocate(lat_ne%f(mesh%ne))
     allocate(lon_nv%f(mesh%nv))
     allocate(lat_nv%f(mesh%nv))

     scalar_topo_at_prime_cell%pos      = 0
     scalar_area_of_prime_cell%pos      = 0
     scalar_area_of_dual_cell%pos       = 1
     scalar_leng_of_hx_edge%pos         = 6
     scalar_leng_of_tr_edge%pos         = 6
     lon_nt%pos                         = 1
     lat_nt%pos                         = 1
     lon_ne%pos                         = 6
     lat_ne%pos                         = 6
     lon_nv%pos                         = 0
     lat_nv%pos                         = 0

     do it = 1, mesh%nt
        scalar_area_of_dual_cell%f(it)  =(rearth**2)*mesh%tri(it)%areag
        lon_nt%f(it)                    = mesh%tri(it)%c%lon*180._r8/pi
        lat_nt%f(it)                    = mesh%tri(it)%c%lat*180._r8/pi
     end do

     do ie = 1, mesh%ne
        scalar_leng_of_hx_edge%f(ie)    = rearth*mesh%edp(ie)%leng
        scalar_leng_of_tr_edge%f(ie)    = rearth*mesh%edt(ie)%leng
        lon_ne%f(ie)                    = mesh%edp(ie)%c%lon*180._r8/pi
        lat_ne%f(ie)                    = mesh%edp(ie)%c%lat*180._r8/pi
     end do

     do iv = 1, mesh%nv
        scalar_area_of_prime_cell%f(iv) =(rearth**2)*mesh%plg(iv)%areag
        lon_nv%f(iv)                    = mesh%vtx(iv)%lon*180._r8/pi
        lat_nv%f(iv)                    = mesh%vtx(iv)%lat*180._r8/pi
     end do

!================================================
!                 Test Dual cell
!================================================

     allocate(scalar_normal_velocity_at_tr_edge%f(1:mesh%ne))
     allocate(scalar_height_at_dual_cell%f(1:mesh%nt))
     allocate(scalar_normal_velocity_at_tr_edge_next%f(1:mesh%ne))
     allocate(scalar_height_at_dual_cell_next%f(1:mesh%nt))
     
     scalar_normal_velocity_at_tr_edge%pos      = 6
     scalar_height_at_dual_cell%pos             = 1
     scalar_normal_velocity_at_tr_edge_next%pos = 6
     scalar_height_at_dual_cell_next%pos        = 1

     return
    end subroutine swe_vars_init

    subroutine swe_vars_clean()
! prog
      deallocate(tend_normal_velocity_at_edge_nm2%f)
      deallocate(tend_normal_velocity_at_edge_nm1%f)
      deallocate(tend_normal_velocity_at_edge%f)
      deallocate(scalar_normal_velocity_at_edge_nm2%f)
      deallocate(scalar_normal_velocity_at_edge_nm1%f)
      deallocate(scalar_normal_velocity_at_edge%f)
      deallocate(scalar_tangent_velocity_at_edge%f)
      deallocate(scalar_normal_flux_at_edge%f)
      deallocate(scalar_normal_velocity_at_edge_next%f)
      deallocate(tend_height_at_prime_cell_nm2%f)
      deallocate(tend_height_at_prime_cell_nm1%f)
      deallocate(tend_height_at_prime_cell%f)
      deallocate(scalar_height_at_prime_cell_nm2%f)
      deallocate(scalar_height_at_prime_cell_nm1%f)
      deallocate(scalar_height_at_prime_cell%f)
      deallocate(scalar_height_at_prime_cell_next%f)
! diag
      deallocate(ue_init%f)
      deallocate(ve_init%f)
      deallocate(ue_reconst%f)
      deallocate(ve_reconst%f)
      deallocate(up_init%f)
      deallocate(vp_init%f)
      deallocate(up_reconst%f)
      deallocate(vp_reconst%f)
      deallocate(scalar_hb_at_prime_cell%f)
      deallocate(scalar_divergence_at_prime_cell%f)
      deallocate(scalar_absolute_vorticity_at_dual_cell%f)
      deallocate(scalar_relative_vorticity_at_dual_cell%f)
      deallocate(scalar_potential_vorticity_at_dual_cell%f)
! cons
      deallocate(scalar_topo_at_prime_cell%f)
      deallocate(scalar_area_of_prime_cell%f)
      deallocate(scalar_area_of_dual_cell%f)
      deallocate(scalar_leng_of_hx_edge%f)
      deallocate(scalar_leng_of_tr_edge%f)
      deallocate(lon_nt%f)
      deallocate(lat_nt%f)
      deallocate(lon_ne%f)
      deallocate(lat_ne%f)
      deallocate(lon_nv%f)
      deallocate(lat_nv%f)
! dual
      deallocate(scalar_normal_velocity_at_tr_edge%f)
      deallocate(scalar_height_at_dual_cell%f)
      deallocate(scalar_normal_velocity_at_tr_edge_next%f)
      deallocate(scalar_height_at_dual_cell_next%f)

      return
    end subroutine swe_vars_clean

   end module swe_vars_module
