
!----------------------------------------------------------------------------
! Created on 2016
! Version 1.0
! Description: Module for computation needed for SWE
! Revision history:
!----------------------------------------------------------------------------

 module grist_hori_swe_module

   use grist_domain_types  , only: global_domain
   use grist_constants,      only: i4, r8, rearth
   use grist_nml_module,     only: pv_order, mass_wtd_vor,stencil_width
   use grist_flux_operators, only: calc_normal_flux_tr_edge
   use grist_recon_module,   only: calc_coriolis_term, calc_scalar_at_dual_cell
   use grist_gcd_module,     only: calc_vorticity_at_dual_cell

   implicit none

   private

   public :: calc_tend_velocity_coriolis   , &
             calc_kinetic_energy

   contains

  subroutine calc_tend_velocity_coriolis(mesh,scalar_height_at_prime_cell   ,&
                                              scalar_normal_velocity_at_edge,&
                                              scalar_normal_flux_at_edge    ,&
                                              tend_velocity_coriolis        ,&
                                              dtime, irk_step)
!io
   type(global_domain),   intent(inout) :: mesh
   real(r8), allocatable, intent(in)    :: scalar_height_at_prime_cell(:)    ! initial condition
   real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge(:) ! initial condition
   real(r8), allocatable, intent(in)    :: scalar_normal_flux_at_edge(:)     ! initial condition
   real(r8), allocatable, intent(inout) :: tend_velocity_coriolis(:)         ! what we send back 
   real(r8)             , intent(in)    :: dtime
   integer(i4)          , intent(in)    :: irk_step
!local
   real(r8), allocatable                :: scalar_height_at_dual_cell(:)
   real(r8), allocatable                :: scalar_absolute_vorticity_at_dual_cell(:)
   real(r8), allocatable                :: scalar_relative_vorticity_at_dual_cell(:)
   real(r8), allocatable                :: scalar_pv_at_dual_cell(:)
   real(r8), allocatable                :: scalar_pv_at_edge(:)
   real(r8), allocatable                :: scalar_tangent_velocity_at_edge(:)
   real(r8), allocatable                :: scalar_coriolis_term_at_edge(:)   ! what we need
   real(r8), allocatable                :: tmp(:)
   integer(i4)                          :: it
   integer(i4)                          :: ie
   integer(i4)                          :: iv

   if(.not.allocated(scalar_height_at_dual_cell))      allocate(scalar_height_at_dual_cell(1:mesh%nt_full))
   if(.not.allocated(scalar_absolute_vorticity_at_dual_cell))&
                                                       allocate(scalar_absolute_vorticity_at_dual_cell(1:mesh%nt_full))
   if(.not.allocated(scalar_relative_vorticity_at_dual_cell))&
                                                       allocate(scalar_relative_vorticity_at_dual_cell(1:mesh%nt_full))
   if(.not.allocated(scalar_pv_at_edge))               allocate(scalar_pv_at_edge(1:mesh%ne_full))
   if(.not.allocated(scalar_pv_at_dual_cell))          allocate(scalar_pv_at_dual_cell(1:mesh%nt_full))
   if(.not.allocated(scalar_coriolis_term_at_edge))    allocate(scalar_coriolis_term_at_edge(1:mesh%ne_full))
   if(.not.allocated(tmp))                             allocate(tmp(1:mesh%ne_full))
   if(.not.allocated(scalar_tangent_velocity_at_edge)) allocate(scalar_tangent_velocity_at_edge(1:mesh%ne_full))

! calculate vorticity at dual cell using Stokes' Theorem
     mesh%ne = mesh%ne_halo(2)
     mesh%nt = mesh%nt_halo(2)
     call calc_vorticity_at_dual_cell(mesh, scalar_normal_velocity_at_edge        ,&
                                            scalar_absolute_vorticity_at_dual_cell,&
                                            scalar_relative_vorticity_at_dual_cell) ! for diagnose

     if(mass_wtd_vor)then
        call calc_scalar_at_dual_cell(mesh,scalar_height_at_prime_cell, scalar_height_at_dual_cell)
        scalar_pv_at_dual_cell = scalar_absolute_vorticity_at_dual_cell/scalar_height_at_dual_cell
     else
        scalar_pv_at_dual_cell = scalar_absolute_vorticity_at_dual_cell
     end if
!
! reconstruct tangent wind
!
     tmp(:) = 1._r8
     call calc_coriolis_term(mesh,tmp,scalar_normal_velocity_at_edge, scalar_tangent_velocity_at_edge)
     mesh%ne = mesh%ne_compute
     mesh%nt = mesh%nt_compute
!
! tangent flux along hx's edge
!
     mesh%ne = mesh%ne_halo(1)
     mesh%nt = mesh%nt_halo(1)
     call calc_normal_flux_tr_edge(mesh, scalar_pv_at_dual_cell         ,&
                                         scalar_tangent_velocity_at_edge,&
                                         scalar_normal_velocity_at_edge ,&
                                         scalar_pv_at_edge              ,&
                                         .false.                        ,&
                                         pv_order(irk_step), dtime)

     mesh%ne = mesh%ne_compute
     mesh%nt = mesh%nt_compute

!================================================
!       The Final nonlinear Coriolis tendency
!================================================

     if(mass_wtd_vor)then
        call calc_coriolis_term(mesh,scalar_pv_at_edge,&
                                     scalar_normal_flux_at_edge,&
                                     scalar_coriolis_term_at_edge)

     else
        call calc_coriolis_term(mesh,scalar_pv_at_edge,&
                                     scalar_normal_velocity_at_edge,&
                                     scalar_coriolis_term_at_edge)
     end if

! send back
      tend_velocity_coriolis(:)  = -1._r8*scalar_coriolis_term_at_edge(:)

! clean 
      deallocate(scalar_height_at_dual_cell)
      deallocate(scalar_absolute_vorticity_at_dual_cell)
      deallocate(scalar_relative_vorticity_at_dual_cell)
      deallocate(scalar_pv_at_edge)
      deallocate(scalar_pv_at_dual_cell)
      deallocate(scalar_coriolis_term_at_edge)
      deallocate(tmp)
      deallocate(scalar_tangent_velocity_at_edge)

    return
  end subroutine calc_tend_velocity_coriolis

!==============================================
!   Calculate kinetic energy at prime cell
!==============================================

  subroutine calc_kinetic_energy(mesh,scalar_normal_velocity_at_edge,kinetic_energy)
! use
  use grist_nml_module,  only: ke_method, ke_alpha
! io
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge(:)
   real(r8), allocatable, intent(inout) :: kinetic_energy(:)      ! at cell
! local
   real(r8), allocatable                :: kinetic_energy_hx(:)   ! at hx
   real(r8), allocatable                :: kinetic_energy_tr(:)   ! at tr
   real(r8)                             :: length_of_voronoi_edge
   real(r8)                             :: length_of_triangle_edge
   real(r8)                             :: tmp_sum
   real(r8)                             :: wind
   integer                              :: i, it, ie, iv
   integer(i4)                          :: index_edge

   select case(ke_method)

   case (0)   ! R10's formulation

        !do iv = 1, mesh%nv
        !   tmp_sum  = 0._r8
        !   do ie = 1, mesh%vtx(iv)%nnb
        !      index_edge              = mesh%vtx(iv)%ed(ie)   ! global index of this edge pair 
        !      length_of_triangle_edge = rearth*mesh%edt(index_edge)%leng
        !      length_of_voronoi_edge  = rearth*mesh%edp(index_edge)%leng
        !      wind                    = scalar_normal_velocity_at_edge(index_edge)
        !      tmp_sum                 = tmp_sum+0.25_r8*length_of_voronoi_edge*length_of_triangle_edge*wind**2
        !   end do
        !   kinetic_energy(iv) = tmp_sum/((rearth**2)*mesh%plg(iv)%areag)
        !end do

        do iv = 1, mesh%nv
           tmp_sum  = 0._r8
           do ie = 1, mesh%vtx_nnb(iv)
              index_edge              = mesh%vtx_ed(ie, iv)   ! global index of this edge pair 
              length_of_triangle_edge = rearth*mesh%edt_leng(index_edge)
              length_of_voronoi_edge  = rearth*mesh%edp_leng(index_edge)
              wind                    = scalar_normal_velocity_at_edge(index_edge)
              tmp_sum                 = tmp_sum+0.25_r8*length_of_voronoi_edge*length_of_triangle_edge*wind**2
           end do
           kinetic_energy(iv) = tmp_sum/((rearth**2)*mesh%plg_areag(iv))
        end do
   case(1)   ! s12/g13's formulation
!
! r10's formulation, again
!
        !do iv = 1, mesh%nv
        !   tmp_sum  = 0._r8
        !   do ie = 1, mesh%vtx(iv)%nnb
        !      index_edge              = mesh%vtx(iv)%ed(ie)   ! global index of this edge pair 
        !      length_of_triangle_edge = rearth*mesh%edt(index_edge)%leng
        !      length_of_voronoi_edge  = rearth*mesh%edp(index_edge)%leng
        !      wind                    = scalar_normal_velocity_at_edge(index_edge)
        !      tmp_sum                 = tmp_sum+0.25_r8*length_of_voronoi_edge*length_of_triangle_edge*wind**2
        !   end do
        !   kinetic_energy(iv) = tmp_sum/((rearth**2)*mesh%plg(iv)%areag)
        !end do

        do iv = 1, mesh%nv
           tmp_sum  = 0._r8
           do ie = 1, mesh%vtx_nnb(iv)
              index_edge              = mesh%vtx_ed(ie, iv)   ! global index of this edge pair 
              length_of_triangle_edge = rearth*mesh%edt_leng(index_edge)
              length_of_voronoi_edge  = rearth*mesh%edp_leng(index_edge)
              wind                    = scalar_normal_velocity_at_edge(index_edge)
              tmp_sum                 = tmp_sum+0.25_r8*length_of_voronoi_edge*length_of_triangle_edge*wind**2
           end do
           kinetic_energy(iv) = tmp_sum/((rearth**2)*mesh%plg_areag(iv))
        end do

       if(.not.allocated(kinetic_energy_hx)) allocate(kinetic_energy_hx(mesh%nv))
       if(.not.allocated(kinetic_energy_tr)) allocate(kinetic_energy_tr(mesh%nt))
!
! compute ke at cell vertex/triangular center
!
        !do it = 1, mesh%nt
        !   tmp_sum  = 0._r8
        !   do ie = 1, 3
        !      index_edge              = mesh%tri(it)%ed(ie)
        !      length_of_voronoi_edge  = mesh%edp(index_edge)%leng
        !      length_of_triangle_edge = mesh%edt(index_edge)%leng
        !      wind                    = scalar_normal_velocity_at_edge(index_edge)
        !      tmp_sum                 = tmp_sum + 0.25_r8*length_of_voronoi_edge*length_of_triangle_edge*(wind**2)
        !   end do
        !      kinetic_energy_tr(it) = tmp_sum/mesh%tri(it)%areag
        !end do 
        do it = 1, mesh%nt
           tmp_sum  = 0._r8
           do ie = 1, mesh%tri_nnb(it) 
              index_edge              = mesh%tri_ed(ie, it)
              length_of_voronoi_edge  = mesh%edp_leng(index_edge)
              length_of_triangle_edge = mesh%edt_leng(index_edge)
              wind                    = scalar_normal_velocity_at_edge(index_edge)
              tmp_sum                 = tmp_sum + 0.25_r8*length_of_voronoi_edge*length_of_triangle_edge*(wind**2)
           end do
              kinetic_energy_tr(it) = tmp_sum/mesh%tri_areag(it)
        end do 
!
! remap ke from cell vertex to cell center
!
        !do iv = 1, mesh%nv
        !   tmp_sum = 0._r8
        !   do ie = 1, mesh%vtx(iv)%nnb  ! local index of each corner/tr
        !      tmp_sum = tmp_sum+kinetic_energy_tr(mesh%vtx(iv)%tr(ie))*mesh%plg(iv)%kite_area(ie)
        !   end do
        !   kinetic_energy_hx(iv) = tmp_sum
! merge
        !   kinetic_energy(iv) = ke_alpha*kinetic_energy(iv)+(1._r8-ke_alpha)*kinetic_energy_hx(iv)
        !end do
        do iv = 1, mesh%nv
           tmp_sum = 0._r8
           do ie = 1, mesh%vtx_nnb(iv)  ! local index of each corner/tr
              tmp_sum = tmp_sum+kinetic_energy_tr(mesh%vtx_tr(ie, iv))*mesh%plg_kite_area(ie, iv)
           end do
           kinetic_energy_hx(iv) = tmp_sum
! merge
           kinetic_energy(iv) = ke_alpha*kinetic_energy(iv)+(1._r8-ke_alpha)*kinetic_energy_hx(iv)
        end do

        deallocate(kinetic_energy_hx)
        deallocate(kinetic_energy_tr)

   case default

       print*, "you must select a ke_method, stop"
       stop
   end select

   return
  end subroutine calc_kinetic_energy

  end module grist_hori_swe_module
