
!----------------------------------------------------------------------------
! Created on 2016
! Version 1.0
! Description: Incremental Remapping Methods or Flux-Form Semi-Lagrangian
!              1) mr07
!              2) sm10 described in MK13
!              3) tr version mr07
!              4) tr version sm10
! Revision history:
!----------------------------------------------------------------------------

 module grist_ffsl_flux_module

  use grist_domain_types    , only: global_domain
  use grist_constants       , only: i4, r8, rearth
  use grist_math_module     , only: norm
  use grist_recon_module    , only: calc_mr07_area_at_hx, &
                                    calc_sm10_area_at_hx, &
                                    calc_mr07_area_at_tr, &
                                    calc_sm10_area_at_tr

  implicit none

  private

  public ::  calc_hx_normal_flux_mr07 , & ! use b-weight formed by X-axis connects cell0 and 1st nb
             calc_hx_normal_flux_sm10 , &
             calc_tr_normal_flux_mr07 , &
             calc_tr_normal_flux_sm10

  real(r8), parameter  :: one = 1._r8

  contains

!================================================
!       MR07 FLUX AT THE HEXAGONAL EDGE
!================================================

  subroutine calc_hx_normal_flux_mr07(mesh,scalar_normal_velocity_at_edge , &
                                        scalar_tangent_velocity_at_edge, &
                                        scalar_height_at_prime_cell    , &
                                        scalar_normal_flux_at_edge     , &
                                        adv_flag                       , &
                                        dtime)
! io
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge(:)
   real(r8), allocatable, intent(in)    :: scalar_tangent_velocity_at_edge(:)
   real(r8), allocatable, intent(in)    :: scalar_height_at_prime_cell(:)
   real(r8), allocatable, intent(inout) :: scalar_normal_flux_at_edge(:)
   integer(i4),           intent(in)    :: adv_flag
   real(r8)   ,           intent(in)    :: dtime
! local
   real(r8),dimension(1:3)              :: v0v1, v2v3 ! 0->1, 2->3
   real(r8)                             :: nr_wind_edge
   real(r8)                             :: tg_wind_edge
   real(r8)                             :: g1(3), g2(3)
   real(r8)                             :: g1_lob_x
   real(r8)                             :: g1_lob_y
   real(r8)                             :: scalar_upwind
   integer(i4)                          :: v0, v1, v2, v3
   integer(i4)                          :: v_upwind
   integer(i4)                          :: ie
   integer(i4)                          :: flag, flag0, flag1,flag2,flag3

!================================================
! normal wind follows edp's nr positive
! tangent wind follows edp's tg positive
!================================================

     do ie = 1, mesh%ne

        v0            = mesh%edt(ie)%v(1)
        v1            = mesh%edt(ie)%v(2)
        v0v1          = mesh%vtx(v1)%p-mesh%vtx(v0)%p

        v2            = mesh%edp(ie)%v(1)
        v3            = mesh%edp(ie)%v(2)
        v2v3          = mesh%tri(v3)%c%p-mesh%tri(v2)%c%p

        flag          = sign(1._r8,dot_product(v0v1/norm(v0v1),mesh%edp(ie)%nr)) ! v0v1 direction relative to nr
        flag0         = sign(1._r8,dot_product(v2v3/norm(v2v3),mesh%edp(ie)%tg)) ! v2v3 direction relative to tg

        flag1         = sign(1._r8,scalar_normal_velocity_at_edge(ie))           ! sign of normal  wind
        flag3         = sign(1._r8,scalar_tangent_velocity_at_edge(ie))          ! sign of tangent wind
        
        nr_wind_edge  = scalar_normal_velocity_at_edge(ie)
        tg_wind_edge  = scalar_tangent_velocity_at_edge(ie)

!================================================
! check upwind side
! initialize v0 as the upwind point
!================================================

        v_upwind      = v0
        flag2         = 0
        if(flag*flag1.eq.-1)then ! v1 is upwind point
           v_upwind   = v1
           flag2      = 1
        end if

!================================================
! Compute Gauss quadrature points
!================================================

        g2       = (mesh%tri(mesh%edp(ie)%v(1))%c%p+mesh%tri(mesh%edp(ie)%v(2))%c%p)*0.5_r8
        g2       = g2/norm(g2)
        g1       = g2-abs(nr_wind_edge*dtime*0.5_r8)/rearth*v0v1/norm(v0v1)*flag*flag1
        g1       = g1-abs(tg_wind_edge*dtime*0.5_r8)/rearth*v2v3/norm(v2v3)*flag0*flag3
        g1       = g1-mesh%vtx(v_upwind)%p  ! vector points from upwind cell to centroid

        g1_lob_x = dot_product(g1,mesh%vtx(v_upwind)%lob_nx(1,:))
        g1_lob_y = dot_product(g1,mesh%vtx(v_upwind)%lob_ny(1,:))

        call calc_mr07_area_at_hx(mesh, scalar_height_at_prime_cell, v_upwind, &
                                        g1_lob_x, g1_lob_y, scalar_upwind)

        scalar_normal_flux_at_edge(ie) = scalar_upwind ! *scalar_normal_velocity_at_edge(ie)

     end do

   return
  end subroutine calc_hx_normal_flux_mr07

!================================================
!       SM10 FLUX AT THE HEXAGONAL EDGE
!================================================

  subroutine calc_hx_normal_flux_sm10(mesh,scalar_normal_velocity_at_edge , &
                                           scalar_tangent_velocity_at_edge, &
                                           scalar_height_at_prime_cell    , &
                                           scalar_normal_flux_at_edge     , &
                                           adv_flag                       , &
                                           dtime)
! io
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge(:)
   real(r8), allocatable, intent(in)    :: scalar_tangent_velocity_at_edge(:)
   real(r8), allocatable, intent(in)    :: scalar_height_at_prime_cell(:)
   real(r8), allocatable, intent(inout) :: scalar_normal_flux_at_edge(:)
   integer(i4),           intent(in)    :: adv_flag
   real(r8)   ,           intent(in)    :: dtime
! local
   real(r8),dimension(1:3)              :: v0v1, v2v3 ! 0->1, 2->3
   real(r8)                             :: nr_wind_edge
   real(r8)                             :: tg_wind_edge
   real(r8)                             :: edge_point(3),t1_point(3),t2_point(3)
   real(r8)                             :: g1(3),g2(3),g3(3),g4(3),g5(3) ! see ms13's FIG1
   real(r8)                             :: g1_lob_x, g2_lob_x, g3_lob_x, g4_lob_x, g5_lob_x
   real(r8)                             :: g1_lob_y, g2_lob_y, g3_lob_y, g4_lob_y, g5_lob_y
   real(r8)                             :: scalar_upwind
   integer(i4)                          :: v0, v1, v2, v3
   integer(i4)                          :: v_upwind, alpha
   integer(i4)                          :: ie
   integer(i4)                          :: flag, flag0, flag1,flag2,flag3

!================================================
! normal wind follows edp's nr positive
! tangent wind follows edp's tg positive
!================================================

     do ie = 1, mesh%ne

        v0            = mesh%edt(ie)%v(1)
        v1            = mesh%edt(ie)%v(2)
        v0v1          = mesh%vtx(v1)%p-mesh%vtx(v0)%p

        v2            = mesh%edp(ie)%v(1)
        v3            = mesh%edp(ie)%v(2)
        v2v3          = mesh%tri(v3)%c%p-mesh%tri(v2)%c%p

        flag          = sign(1._r8,dot_product(v0v1/norm(v0v1),mesh%edp(ie)%nr)) ! v0v1 direction relative to nr
        flag0         = sign(1._r8,dot_product(v2v3/norm(v2v3),mesh%edp(ie)%tg)) ! v2v3 direction relative to tg

        flag1         = sign(1._r8,scalar_normal_velocity_at_edge(ie))           ! sign of normal  wind
        flag3         = sign(1._r8,scalar_tangent_velocity_at_edge(ie))          ! sign of tangent wind
        
        nr_wind_edge  = scalar_normal_velocity_at_edge(ie)
        tg_wind_edge  = scalar_tangent_velocity_at_edge(ie)

!================================================
! check upwind side
! initialize v0 as the upwind point
!================================================

        v_upwind      = v0
        flag2         = 0
        if(flag*flag1.eq.-1)then ! v1 is upwind point
           v_upwind   = v1
           flag2      = 1
        end if

!================================================
! compute Gauss quadrature points
! alpha=1 recovers de original formulation
!================================================
        alpha      = 1

        t1_point   = mesh%tri(mesh%edp(ie)%v(1))%c%p
        t2_point   = mesh%tri(mesh%edp(ie)%v(2))%c%p
        edge_point =(t1_point+t2_point)*0.5_r8
        edge_point = edge_point/norm(edge_point)

        g2   = edge_point+abs(nr_wind_edge*dtime)*(1-alpha)/rearth*v0v1/norm(v0v1)*flag*flag1 
        g2   =         g2+abs(tg_wind_edge*dtime)*(1-alpha)/rearth*v2v3/norm(v2v3)*flag0*flag3

        g3   = edge_point-abs(nr_wind_edge*dtime)*alpha/rearth*v0v1/norm(v0v1)*flag*flag1
        g3   =         g3-abs(tg_wind_edge*dtime)*alpha/rearth*v2v3/norm(v2v3)*flag0*flag3

        g1   = edge_point+0.5_r8*abs(nr_wind_edge*dtime)*(1._r8-2._r8*alpha)/rearth*v0v1/norm(v0v1)*flag*flag1
        g1   =         g1+0.5_r8*abs(tg_wind_edge*dtime)*(1._r8-2._r8*alpha)/rearth*v2v3/norm(v2v3)*flag0*flag3

        g4   = t1_point+0.5_r8*abs(nr_wind_edge*dtime)*(1._r8-2._r8*alpha)/rearth*v0v1/norm(v0v1)*flag*flag1
        g4   =       g4+0.5_r8*abs(tg_wind_edge*dtime)*(1._r8-2._r8*alpha)/rearth*v2v3/norm(v2v3)*flag0*flag3

        g5   = t2_point+0.5_r8*abs(nr_wind_edge*dtime)*(1._r8-2._r8*alpha)/rearth*v0v1/norm(v0v1)*flag*flag1
        g5   =       g5+0.5_r8*abs(tg_wind_edge*dtime)*(1._r8-2._r8*alpha)/rearth*v2v3/norm(v2v3)*flag0*flag3

! only upwind set

        g1   = g1-mesh%vtx(v_upwind)%p    ! vector points from upwind cell to centroid
        g2   = g2-mesh%vtx(v_upwind)%p    ! vector points from upwind cell to centroid
        g3   = g3-mesh%vtx(v_upwind)%p    ! vector points from upwind cell to centroid
        g4   = g4-mesh%vtx(v_upwind)%p    ! vector points from upwind cell to centroid
        g5   = g5-mesh%vtx(v_upwind)%p    ! vector points from upwind cell to centroid

        g1_lob_x = dot_product(g1,mesh%vtx(v_upwind)%lob_nx(1,:))
        g1_lob_y = dot_product(g1,mesh%vtx(v_upwind)%lob_ny(1,:))
        g2_lob_x = dot_product(g2,mesh%vtx(v_upwind)%lob_nx(1,:))
        g2_lob_y = dot_product(g2,mesh%vtx(v_upwind)%lob_ny(1,:))
        g3_lob_x = dot_product(g3,mesh%vtx(v_upwind)%lob_nx(1,:))
        g3_lob_y = dot_product(g3,mesh%vtx(v_upwind)%lob_ny(1,:))
        g4_lob_x = dot_product(g4,mesh%vtx(v_upwind)%lob_nx(1,:))
        g4_lob_y = dot_product(g4,mesh%vtx(v_upwind)%lob_ny(1,:))
        g5_lob_x = dot_product(g5,mesh%vtx(v_upwind)%lob_nx(1,:))
        g5_lob_y = dot_product(g5,mesh%vtx(v_upwind)%lob_ny(1,:))


        call calc_sm10_area_at_hx(mesh, scalar_height_at_prime_cell, v_upwind, &
                                        g1_lob_x, g1_lob_y, &
                                        g2_lob_x, g2_lob_y, &
                                        g3_lob_x, g3_lob_y, &
                                        g4_lob_x, g4_lob_y, &
                                        g5_lob_x, g5_lob_y, &
                                        scalar_upwind)

        scalar_normal_flux_at_edge(ie) = scalar_upwind !*scalar_normal_velocity_at_edge(ie)

     end do

   return
  end subroutine calc_hx_normal_flux_sm10

!================================================
!       MR07 FLUX AT THE TRIANGULAR EDGE
!================================================

  subroutine calc_tr_normal_flux_mr07(mesh,scalar_normal_velocity_at_edge , &
                                           scalar_tangent_velocity_at_edge, &
                                           scalar_height_at_dual_cell    , &
                                           scalar_normal_flux_at_edge     , &
                                           called_by_tr                   , &
                                           dtime)
! io
      type(global_domain),   intent(in)    :: mesh
      real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge(:)
      real(r8), allocatable, intent(in)    :: scalar_tangent_velocity_at_edge(:)
      real(r8), allocatable, intent(in)    :: scalar_height_at_dual_cell(:)
      real(r8), allocatable, intent(inout) :: scalar_normal_flux_at_edge(:)
      logical,               intent(in)    :: called_by_tr
      real(r8)   ,           intent(in)    :: dtime
! local
      real(r8),dimension(1:3)              :: v0v1, v2v3 ! 0->1, 2->3
      real(r8)                             :: nr_wind_edge
      real(r8)                             :: tg_wind_edge
      real(r8)                             :: g1(3), g2(3)
      real(r8)                             :: g1_lob_x
      real(r8)                             :: g1_lob_y
      real(r8)                             :: scalar_upwind
      integer(i4)                          :: v0, v1, v2, v3
      integer(i4)                          :: v_upwind
      integer(i4)                          :: ie
      integer(i4)                          :: flag, flag0, flag1,flag2,flag3

!================================================
! normal wind follows edp's nr positive
! tangent wind follows edp's tg positive
!================================================

      do ie = 1, mesh%ne

         v0            = mesh%edp(ie)%v(1)
         v1            = mesh%edp(ie)%v(2)
         v0v1          = mesh%tri(v1)%c%p-mesh%tri(v0)%c%p

         v2            = mesh%edt(ie)%v(1)
         v3            = mesh%edt(ie)%v(2)
         v2v3          = mesh%vtx(v3)%p-mesh%vtx(v2)%p

         flag          = sign(1._r8,dot_product(v0v1/norm(v0v1),mesh%edt(ie)%nr)) ! v0v1 direction relative to nr
         flag0         = sign(1._r8,dot_product(v2v3/norm(v2v3),mesh%edt(ie)%tg)) ! v2v3 direction relative to tg

         flag1         = sign(1._r8,scalar_normal_velocity_at_edge(ie))           ! sign of normal  wind
         flag3         = sign(1._r8,scalar_tangent_velocity_at_edge(ie))          ! sign of tangent wind
        
         nr_wind_edge  = scalar_normal_velocity_at_edge(ie)
         tg_wind_edge  = scalar_tangent_velocity_at_edge(ie)

!================================================
! check upwind side
! initialize v0 as the upwind point
!================================================

         v_upwind      = v0
         flag2         = 0
         if(flag*flag1.eq.-1)then ! v1 is upwind point
            v_upwind   = v1
            flag2      = 1
         end if

!================================================
! Compute Gauss quadrature points
!================================================

         g2       = (mesh%vtx(v2)%p+mesh%vtx(v3)%p)*0.5_r8
         g2       = g2/norm(g2)
         g1       = g2-abs(nr_wind_edge*dtime*0.5_r8)/rearth*v0v1/norm(v0v1)*flag*flag1
         g1       = g1-abs(tg_wind_edge*dtime*0.5_r8)/rearth*v2v3/norm(v2v3)*flag0*flag3
         g1       = g1-mesh%tri(v_upwind)%c%p   ! vector points from upwind cell to centroid

         g1_lob_x = dot_product(g1,mesh%tri(v_upwind)%lob_nx(1,:))
         g1_lob_y = dot_product(g1,mesh%tri(v_upwind)%lob_ny(1,:))

         call calc_mr07_area_at_tr(mesh, scalar_height_at_dual_cell, v_upwind, &
                                         g1_lob_x, g1_lob_y, scalar_upwind)

         scalar_normal_flux_at_edge(ie) = scalar_upwind  ! only return q

      end do

     return
  end subroutine calc_tr_normal_flux_mr07

!================================================
!       SM10 FLUX AT THE TRIANGULAR EDGE
!================================================

  subroutine calc_tr_normal_flux_sm10(mesh,scalar_normal_velocity_at_edge , &
                                           scalar_tangent_velocity_at_edge, &
                                           scalar_height_at_dual_cell     , &
                                           scalar_normal_flux_at_edge     , &
                                           called_by_tr                   , &
                                           dtime)
! io
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge(:)
   real(r8), allocatable, intent(in)    :: scalar_tangent_velocity_at_edge(:)
   real(r8), allocatable, intent(in)    :: scalar_height_at_dual_cell(:)
   real(r8), allocatable, intent(inout) :: scalar_normal_flux_at_edge(:)
   logical,               intent(in)    :: called_by_tr
   real(r8),              intent(in)    :: dtime
! local
   real(r8),dimension(1:3)              :: v0v1, v2v3 ! 0->1, 2->3
   real(r8)                             :: nr_wind_edge
   real(r8)                             :: tg_wind_edge
   real(r8)                             :: edge_point(3),t1_point(3),t2_point(3)
   real(r8)                             :: g1(3),g2(3),g3(3),g4(3),g5(3) ! see ms13's FIG1
   real(r8)                             :: g1_lob_x, g2_lob_x, g3_lob_x, g4_lob_x, g5_lob_x
   real(r8)                             :: g1_lob_y, g2_lob_y, g3_lob_y, g4_lob_y, g5_lob_y
   real(r8)                             :: scalar_upwind
   integer(i4)                          :: v0, v1, v2, v3
   integer(i4)                          :: v_upwind
   integer(i4)                          :: ie
   integer(i4)                          :: flag, flag0, flag1,flag2,flag3

!================================================
! normal wind follows edt's nr positive
! tangent wind follows edt's tg positive
!================================================

     do ie = 1, mesh%ne

        v0            = mesh%edp(ie)%v(1)
        v1            = mesh%edp(ie)%v(2)
        v0v1          = mesh%tri(v1)%c%p-mesh%tri(v0)%c%p

        v2            = mesh%edt(ie)%v(1)
        v3            = mesh%edt(ie)%v(2)
        v2v3          = mesh%vtx(v3)%p-mesh%vtx(v2)%p

        if(called_by_tr)then
           flag          = sign(1._r8,dot_product(v0v1/norm(v0v1),mesh%edt(ie)%nr)) ! v0v1 direction relative to edt's nr
           flag0         = sign(1._r8,dot_product(v2v3/norm(v2v3),mesh%edt(ie)%tg)) ! v2v3 direction relative to edt's tg
        else
           flag          = sign(1._r8,dot_product(v0v1/norm(v0v1),mesh%edp(ie)%tg)) ! v0v1 direction relative to edt's nr
           flag0         = sign(1._r8,dot_product(v2v3/norm(v2v3),mesh%edp(ie)%nr)) ! v2v3 direction relative to edt's tg
        end if

        flag1         = sign(1._r8,scalar_normal_velocity_at_edge(ie))           ! sign of normal  wind
        flag3         = sign(1._r8,scalar_tangent_velocity_at_edge(ie))          ! sign of tangent wind
        
        nr_wind_edge  = scalar_normal_velocity_at_edge(ie)
        tg_wind_edge  = scalar_tangent_velocity_at_edge(ie)

!================================================
! check upwind side
! initialize v0 as the upwind point
!================================================

        v_upwind      = v0
        flag2         = 0
        if(flag*flag1.eq.-1)then ! v1 is upwind point
           v_upwind   = v1
           flag2      = 1
        end if

!================================================
! compute Gauss quadrature points
!================================================

! only upwind set

        g2   = (mesh%vtx(v2)%p+mesh%vtx(v3)%p)*0.5_r8
        g2   = g2/norm(g2)

        g1   = g2-abs(nr_wind_edge*dtime*0.5_r8)/rearth*v0v1/norm(v0v1)*flag*flag1
        g1   = g1-abs(tg_wind_edge*dtime*0.5_r8)/rearth*v2v3/norm(v2v3)*flag0*flag3
        g3   = g2-abs(nr_wind_edge*dtime)/rearth*v0v1/norm(v0v1)*flag*flag1
        g3   = g3-abs(tg_wind_edge*dtime)/rearth*v2v3/norm(v2v3)*flag0*flag3
        
        g4   = mesh%vtx(mesh%edt(ie)%v(1))%p-abs(nr_wind_edge*dtime*0.5_r8)/rearth*v0v1/norm(v0v1)*flag*flag1
        g4   = g4-abs(tg_wind_edge*dtime*0.5_r8)/rearth*v2v3/norm(v2v3)*flag0*flag3

        g5   = mesh%vtx(mesh%edt(ie)%v(2))%p-abs(nr_wind_edge*dtime*0.5_r8)/rearth*v0v1/norm(v0v1)*flag*flag1
        g5   = g5-abs(tg_wind_edge*dtime*0.5_r8)/rearth*v2v3/norm(v2v3)*flag0*flag3

! only upwind set

        g1   = g1-mesh%tri(v_upwind)%c%p    ! vector points from upwind cell to centroid
        g2   = g2-mesh%tri(v_upwind)%c%p    ! vector points from upwind cell to centroid
        g3   = g3-mesh%tri(v_upwind)%c%p    ! vector points from upwind cell to centroid
        g4   = g4-mesh%tri(v_upwind)%c%p    ! vector points from upwind cell to centroid
        g5   = g5-mesh%tri(v_upwind)%c%p    ! vector points from upwind cell to centroid

        g1_lob_x = dot_product(g1,mesh%tri(v_upwind)%lob_nx(1,:))
        g1_lob_y = dot_product(g1,mesh%tri(v_upwind)%lob_ny(1,:))
        g2_lob_x = dot_product(g2,mesh%tri(v_upwind)%lob_nx(1,:))
        g2_lob_y = dot_product(g2,mesh%tri(v_upwind)%lob_ny(1,:))
        g3_lob_x = dot_product(g3,mesh%tri(v_upwind)%lob_nx(1,:))
        g3_lob_y = dot_product(g3,mesh%tri(v_upwind)%lob_ny(1,:))
        g4_lob_x = dot_product(g4,mesh%tri(v_upwind)%lob_nx(1,:))
        g4_lob_y = dot_product(g4,mesh%tri(v_upwind)%lob_ny(1,:))
        g5_lob_x = dot_product(g5,mesh%tri(v_upwind)%lob_nx(1,:))
        g5_lob_y = dot_product(g5,mesh%tri(v_upwind)%lob_ny(1,:))


        call calc_sm10_area_at_tr(mesh, scalar_height_at_dual_cell, v_upwind, &
                                        g1_lob_x, g1_lob_y, &
                                        g2_lob_x, g2_lob_y, &
                                        g3_lob_x, g3_lob_y, &
                                        g4_lob_x, g4_lob_y, &
                                        g5_lob_x, g5_lob_y, &
                                        scalar_upwind)

        scalar_normal_flux_at_edge(ie) = scalar_upwind

     end do

   return
  end subroutine calc_tr_normal_flux_sm10

  end module grist_ffsl_flux_module
