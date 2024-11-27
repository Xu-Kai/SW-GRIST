
!----------------------------------------------------------------------------
! Created on 2016
! Version 1.0
! Description: Module for reconsructing:
!              1) derivative at cell point (polynomial)
!              2) scalar value at edge point (polynomial)
!              3) tangent coriolis term (vector recons)
!              4) remap from hx to tr
!              5) U, V recovered
! Revision history:
!----------------------------------------------------------------------------

 module grist_recon_module

  use grist_constants,      only: i4, r8, rearth
  use grist_domain_types  , only: global_domain
  use grist_math_module,    only: convert_vector_cart2sph,norm,cross_product
  use grist_mesh_weight_icosh,only: project_sphere2tplane

  implicit none

  private

  public :: calc_der_1st_at_hx  ,&
            calc_double_der_2nd_at_hx, &
            calc_der_2nd_at_hx  ,&
            calc_der_4th_at_hx  ,&
            calc_der_2nd_at_tr  ,&
            calc_der_4th_at_tr  ,&
            calc_mid_2nd_at_hx  ,&
            calc_mid_2nd_at_tr  ,&
            calc_mr07_area_at_hx,&
            calc_sm10_area_at_hx,&
            calc_mr07_area_at_tr,&
            calc_sm10_area_at_tr,&
            calc_coriolis_term  ,&
            project_uv          ,&
            vector_recon_perot_edge2cell_uv,  &
            vector_recon_perot_edge2trcell_uv,&
            calc_scalar_at_dual_cell

  contains

   subroutine calc_der_1st_at_hx(mesh,scalar_height_at_prime_cell,iv,ie,flag,der)
! io
    type(global_domain),  intent(in) :: mesh
    real(r8),allocatable, intent(in) :: scalar_height_at_prime_cell(:)
    integer(i4)         , intent(in) :: iv
    integer(i4)         , intent(in) :: ie
    integer(i4)         , intent(in) :: flag  ! 0 or 1
    real(r8)            , intent(out):: der
! local
    real(r8)            , allocatable :: s_array(:)
    real(r8)                          :: f_array

      allocate(s_array(1:mesh%plg(iv)%stencil_number_2nd))
      s_array(:) = scalar_height_at_prime_cell(mesh%plg(iv)%stencil_index_2nd(:))
      s_array(:) = s_array(:)-scalar_height_at_prime_cell(iv)
      f_array    = dot_product(mesh%edt(ie)%v_weight_1st_prime_cell(flag+1,1:mesh%plg(iv)%stencil_number_2nd), s_array(:))
      der        = f_array

      deallocate(s_array)
      return
   end subroutine calc_der_1st_at_hx

   subroutine calc_double_der_2nd_at_hx(mesh,local_der_2nd_center,local_der_2nd,nnb,iv,ie,flag,der)
! io
    type(global_domain),  intent(in)  :: mesh
    real(r8),             intent(in)  :: local_der_2nd_center
    real(r8),             intent(in)  :: local_der_2nd(10)
    integer(i4)         , intent(in)  :: nnb
    integer(i4)         , intent(in)  :: iv, ie
    integer(i4)         , intent(in)  :: flag  ! 0 or 1
    real(r8)            , intent(out) :: der
! local
    real(r8)            , allocatable :: s_array(:)
    real(r8)                          :: f_array

      !allocate(s_array(1:mesh%plg(iv)%stencil_number_2nd))
      !s_array(:) = local_der_2nd(1:mesh%plg(iv)%stencil_number_2nd)-local_der_2nd_center

      !f_array = dot_product(mesh%edt(ie)%v_weight_prime_cell(flag+1,1:mesh%plg(iv)%stencil_number_2nd), s_array(:))
      !der     = 2._r8*f_array

      !deallocate(s_array)

      allocate(s_array(1:mesh%plg_stencil_number_2nd(iv)))
      s_array(:) = local_der_2nd(1:mesh%plg_stencil_number_2nd(iv))-local_der_2nd_center

      f_array = dot_product(mesh%edt_v_weight_prime_cell(1:mesh%plg_stencil_number_2nd(iv), flag+1, ie), s_array(:))
      der     = 2._r8*f_array

      deallocate(s_array)
      return
   end subroutine calc_double_der_2nd_at_hx

   subroutine calc_der_2nd_at_hx(mesh,scalar_height_at_prime_cell,iv,ie,flag,der)
! io
   type(global_domain),  intent(in) :: mesh
   real(r8),allocatable, intent(in) :: scalar_height_at_prime_cell(:)
   integer(i4)         , intent(in) :: iv
   integer(i4)         , intent(in) :: ie
   integer(i4)         , intent(in) :: flag  ! 0 or 1
   real(r8)            , intent(out):: der
! local
   real(r8)            , allocatable :: s_array(:)
   real(r8)                          :: f_array

      !allocate(s_array(1:mesh%plg(iv)%stencil_number_2nd))

      !s_array(:) = scalar_height_at_prime_cell(mesh%plg(iv)%stencil_index_2nd(:))
      !s_array(:) = s_array(:)-scalar_height_at_prime_cell(iv)
      !f_array = dot_product(mesh%edt(ie)%v_weight_prime_cell(flag+1,1:mesh%plg(iv)%stencil_number_2nd), s_array(:))
      !der     = 2._r8*f_array

      !deallocate(s_array)
      allocate(s_array(1:mesh%plg_stencil_number_2nd(iv)))

      s_array(:) = scalar_height_at_prime_cell(mesh%plg_stencil_index_2nd(1:mesh%plg_stencil_number_2nd(iv), iv))
      s_array(:) = s_array(:)-scalar_height_at_prime_cell(iv)
      f_array = dot_product(mesh%edt_v_weight_prime_cell(1:mesh%plg_stencil_number_2nd(iv), flag+1, ie), s_array(:))
      der     = 2._r8*f_array

      deallocate(s_array)
      return
   end subroutine calc_der_2nd_at_hx

   subroutine calc_der_4th_at_hx(mesh,scalar_height_at_prime_cell,iv,ie,flag,der4,der2)
! io 
   type(global_domain),   intent(in) :: mesh
   real(r8), allocatable, intent(in) :: scalar_height_at_prime_cell(:)
   integer(i4)         ,  intent(in) :: iv
   integer(i4)         ,  intent(in) :: ie
   integer(i4)         ,  intent(in) :: flag  ! 0 or 1
   real(r8)            ,  intent(out):: der4
   real(r8)            ,  intent(out):: der2
! local
   real(r8)   , allocatable          :: s_array(:)
   real(r8)                          :: f_array(2)

      !allocate(s_array(1:mesh%plg(iv)%stencil_number_4th))

      !s_array(:) = scalar_height_at_prime_cell(mesh%plg(iv)%stencil_index_4th(:))
      !s_array(:) = s_array(:)-scalar_height_at_prime_cell(iv)

      !f_array(1) = dot_product(mesh%edt(ie)%v_weight_4th_prime_cell(flag+1,1:mesh%plg(iv)%stencil_number_4th), s_array(:))
      !f_array(2) = dot_product(mesh%edt(ie)%v_weight_4th_2nd_prime_cell(flag+1,1:mesh%plg(iv)%stencil_number_4th), s_array(:))

      !der4       = 24._r8*f_array(1)
      !der2       = 2._r8*f_array(2)

      !deallocate(s_array)

      allocate(s_array(1:mesh%plg_stencil_number_4th(iv)))

      s_array(:) = scalar_height_at_prime_cell(mesh%plg_stencil_index_4th(1:mesh%plg_stencil_number_4th(iv), iv))
      s_array(:) = s_array(:)-scalar_height_at_prime_cell(iv)

      f_array(1) = dot_product(mesh%edt_v_weight_4th_prime_cell(1:mesh%plg_stencil_number_4th(iv),flag+1,ie), s_array(:))
      f_array(2) = dot_product(mesh%edt_v_weight_4th_2nd_prime_cell(1:mesh%plg_stencil_number_4th(iv),flag+1,ie), s_array(:))

      der4       = 24._r8*f_array(1)
      der2       = 2._r8*f_array(2)
     return
   end subroutine calc_der_4th_at_hx

   subroutine calc_der_2nd_at_tr(mesh,scalar_at_dual_cell,it,ie,flag,der)
! io
   type(global_domain),   intent(in) :: mesh
   real(r8), allocatable, intent(in) :: scalar_at_dual_cell(:)
   integer(i4),           intent(in) :: it
   integer(i4),           intent(in) :: ie
   integer(i4),           intent(in) :: flag
   real(r8),              intent(out):: der
! local
   real(r8)   , allocatable          :: s_array(:)
   real(r8)                          :: f_array
  
      !allocate(s_array(1:mesh%tri(it)%stencil_number_2nd))
      !s_array(:) = scalar_at_dual_cell(mesh%tri(it)%stencil_index_2nd(:))
      !s_array(:) = s_array-scalar_at_dual_cell(it)

      !f_array    = dot_product(mesh%edp(ie)%v_weight_dual_cell(flag+1,1:mesh%tri(it)%stencil_number_2nd), s_array(:))
      !der        = 2._r8*f_array

      !deallocate(s_array)

      allocate(s_array(1:mesh%tri_stencil_number_2nd(it)))
      s_array(:) = scalar_at_dual_cell(mesh%tri_stencil_index_2nd(1:mesh%tri_stencil_number_2nd(it), it))
      s_array(:) = s_array-scalar_at_dual_cell(it)

      f_array    = dot_product(mesh%edp_v_weight_dual_cell(1:mesh%tri_stencil_number_2nd(it),flag+1,ie), s_array(:))
      der        = 2._r8*f_array

      deallocate(s_array)
     return
   end subroutine calc_der_2nd_at_tr

   subroutine calc_der_4th_at_tr(mesh,scalar_at_dual_cell,it,ie,flag,der)
! io
   type(global_domain),   intent(in) :: mesh
   real(r8), allocatable, intent(in) :: scalar_at_dual_cell(:)
   integer(i4)         ,  intent(in) :: it
   integer(i4)         ,  intent(in) :: ie
   integer(i4)         ,  intent(in) :: flag  ! 0 or 1
   real(r8)            ,  intent(out):: der
! local
   real(r8)   , allocatable          :: s_array(:)
   real(r8)                          :: f_array

      allocate(s_array(1:mesh%tri(it)%stencil_number_4th))

      s_array(:) = scalar_at_dual_cell(mesh%tri(it)%stencil_index_4th(:))
      s_array(:) = s_array(:)-scalar_at_dual_cell(it)

      f_array    = dot_product(mesh%edp(ie)%v_weight_4th_dual_cell(flag+1,1:mesh%tri(it)%stencil_number_4th), s_array(:))
      der = 24._r8*f_array

      deallocate(s_array)
     return
   end subroutine calc_der_4th_at_tr

!================================================
!           reconstruct at mid point
!================================================

   subroutine calc_mid_2nd_at_hx(mesh,scalar_at_prime_cell,iv,ie,flag,mid)
! io
    type(global_domain),   intent(in) :: mesh
    real(r8), allocatable, intent(in) :: scalar_at_prime_cell(:)
    integer(i4)         ,  intent(in) :: iv
    integer(i4)         ,  intent(in) :: ie
    integer(i4)         ,  intent(in) :: flag  ! 0 or 1
    real(r8)            ,  intent(out):: mid
! local
    real(r8), allocatable             :: b_array(:,:)
    real(r8), allocatable             :: s_array(:)
    integer(i4), allocatable          :: n_index(:)  ! index in the stencil
    real(r8)                          :: f_array(5)
    real(r8)                          :: local_nx 
    real(r8)                          :: local_ny 

      allocate(b_array(5,mesh%vtx(iv)%nnb))
      allocate(s_array(1:mesh%vtx(iv)%nnb))
      allocate(n_index(1:mesh%vtx(iv)%nnb))

!================================================
!      Compute F Matrix and 2nd-derivatives
!================================================

      local_nx   = mesh%edp(ie)%v_mid_point_lob_prime_cell(flag+1,1)
      local_ny   = mesh%edp(ie)%v_mid_point_lob_prime_cell(flag+1,2)

      b_array      = mesh%plg(iv)%blocal(1,:,:)
      n_index(:)   = mesh%vtx(iv)%nb(:)
      s_array(:)   = scalar_at_prime_cell(n_index(:))
      s_array(:)   = s_array(:)-scalar_at_prime_cell(iv)

      f_array      = matmul(b_array,s_array)

      mid   = scalar_at_prime_cell(iv)+&
                     f_array(1)*local_nx+&
                     f_array(2)*local_ny+&
                     f_array(3)*(local_nx**2)+&
                     f_array(4)*local_nx*local_ny+&
                     f_array(5)*(local_ny**2)

      return
   end subroutine calc_mid_2nd_at_hx

   subroutine calc_mid_2nd_at_tr(mesh,scalar_at_dual_cell,it,ie,flag,mid)
! io
    type(global_domain),   intent(in) :: mesh
    real(r8), allocatable, intent(in) :: scalar_at_dual_cell(:)
    integer(i4)         ,  intent(in) :: it
    integer(i4)         ,  intent(in) :: ie
    integer(i4)         ,  intent(in) :: flag  ! 0 or 1
    real(r8)            ,  intent(out):: mid
! local
    real(r8), allocatable             :: b_array(:,:)
    real(r8), allocatable             :: s_array(:)
    real(r8)                          :: f_array(5)
    real(r8)                          :: local_nx 
    real(r8)                          :: local_ny


       allocate(b_array(5,mesh%tri(it)%stencil_number_2nd))
       allocate(s_array(1:mesh%tri(it)%stencil_number_2nd))

!================================================
!      Compute F Matrix and 2nd-derivatives
!================================================

       local_nx   = mesh%edt(ie)%v_mid_point_lob_dual_cell(flag+1,1)
       local_ny   = mesh%edt(ie)%v_mid_point_lob_dual_cell(flag+1,2)

       b_array    = mesh%tri(it)%blocal(1,:,:)
       s_array(:) = scalar_at_dual_cell(mesh%tri(it)%stencil_index_2nd(:))
       s_array(:) = s_array-scalar_at_dual_cell(it)

       f_array    = matmul(b_array,s_array)

       mid        = scalar_at_dual_cell(it)+&
                           f_array(1)*local_nx+&
                           f_array(2)*local_ny+&
                           f_array(3)*(local_nx**2)+&
                           f_array(4)*(local_nx*local_ny)+&
                           f_array(5)*(local_ny**2)

       deallocate(b_array)
       deallocate(s_array)

      return

   end subroutine calc_mid_2nd_at_tr


  subroutine calc_coriolis_term(mesh,scalar_pv_at_edge,&
                                     scalar_normal_flux_at_edge,&
                                     scalar_coriolis_term_at_edge)

  use grist_nml_module, only: conserve_scheme
! io 
   type(global_domain),     intent(in)    :: mesh
   real(r8), allocatable  , intent(in)    :: scalar_pv_at_edge(:)
   real(r8), allocatable  , intent(in)    :: scalar_normal_flux_at_edge(:)
   real(r8), allocatable  , intent(inout) :: scalar_coriolis_term_at_edge(:)
! local
   integer(i4)                         :: ie
   integer(i4)                         :: kk
   real(r8)                            :: length_of_triangle
   real(r8)                            :: pv_at_edge(16)
   real(r8)                            :: cell_sum(16)

!
! te conserving scheme
!
   !if(trim(conserve_scheme).eq.'te')then
   !  do ie = 1, mesh%ne  ! for each edge e
   !     length_of_triangle = rearth*mesh%edt(ie)%leng
   !     pv_at_edge(1:mesh%edp(ie)%nedge) = scalar_pv_at_edge(ie)
   !     cell_sum(1:mesh%edp(ie)%nedge) = mesh%edp(ie)%trsk_on_edge(1:mesh%edp(ie)%nedge)*&
   !                  mesh%edp(ie)%edpl_on_edge(1:mesh%edp(ie)%nedge)*&
   !                 scalar_normal_flux_at_edge(mesh%edp(ie)%edge_on_edge(1:mesh%edp(ie)%nedge))*&
   !                 (pv_at_edge(1:mesh%edp(ie)%nedge)+scalar_pv_at_edge(mesh%edp(ie)%edge_on_edge(1:mesh%edp(ie)%nedge)))/2._r8
   !     scalar_coriolis_term_at_edge(ie) = sum(cell_sum(1:mesh%edp(ie)%nedge))/length_of_triangle
   !  end do
   !end if
   if(trim(conserve_scheme).eq.'te')then
     do ie = 1, mesh%ne  ! for each edge e
        length_of_triangle = rearth*mesh%edt_leng(ie)
        pv_at_edge(1:mesh%edp_nedge(ie)) = scalar_pv_at_edge(ie)
        cell_sum(1:mesh%edp_nedge(ie)) = mesh%edp_trsk_on_edge(1:mesh%edp_nedge(ie), ie)*&
                     mesh%edp_edpl_on_edge(1:mesh%edp_nedge(ie), ie)*&
                     scalar_normal_flux_at_edge(mesh%edp_edge_on_edge(1:mesh%edp_nedge(ie), ie))*&
                    (pv_at_edge(1:mesh%edp_nedge(ie))+scalar_pv_at_edge(mesh%edp_edge_on_edge(1:mesh%edp_nedge(ie), ie)))/2._r8
        scalar_coriolis_term_at_edge(ie) = sum(cell_sum(1:mesh%edp_nedge(ie)))/length_of_triangle
     end do
   end if
!
! pe conserving scheme
!
   !if(trim(conserve_scheme).eq.'pe')then
   !  do ie = 1, mesh%ne  ! for each edge e
   !     length_of_triangle = rearth*mesh%edt(ie)%leng
   !     pv_at_edge(1:mesh%edp(ie)%nedge) = scalar_pv_at_edge(ie)
   !     cell_sum(1:mesh%edp(ie)%nedge) = mesh%edp(ie)%trsk_on_edge(1:mesh%edp(ie)%nedge)*&
   !                mesh%edp(ie)%edpl_on_edge(1:mesh%edp(ie)%nedge)*&
   !                scalar_normal_flux_at_edge(mesh%edp(ie)%edge_on_edge(1:mesh%edp(ie)%nedge))*&
   !                pv_at_edge(1:mesh%edp(ie)%nedge)
   !     scalar_coriolis_term_at_edge(ie) = sum(cell_sum(1:mesh%edp(ie)%nedge))/length_of_triangle
   !  end do
   !end if

   if(trim(conserve_scheme).eq.'pe')then
     do ie = 1, mesh%ne  ! for each edge e
        length_of_triangle = rearth*mesh%edt_leng(ie)
        pv_at_edge(1:mesh%edp_nedge(ie)) = scalar_pv_at_edge(ie)
        cell_sum(1:mesh%edp_nedge(ie)) = mesh%edp_trsk_on_edge(1:mesh%edp_nedge(ie), ie)*&
                   mesh%edp_edpl_on_edge(1:mesh%edp_nedge(ie), ie)*&
                   scalar_normal_flux_at_edge(mesh%edp_edge_on_edge(1:mesh%edp_nedge(ie), ie))*&
                   pv_at_edge(1:mesh%edp_nedge(ie))
        scalar_coriolis_term_at_edge(ie) = sum(cell_sum(1:mesh%edp_nedge(ie)))/length_of_triangle
     end do
   end if
   return
  end subroutine calc_coriolis_term

!==================================================
!  Given a normal component of vector, reconstruct
! U and V components
!==================================================

  subroutine project_uv(mesh,scalar_normal_velocity_at_edge,scalar_u,scalar_v)

! io
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge(:)
   real(r8), allocatable, intent(inout) :: scalar_u(:)
   real(r8), allocatable, intent(inout) :: scalar_v(:)
! local
   real(r8), allocatable                :: scalar_tangent_velocity_at_edge(:)
   real(r8), allocatable                :: tmp(:)
   real(r8)                             :: cos_theta_east
   real(r8)                             :: cos_theta_north
   real(r8)                             :: vector(1:3)
   real(r8)                             :: utmp1
   real(r8)                             :: vtmp1
   integer(i4)                          :: ie

   allocate(scalar_tangent_velocity_at_edge(mesh%ne_full))
   allocate(tmp(mesh%ne_full))

   tmp(:) = 1._r8
!
!    We let pv_at_edge(tmp)=1, and input Normal 
!  V, routine calc_coriolis_term can be directly
!  used to reconstruct tangent wind.
!

   call calc_coriolis_term(mesh,tmp,&
                                scalar_normal_velocity_at_edge,&
                                scalar_tangent_velocity_at_edge)

   tmp = scalar_normal_velocity_at_edge

   do ie  = 1, mesh%ne

       vector(:) = scalar_normal_velocity_at_edge(ie)*mesh%edp(ie)%nr+&
                  scalar_tangent_velocity_at_edge(ie)*mesh%edp(ie)%tg

       call convert_vector_cart2sph(mesh%edp(ie)%c%p, vector , utmp1, vtmp1)

       scalar_u(ie) = utmp1
       scalar_v(ie) = vtmp1
!================================================
!   Alternative way to get U&V direction
!   Not used and not completed, lack tangent. 
!================================================
!  cos_theta_east = dot_product(mesh%edp(ie)%c%east,mesh%edp(ie)%nr)
!  cos_theta_north= dot_product(mesh%edp(ie)%c%north,mesh%edp(ie)%nr)
!  scalar_u(ie) = scalar_normal_velocity_at_edge(ie)*cos_theta_east
!  scalar_v(ie) = scalar_normal_velocity_at_edge(ie)*cos_theta_north
!================================================
    end do

    deallocate(scalar_tangent_velocity_at_edge)
    deallocate(tmp)

    return
   end subroutine project_uv

!
! use perot reconstrucion to recover UV at centroid
! based on normal wind at edges.
!
  subroutine vector_recon_perot_edge2cell_uv(mesh,scalar_normal_velocity_at_edge,scalar_u,scalar_v)

! io
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge(:)
   real(r8), allocatable, intent(inout) :: scalar_u(:)  ! at Voronoi centroid
   real(r8), allocatable, intent(inout) :: scalar_v(:)  ! at Voronoi centroid
! local
   real(r8)                             :: vector_sum(1:3),vector_tmp(1:3),local_basis(1:3)
   real(r8)                             :: utmp, vtmp
   integer(i4)                          :: iv ,inb, ie


      do iv  = 1, mesh%nv
         vector_sum = 0._r8
         do inb  = 1, mesh%vtx(iv)%nnb
            ie = mesh%vtx(iv)%ed(inb)
            call project_sphere2tplane(mesh%edt(ie)%c%p,mesh%vtx(iv)%p,local_basis) ! unit sphere
            vector_tmp = dot_product(local_basis/norm(local_basis),mesh%edp(ie)%nr)*local_basis*mesh%edp(ie)%leng/(mesh%plg(iv)%areag)
            vector_sum = vector_sum + vector_tmp*scalar_normal_velocity_at_edge(ie)
         end do
         call convert_vector_cart2sph(mesh%vtx(iv)%p,vector_sum,utmp,vtmp)
         scalar_u(iv)  = utmp
         scalar_v(iv)  = vtmp
      end do

      return
  
  end subroutine vector_recon_perot_edge2cell_uv

!
! use perot reconstrucion to recover UV at tr's center
! based on normal wind at hexaonal edges/tangent of tr edges
!
  subroutine vector_recon_perot_edge2trcell_uv(mesh,scalar_normal_velocity_at_edge,scalar_u,scalar_v)
! io
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge(:)
   real(r8), allocatable, intent(inout) :: scalar_u(:)  ! at Voronoi centroid
   real(r8), allocatable, intent(inout) :: scalar_v(:)  ! at Voronoi centroid
! local
   real(r8)                             :: vector_sum(1:3),vector_tmp(1:3),tmp,local_basis(1:3)
   real(r8)                             :: utmp, vtmp
   integer(i4)                          :: it ,inb, ie


      do it  = 1, mesh%nt
         vector_sum = 0._r8
         do inb  = 1, mesh%tri_nnb(it) 
            ie = mesh%tri(it)%ed(inb)
            call project_sphere2tplane(mesh%edt(ie)%c%p,mesh%tri(it)%c%p,local_basis) ! unit sphere
            tmp        = dot_product(mesh%edp(ie)%nr,cross_product(mesh%edt(ie)%c%p,local_basis)/norm(cross_product(mesh%edt(ie)%c%p,local_basis)))
            vector_tmp = cross_product(mesh%edt(ie)%c%p,local_basis)*tmp*mesh%edt(ie)%leng/(mesh%tri(it)%areag)
            vector_sum = vector_sum + vector_tmp*scalar_normal_velocity_at_edge(ie)
         end do
         call convert_vector_cart2sph(mesh%tri(it)%c%p,vector_sum,utmp,vtmp)
         scalar_u(it)  = utmp
         scalar_v(it)  = vtmp
      end do

      return
   end subroutine vector_recon_perot_edge2trcell_uv

!================================================
!   Remap scalar from prime cell to dual cell
!================================================

  subroutine calc_scalar_at_dual_cell(mesh,scalar_height_at_prime_cell,&
                                           scalar_height_at_dual_cell)
! io
   type(global_domain),   intent(in)   :: mesh
   real(r8), allocatable, intent(in)   :: scalar_height_at_prime_cell(:)
   real(r8), allocatable, intent(inout):: scalar_height_at_dual_cell(:)
! local
   integer(i4)                         :: it
   integer(i4)                         :: ie
   integer(i4)                         :: index_of_cell
   real(r8)                            :: area_of_cell
   real(r8)                            :: height_at_cell
   real(r8)                            :: Riv_of_cell
   real(r8)                            :: tmp


   do it = 1, mesh%nt
     tmp = 0._r8
!================================================
! index of prime cells that form this triangle
!================================================
     do ie = 1, mesh%tri_nnb(it)
        index_of_cell   = mesh%tri(it)%v(ie)
        height_at_cell  = scalar_height_at_prime_cell(index_of_cell)
        tmp             = tmp+(mesh%tri(it)%kite_area(ie))*height_at_cell
     end do
        scalar_height_at_dual_cell(it) = tmp/(mesh%tri(it)%areag)
   end do

   return
  end subroutine calc_scalar_at_dual_cell
    
!================================================
! use same B-local for each edge
!================================================

  subroutine calc_mr07_area_at_hx(mesh,scalar_at_prime_cell,iv,lob_x,lob_y,scalar_upwind)
! io
   type(global_domain),   intent(in) :: mesh
   real(r8), allocatable, intent(in) :: scalar_at_prime_cell(:)
   integer(i4)         ,  intent(in) :: iv
   real(r8)            ,  intent(in) :: lob_x, lob_y
   real(r8)            ,  intent(out):: scalar_upwind
! local
   real(r8)            ,  allocatable:: b_array(:,:)
   real(r8)            ,  allocatable:: s_array(:)
   real(r8)                          :: f_array(5)


      allocate(b_array(5,mesh%plg(iv)%stencil_number_2nd))
      allocate(s_array(1:mesh%plg(iv)%stencil_number_2nd))

      b_array       = mesh%plg(iv)%blocal(1,:,:)
      s_array(:)    = scalar_at_prime_cell(mesh%plg(iv)%stencil_index_2nd(:))
      s_array(:)    = s_array(:)-scalar_at_prime_cell(iv)
      f_array       = matmul(b_array,s_array)
      scalar_upwind = polyfit2d(scalar_at_prime_cell(iv),f_array,lob_x,lob_y)

      deallocate(b_array)
      deallocate(s_array)

      return
  end subroutine calc_mr07_area_at_hx

  subroutine calc_sm10_area_at_hx(mesh, scalar_at_prime_cell, iv, & 
                                        g1_lob_x, g1_lob_y,&
                                        g2_lob_x, g2_lob_y,&
                                        g3_lob_x, g3_lob_y,&
                                        g4_lob_x, g4_lob_y,&
                                        g5_lob_x, g5_lob_y,&
                                        scalar_upwind)
! io
   type(global_domain),   intent(in) :: mesh
   real(r8), allocatable, intent(in) :: scalar_at_prime_cell(:)
   integer(i4)         ,  intent(in) :: iv
   real(r8)            ,  intent(in) :: g1_lob_x, g1_lob_y
   real(r8)            ,  intent(in) :: g2_lob_x, g2_lob_y
   real(r8)            ,  intent(in) :: g3_lob_x, g3_lob_y
   real(r8)            ,  intent(in) :: g4_lob_x, g4_lob_y
   real(r8)            ,  intent(in) :: g5_lob_x, g5_lob_y
   real(r8)            ,  intent(out):: scalar_upwind
! local
   real(r8)            , allocatable :: b_array(:,:)
   real(r8)            , allocatable :: s_array(:)
   real(r8)                          :: f_array(5)
   real(r8)                          :: sub_gauss_inte
   real(r8)                          :: total_sum
   real(r8)                          :: correction
   integer(i4)                       :: ii,inb

      allocate(b_array(5,mesh%plg(iv)%stencil_number_2nd))
      allocate(s_array(1:mesh%plg(iv)%stencil_number_2nd))

      b_array    = mesh%plg(iv)%blocal(1,:,:)
      s_array(:) = scalar_at_prime_cell(mesh%plg(iv)%stencil_index_2nd(:))
      s_array(:) = s_array(:)-scalar_at_prime_cell(iv)
      f_array    = matmul(b_array,s_array)
!
! compute correction
!
      total_sum = 0._r8

      do inb = 1, mesh%vtx(iv)%nnb
         sub_gauss_inte    = 0._r8
         do ii = 1, 3
            sub_gauss_inte = sub_gauss_inte+polyfit2d(scalar_at_prime_cell(iv), &
                                                      f_array,mesh%plg(iv)%sub_triangle_midp_lob(inb,ii,1),&
                                                              mesh%plg(iv)%sub_triangle_midp_lob(inb,ii,2))
         end do
         sub_gauss_inte = sub_gauss_inte*mesh%plg(iv)%sub_triangle_area(inb)/3._r8
         total_sum      = total_sum+sub_gauss_inte
      end do
      total_sum     = total_sum/mesh%plg(iv)%areag
      correction    = total_sum-scalar_at_prime_cell(iv)

      scalar_upwind = (2*polyfit2d(scalar_at_prime_cell(iv), f_array,g1_lob_x,g1_lob_y)+&
                         polyfit2d(scalar_at_prime_cell(iv), f_array,g2_lob_x,g2_lob_y)+& 
                         polyfit2d(scalar_at_prime_cell(iv), f_array,g3_lob_x,g3_lob_y)+& 
                         polyfit2d(scalar_at_prime_cell(iv), f_array,g4_lob_x,g4_lob_y)+& 
                         polyfit2d(scalar_at_prime_cell(iv), f_array,g5_lob_x,g5_lob_y))/6._r8-correction

      return
  end subroutine calc_sm10_area_at_hx

  subroutine calc_mr07_area_at_tr(mesh,scalar_at_prime_cell,it,lob_x,lob_y,scalar_upwind)
! io
   type(global_domain),   intent(in) :: mesh
   real(r8), allocatable, intent(in) :: scalar_at_prime_cell(:)
   integer(i4)         ,  intent(in) :: it
   real(r8)            ,  intent(in) :: lob_x, lob_y
   real(r8)            ,  intent(out):: scalar_upwind
! local
   real(r8)            ,  allocatable:: b_array(:,:)
   real(r8)            ,  allocatable:: s_array(:)
   real(r8)                          :: f_array(5)

      allocate(b_array(5,mesh%tri(it)%stencil_number_2nd))
      allocate(s_array(1:mesh%tri(it)%stencil_number_2nd))

      b_array       = mesh%tri(it)%blocal(1,:,:)
      s_array(:)    = scalar_at_prime_cell(mesh%tri(it)%stencil_index_2nd(:))
      s_array(:)    = s_array(:)-scalar_at_prime_cell(it)
      f_array       = matmul(b_array,s_array)
      scalar_upwind = polyfit2d(scalar_at_prime_cell(it),f_array,lob_x,lob_y)

      deallocate(b_array)
      deallocate(s_array)

      return
  end subroutine calc_mr07_area_at_tr
  
  subroutine calc_sm10_area_at_tr(mesh, scalar_at_dual_cell, it, & 
                                        g1_lob_x, g1_lob_y,&
                                        g2_lob_x, g2_lob_y,&
                                        g3_lob_x, g3_lob_y,&
                                        g4_lob_x, g4_lob_y,&
                                        g5_lob_x, g5_lob_y,&
                                        scalar_upwind)
! io
   type(global_domain),   intent(in) :: mesh
   real(r8), allocatable, intent(in) :: scalar_at_dual_cell(:)
   integer(i4)         ,  intent(in) :: it
   real(r8)            ,  intent(in) :: g1_lob_x, g1_lob_y
   real(r8)            ,  intent(in) :: g2_lob_x, g2_lob_y
   real(r8)            ,  intent(in) :: g3_lob_x, g3_lob_y
   real(r8)            ,  intent(in) :: g4_lob_x, g4_lob_y
   real(r8)            ,  intent(in) :: g5_lob_x, g5_lob_y
   real(r8)            ,  intent(out):: scalar_upwind
! local
   real(r8)            , allocatable :: b_array(:,:)
   real(r8)            , allocatable :: s_array(:)
   real(r8)                          :: f_array(5)
   real(r8)                          :: sub_gauss_inte
   real(r8)                          :: gauss(3), gauss_lob_x, gauss_lob_y
   real(r8)                          :: total_sum
   real(r8)                          :: correction
   integer(i4)                       :: ii,inb

      allocate(b_array(5,mesh%tri(it)%stencil_number_2nd))
      allocate(s_array(1:mesh%tri(it)%stencil_number_2nd))

      b_array    = mesh%tri(it)%blocal(1,:,:)
      s_array(:) = scalar_at_dual_cell(mesh%tri(it)%stencil_index_2nd(:))
      s_array(:) = s_array(:)-scalar_at_dual_cell(it)
      f_array    = matmul(b_array,s_array)
!
! compute correction
!
       sub_gauss_inte    = 0._r8
       do ii = 1, 3
          gauss(:)    = mesh%edt(mesh%tri(it)%ed(ii))%c%p-mesh%tri(it)%c%p
          gauss_lob_x = dot_product(gauss,mesh%tri(it)%lob_nx(1,:))
          gauss_lob_y = dot_product(gauss,mesh%tri(it)%lob_ny(1,:))
          sub_gauss_inte = sub_gauss_inte+polyfit2d(scalar_at_dual_cell(it), &
                                                    f_array,gauss_lob_x,gauss_lob_y)
       end do
       sub_gauss_inte = sub_gauss_inte/3._r8
       correction     = sub_gauss_inte-scalar_at_dual_cell(it)

       scalar_upwind  = (2*polyfit2d(scalar_at_dual_cell(it), f_array,g1_lob_x,g1_lob_y)+&
                           polyfit2d(scalar_at_dual_cell(it), f_array,g2_lob_x,g2_lob_y)+& 
                           polyfit2d(scalar_at_dual_cell(it), f_array,g3_lob_x,g3_lob_y)+& 
                           polyfit2d(scalar_at_dual_cell(it), f_array,g4_lob_x,g4_lob_y)+& 
                           polyfit2d(scalar_at_dual_cell(it), f_array,g5_lob_x,g5_lob_y))/6._r8-correction

       return
  end subroutine calc_sm10_area_at_tr
!================================================
! PRIVATE
!================================================
   real(r8) function polyfit2d(scalar0,f_array,lob_x,lob_y)
! io
   real(r8),  intent(in) :: scalar0
   real(r8),  intent(in) :: f_array(5)
   real(r8),  intent(in) :: lob_x, lob_y

     polyfit2d = scalar0+f_array(1)*lob_x+&
                         f_array(2)*lob_y+&
                         f_array(3)*lob_x**2+&
                         f_array(4)*lob_x*lob_y+&
                         f_array(5)*lob_y**2
    return
   end function polyfit2d
  end module grist_recon_module
