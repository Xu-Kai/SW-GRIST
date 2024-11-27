
!----------------------------------------------------------------------------
! Created on 2016
! Version 1.0
! Description: Gradient, Curl, Divergence operators
! Revision history:
!----------------------------------------------------------------------------

 module grist_gcd_module

  use grist_constants,       only: gravity, i4, r8, rearth, omega
  use grist_domain_types  ,  only: global_domain
  use grist_recon_module,    only: calc_der_1st_at_hx

  implicit none

  private

  public :: gradient_operator     , &
            gradient_operator_quad, &
            curl_operator         , &
            divergence_operator   , &
            divergence_operator_tr, &
            divergence_operator_tr_smth, &
            calc_vorticity_at_dual_cell, &
            calc_vorticity_at_dual_cell_gass

  contains

  subroutine gradient_operator(mesh,scalar_at_prime_cell,gradient_at_prime_edge)
!io
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_at_prime_cell(:)
   real(r8), allocatable, intent(inout) :: gradient_at_prime_edge(:)
!
!local
!
   integer(i4)                          :: ie, v0, v1
   real(r8)                             :: phi0, phi1
   real(r8)                             :: flag
   real(r8)                             :: v0v1(3)

! gradient is defined at prime's edge point

   do ie = 1, mesh%ne ! global index
      !v0                         = mesh%edt(ie)%v(1)
      !v1                         = mesh%edt(ie)%v(2)
      !phi0                       = scalar_at_prime_cell(v0)
      !phi1                       = scalar_at_prime_cell(v1)
      !!v0v1                       = mesh%vtx(v1)%p-mesh%vtx(v0)%p
      !!flag                       = sign(1._r8,dot_product(v0v1,mesh%edp(ie)%nr))
      !!gradient_at_prime_edge(ie) = flag*(phi1-phi0)/(rearth*mesh%edt(ie)%leng)
      !gradient_at_prime_edge(ie) = mesh%edt(ie)%edpNr_edtTg*(phi1-phi0)/(rearth*mesh%edt(ie)%leng)
      v0                         = mesh%edt_v(1, ie)
      v1                         = mesh%edt_v(2, ie)
      phi0                       = scalar_at_prime_cell(v0)
      phi1                       = scalar_at_prime_cell(v1)
      !v0v1                       = mesh%vtx(v1)%p-mesh%vtx(v0)%p
      !flag                       = sign(1._r8,dot_product(v0v1,mesh%edp(ie)%nr))
      !gradient_at_prime_edge(ie) = flag*(phi1-phi0)/(rearth*mesh%edt(ie)%leng)
      gradient_at_prime_edge(ie) = mesh%edt_edpNr_edtTg(ie)*(phi1-phi0)/(rearth*mesh%edt_leng(ie))
   end do

   return
  end subroutine gradient_operator

  subroutine gradient_operator_quad(mesh,scalar_at_prime_cell,gradient_at_prime_edge)
!io
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_at_prime_cell(:)
   real(r8), allocatable, intent(inout) :: gradient_at_prime_edge(:)
!
!local
!
   integer(i4)                         :: ie, mark
   integer(i4)                         :: v0, v1, v_upwind
   real(r8)                            :: phi0, phi1
   real(r8)                            :: v0v1(3)
   real(r8)                            :: flag
   real(r8)                            :: der

! gradient is defined at prime's edge point

    do ie = 1, mesh%ne
       v0        = mesh%edt(ie)%v(1)       ! 1st end of edge
       v1        = mesh%edt(ie)%v(2)       ! 2nd end of edge
       v0v1      = mesh%vtx(v1)%p-mesh%vtx(v0)%p
       flag      = dot_product(v0v1,mesh%edp(ie)%nr)
       mark      = int(-0.5*flag+1.5)
       v_upwind  = mesh%edt(ie)%v(mark)
       call calc_der_1st_at_hx(mesh,scalar_at_prime_cell,v_upwind,ie,mark-1,der)
       gradient_at_prime_edge(ie) = der/rearth
    end do

   return
  end subroutine gradient_operator_quad

  subroutine curl_operator(mesh, scalar_normal_velocity_at_edge, &
                                 scalar_relative_vorticity_at_dual_cell)
! io
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge(:)
   real(r8), allocatable, intent(inout) :: scalar_relative_vorticity_at_dual_cell(:)
! local
   integer(i4)                          :: it
   integer(i4)                          :: ie
   integer(i4)                          :: index_of_edge
   real(r8)                             :: length_of_edge
   real(r8)                             :: scalar_normal_velocity
   integer(i4)                          :: t_edge_v
   real(r8)                             :: tmp
  
   !do it = 1, mesh%nt ! for each dual cell/triangle
   !   tmp = 0._r8
   !   do ie = 1, 3    ! for each edge of triangle
   !      index_of_edge  = mesh%tri(it)%ed(ie)
   !      length_of_edge = rearth*mesh%edt(index_of_edge)%leng
   !      scalar_normal_velocity = scalar_normal_velocity_at_edge(index_of_edge)
! s!ign correction
   !      t_edge_v   = mesh%edt(index_of_edge)%edpTg_edtNr*mesh%tri(it)%nr(ie)
   !      tmp  = tmp+t_edge_v*length_of_edge*scalar_normal_velocity
   !   end do
   !   scalar_relative_vorticity_at_dual_cell(it) = tmp/((rearth**2)*mesh%tri(it)%areag)
   ! end do
   do it = 1, mesh%nt ! for each dual cell/triangle
      tmp = 0._r8
      do ie = 1, mesh%tri_nnb(it)    ! for each edge of triangle
         index_of_edge  = mesh%tri_ed(ie, it)
         length_of_edge = rearth*mesh%edt_leng(index_of_edge)
         scalar_normal_velocity = scalar_normal_velocity_at_edge(index_of_edge)
! sign correction
         t_edge_v   = mesh%edt_edpTg_edtNr(index_of_edge)*mesh%tri_nr(ie, it)
         tmp  = tmp+t_edge_v*length_of_edge*scalar_normal_velocity
      end do
      scalar_relative_vorticity_at_dual_cell(it) = tmp/((rearth**2)*mesh%tri_areag(it))
    end do

    return
  end subroutine curl_operator

  subroutine divergence_operator(mesh, scalar_normal_flux_at_edge, &
                                       tend_scalar)
! io
   type(global_domain)  ,  intent(in)   :: mesh
   real(r8), allocatable,  intent(in)   :: scalar_normal_flux_at_edge(:)
   real(r8), allocatable,  intent(inout):: tend_scalar(:)
! local
   real(r8)                             :: div_sum
   integer(i4)                          :: iv
   integer(i4)                          :: ineigh
   integer(i4)                          :: index_edge

    do iv = 1, mesh%nv     ! for each vertice in the mesh
       !div_sum = 0._r8
       !do ineigh = 1, mesh%vtx(iv)%nnb
       !   index_edge  = mesh%vtx(iv)%ed(ineigh)
       !   div_sum     = div_sum+scalar_normal_flux_at_edge(index_edge)*&
       !                         mesh%plg(iv)%nr(ineigh)*&
       !                         rearth*mesh%edp(index_edge)%leng

       !end do
       !tend_scalar(iv) = div_sum/((rearth**2)*mesh%plg(iv)%areag) ! GRAD.(F)

       div_sum = 0._r8
       do ineigh = 1, mesh%vtx_nnb(iv)
          index_edge  = mesh%vtx_ed(ineigh, iv)
          div_sum     = div_sum+scalar_normal_flux_at_edge(index_edge)*&
                                mesh%plg_nr(ineigh, iv)*&
                                rearth*mesh%edp_leng(index_edge)

       end do

       tend_scalar(iv) = div_sum/((rearth**2)*mesh%plg_areag(iv)) ! GRAD.(F)
    end do

   return
  end subroutine divergence_operator

  subroutine divergence_operator_tr(mesh,&
                                    scalar_normal_flux_at_edge,&
                                    tend_scalar)
! io
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_normal_flux_at_edge(:)
   real(r8), allocatable, intent(inout) :: tend_scalar(:)
! local
   integer(i4)                          :: it
   integer(i4)                          :: ineigh
   integer(i4)                          :: index_edge
   real(r8)                             :: div_sum
  

    do it = 1, mesh%nt     ! for each vertice in the mesh
       div_sum = 0._r8
       do ineigh = 1, mesh%tri_nnb(it) 
          index_edge  = mesh%tri(it)%ed(ineigh)
          div_sum     = div_sum+scalar_normal_flux_at_edge(index_edge)*&
                                mesh%tri(it)%nr(ineigh)*&
                                rearth*mesh%edt(index_edge)%leng
       end do
       tend_scalar(it) = div_sum/((rearth**2)*mesh%tri(it)%areag) ! no sign
    end do
    
   return
  end subroutine divergence_operator_tr

  subroutine divergence_operator_tr_smth(mesh,&
                                    scalar_normal_flux_at_edge,&
                                    tend_scalar)
! io
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_normal_flux_at_edge(:)
   real(r8), allocatable, intent(inout) :: tend_scalar(:)
! local
   integer(i4)                          :: it
   integer(i4)                          :: ineigh
   integer(i4)                          :: index_edge
   integer(i4)                          :: index_cell
   real(r8)                             :: div_sum
  
    do it = 1, mesh%nt     ! for each vertice in the mesh
       div_sum = 0._r8
       do ineigh = 1, mesh%tri_nnb(it) 
          index_edge  = mesh%tri(it)%ed(ineigh)
          div_sum     = div_sum+scalar_normal_flux_at_edge(index_edge)*&
                                mesh%tri(it)%nr(ineigh)*&
                                rearth*mesh%edt(index_edge)%leng
       end do
       tend_scalar(it) = div_sum/((rearth**2)*mesh%tri(it)%areag) ! no sign
    end do

! smooth as ICON
    do it = 1, mesh%nt
       div_sum = 0._r8
       do ineigh = 1, mesh%tri_nnb(it)
          index_cell = mesh%tri(it)%nb(ineigh)
          div_sum    = div_sum+(tend_scalar(index_cell)*mesh%tri(index_cell)%areag+&
                                tend_scalar(it)*mesh%tri(it)%areag)/(mesh%tri(index_cell)%areag+mesh%tri(it)%areag)
       end do
          tend_scalar(it) = div_sum/3._r8
    end do
    
   return
  end subroutine divergence_operator_tr_smth

!================================================
!     Calculate A&P vorticity at dual cell
!================================================

  subroutine calc_vorticity_at_dual_cell(mesh, scalar_normal_velocity_at_edge,&
                                               scalar_absolute_vorticity_at_dual_cell,&
                                               scalar_relative_vorticity_at_dual_cell)
! io
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge(:)
   real(r8), allocatable, intent(inout) :: scalar_absolute_vorticity_at_dual_cell(:)
   real(r8), allocatable, intent(inout) :: scalar_relative_vorticity_at_dual_cell(:)
! local
   integer(i4)                         :: ie
   integer(i4)                         :: it
   integer(i4)                         :: v1
   integer(i4)                         :: v2
   real(r8)                            :: coriolis

   call curl_operator(mesh,scalar_normal_velocity_at_edge,scalar_relative_vorticity_at_dual_cell)

   do it = 1, mesh%nt
      !coriolis = 2._r8*omega*sin(mesh%tri(it)%c%lat)
      coriolis = 2._r8*omega*sin(mesh%tri_c_lat(it))
      scalar_absolute_vorticity_at_dual_cell(it) = scalar_relative_vorticity_at_dual_cell(it)+coriolis
   end do

   return
  end subroutine calc_vorticity_at_dual_cell

!
! evaluate rhombi-based PV, same input/output as dual-cell one
! except an intermediate edge-based PV is defined
!
  subroutine calc_vorticity_at_dual_cell_gass(mesh, scalar_normal_velocity_at_edge,&
                                                    scalar_absolute_vorticity_at_dual_cell,&
                                                    scalar_relative_vorticity_at_dual_cell)
! io
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge(:)
   real(r8), allocatable, intent(inout) :: scalar_absolute_vorticity_at_dual_cell(:)
   real(r8), allocatable, intent(inout) :: scalar_relative_vorticity_at_dual_cell(:)
! local
   integer(i4)                          :: ie, it, v1, v2, inb, edge_index
   real(r8)                             :: coriolis, tmp_sum
   real(r8), allocatable                :: scalar_relative_vorticity_at_edge(:)

   if(.not.allocated(scalar_relative_vorticity_at_edge)) allocate(scalar_relative_vorticity_at_edge(mesh%ne))
!
! evaluate rhombi-based relative vorticity at edge
!
   do ie  = 1, mesh%ne
      tmp_sum = 0._r8
      do inb  = 1, 4
         edge_index = mesh%edt(ie)%edge_on_rhombi(inb)
         tmp_sum = tmp_sum+scalar_normal_velocity_at_edge(edge_index)*&
                           mesh%edt(ie)%eccw_on_rhombi(inb)*&
                           mesh%edt(ie)%edtl_on_rhombi(inb)*rearth
      end do
      scalar_relative_vorticity_at_edge(ie) = tmp_sum/(rearth*rearth*mesh%edt(ie)%area_on_rhombi)
   end do
!
! average to dual cell
!
   do it = 1, mesh%nt
      tmp_sum = 0._r8
      do inb = 1, mesh%tri_nnb(it)
         tmp_sum = tmp_sum+scalar_relative_vorticity_at_edge(mesh%tri(it)%ed(inb))
      end do
      scalar_relative_vorticity_at_dual_cell(it) = tmp_sum/3._r8
!
! absolute vorticity
!
      coriolis = 2._r8*omega*sin(mesh%tri(it)%c%lat)
      scalar_absolute_vorticity_at_dual_cell(it) = scalar_relative_vorticity_at_dual_cell(it)+coriolis
   end do

   deallocate(scalar_relative_vorticity_at_edge)

   return
  end subroutine calc_vorticity_at_dual_cell_gass

 end module grist_gcd_module
