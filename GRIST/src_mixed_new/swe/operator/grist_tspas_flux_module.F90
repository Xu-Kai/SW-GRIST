
!----------------------------------------------------------------------------
! Created on 2016
! Version 1.0
! Description: An efficiently optimized implemetation of TSPAS
!              Only calculate flux using tspas, no div operator
!              same as before to within time truncation error (really fast!)
! Revision history:
!----------------------------------------------------------------------------

 module grist_tspas_flux_module

  use grist_constants,        only: i4, r8, rearth
  use grist_domain_types,     only: global_domain
  use grist_base_flux_module, only: flux_laxwen, flux_upwind_nowind, flux_laxwen_nowind

  implicit none
  private
  public :: calc_hx_normal_flux_tspas

  real(r8), parameter  :: one  = 1._r8
  real(r8), parameter  :: kesi = 1.e-80_r8

  contains

  subroutine calc_hx_normal_flux_tspas(mesh,scalar_normal_velocity_at_edge,&
                                            scalar_at_prime_cell          ,&
                                            scalar_tspas_flux_at_edge     ,&  
                                            dtime                        ) 
   implicit none
! io
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge(:)
   real(r8), allocatable, intent(in)    :: scalar_at_prime_cell(:)
   real(r8), allocatable, intent(inout) :: scalar_tspas_flux_at_edge(:)
   real(r8),              intent(in)    :: dtime
! type
   real(r8), allocatable    :: avalue(:)
   real(r8), allocatable    :: scalar_value_astk(:)
! real
   real(r8), allocatable    :: gama_compo(:)
   real(r8), allocatable    :: lw_flux(:)
   real(r8)                 :: v1v2(1:3)
   real(r8)                 :: wind_edge, c_value, cc_value
   integer(i4)              :: idc_value
   real(r8)                 :: div_sum
   real(r8)                 :: scalar_max
   real(r8)                 :: scalar_min
   real(r8)                 :: gama
   real(r8)                 :: beta
! integer
   integer(i4)              :: v1,v2 
   integer(i4)              :: flag
   integer(i4)              :: index_edge
   integer                  :: inb,kk,iv,ie

   if(.not.allocated(avalue))            allocate(avalue(ubound(scalar_at_prime_cell,1)))  ! Just a initilization
   if(.not.allocated(scalar_value_astk)) allocate(scalar_value_astk(ubound(scalar_at_prime_cell,1)))

 !
 ! Compute LW flux and store at edges
 !
   allocate(lw_flux(1:mesh%ne_halo(1)))

   do ie = 1, mesh%ne_halo(1)
        v1          = mesh%edt(ie)%v(1)       ! 1st end of this edge
        v2          = mesh%edt(ie)%v(2)       ! 2nd end of this edge
        v1v2        = mesh%vtx(v2)%p-mesh%vtx(v1)%p
        flag        = sign(one,dot_product(v1v2,mesh%edp(ie)%nr))
        wind_edge   = flag*scalar_normal_velocity_at_edge(ie)
        call flux_laxwen(wind_edge, scalar_at_prime_cell(v1),scalar_at_prime_cell(v2),&
                         dtime,rearth*(mesh%edt(ie)%leng),lw_flux(ie))
        lw_flux(ie) = flag*lw_flux(ie)
   end do

 ! First integration
 ! for each vertice
 ! get min and max around this vertice
 ! and PRECALCULATE the updated value using LW

   do iv=1, mesh%nv_halo(1)     ! for each vertice in the mesh

      allocate(gama_compo(1:mesh%vtx(iv)%nnb))
 !
 ! calculate beta
 !
      do ie = 1, mesh%vtx(iv)%nnb
           index_edge = mesh%vtx(iv)%ed(ie)
           wind_edge      = scalar_normal_velocity_at_edge(index_edge)*mesh%plg(iv)%nr(ie) ! to ensure outward positive
           gama_compo(ie) = abs(wind_edge)*(1._r8-abs(wind_edge)*dtime/(rearth*mesh%edt(index_edge)%leng))*&
                                       rearth*mesh%edp(index_edge)%leng
      end do
        gama       = gama_compo(1)
        do  inb = 1, mesh%vtx(iv)%nnb
            gama   = max(gama,gama_compo(inb))
        end do
        deallocate(gama_compo)
        kk         = 3   ! empirical parameter
        beta       = 2._r8/(2._r8-kk*dtime*gama/((rearth**2)*mesh%plg(iv)%areag))
        beta       = max(beta,1._r8)
 !
 ! pre-update scalar
 !
      div_sum = 0._r8
      do ie = 1, mesh%vtx(iv)%nnb
         div_sum = div_sum+beta*lw_flux(mesh%vtx(iv)%ed(ie))*(mesh%plg(iv)%nr(ie))*&
                           (rearth*mesh%edp(mesh%vtx(iv)%ed(ie))%leng)
      end do
      scalar_value_astk(iv) = scalar_at_prime_cell(iv)-1._r8*dtime*div_sum/((rearth**2)*mesh%plg(iv)%areag)

        scalar_max       = scalar_at_prime_cell(iv)
        scalar_min       = scalar_at_prime_cell(iv)

        do  inb = 1, mesh%vtx(iv)%nnb
            scalar_max   = max(scalar_max,scalar_at_prime_cell(mesh%vtx(iv)%nb(inb)))
            scalar_min   = min(scalar_min,scalar_at_prime_cell(mesh%vtx(iv)%nb(inb)))
        end do

        avalue(iv)     = (scalar_value_astk(iv)-scalar_max)*(scalar_value_astk(iv)-scalar_min)

   end do
 ! End of First integration

! unchanged above, check unsigned reconstructed scalar at edge (q)

   do ie = 1, mesh%ne

      v1   = mesh%edt(ie)%v(1)
      v2   = mesh%edt(ie)%v(2)
      v1v2 = mesh%vtx(v2)%p-mesh%vtx(v1)%p
      flag = sign(one,dot_product(v1v2,mesh%edp(ie)%nr))  ! actually always 1, but we explicitly code here
      wind_edge = flag*scalar_normal_velocity_at_edge(ie)
      c_value = 1._r8
      if(avalue(v1).lt.0._r8.and.avalue(v2).lt.0._r8) c_value = 0._r8
      cc_value  = c_value+(one-c_value)*abs(wind_edge)*dtime/(rearth*mesh%edt(ie)%leng)

      scalar_tspas_flux_at_edge(ie) = 0.5_r8*(scalar_at_prime_cell(v1)+scalar_at_prime_cell(v2))-&
                                      0.5_r8*sign(1._r8,wind_edge)*(scalar_at_prime_cell(v2)-scalar_at_prime_cell(v1))*cc_value

   end do

   deallocate(scalar_value_astk)
   deallocate(avalue)
   deallocate(lw_flux)
  
   return
  end subroutine calc_hx_normal_flux_tspas

  end module grist_tspas_flux_module
