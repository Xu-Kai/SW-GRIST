
!----------------------------------------------------------------------------
! Created on 2016
! Version 1.0
! Description: Flux limiter for passive advection 
!             [1] Zalesak FCT
!             [2] Thuburn (removed at some time; find in my previous log)
! Revision history:
!----------------------------------------------------------------------------

 module grist_limiter_module

  use grist_constants,      only: i4, r8, rearth,dbleps
  use grist_domain_types  , only: global_domain

  implicit none

  private

  public :: flux_limiter_fct_hx     , &
            flux_limiter_fct_tr

  real(r8), parameter :: one  = 1._r8
  real(r8), parameter :: zero = 0._r8

  contains

!================================================
!            FCT limiter on HX cells
!================================================

  subroutine flux_limiter_fct_hx(mesh,scalar_at_cell ,&
                                      ho_flux_at_edge,&
                                      up_flux_at_edge,&
                                      lm_flux_at_edge,&
                                      dtime)
! io
    type(global_domain), intent(in)        :: mesh
    real(r8), allocatable  , intent(in)    :: scalar_at_cell(:)     !halo(2)! scalar at time n
    real(r8), allocatable  , intent(in)    :: ho_flux_at_edge(:)    !halo(1)! high-order flux at edge, nr positive
    real(r8), allocatable  , intent(in)    :: up_flux_at_edge(:)    !halo(2)! low-order flux from upwind, nr positive
    real(r8), allocatable  , intent(inout) :: lm_flux_at_edge(:)    ! limited flux
    real(r8)               , intent(in)    :: dtime
! local
    real(r8), allocatable               :: df_flux_at_edge(:)    ! ho-up flux
    real(r8), allocatable               :: scalar_td_at_cell(:)  ! TD solution
    real(r8), allocatable               :: r_inflow_at_cell(:)   ! R
    real(r8), allocatable               :: r_outflow_at_cell(:)  ! R
    real(r8), dimension(1:3)            :: v0v1
    real(r8)                            :: div_sum
    real(r8)                            :: p_inflow, p_outflow
    real(r8)                            :: maf_inflow, maf_outflow
    real(r8)                            :: qmax, qmin
    real(r8)                            :: flag, coef
    integer(i4)                         :: cell_number, edge_number
    integer(i4)                         :: ineigh
    integer(i4)                         :: index_edge
    integer(i4)                         :: icell, ie, v0, v1, v_nb

    cell_number = mesh%nv
    edge_number = mesh%ne

    allocate(scalar_td_at_cell(1:mesh%nv_full))
    allocate(r_inflow_at_cell(1:mesh%nv_full))
    allocate(r_outflow_at_cell(1:mesh%nv_full))
    allocate(df_flux_at_edge(1:mesh%ne_full))

    do ie = 1, mesh%ne_halo(1)!edge_number
       df_flux_at_edge(ie) = (ho_flux_at_edge(ie)-up_flux_at_edge(ie))*rearth*mesh%edp(ie)%leng
    end do
!
! [1] preupdate TD solution at cell
!
    do icell = 1, mesh%nv_halo(2)!cell_number     ! for each vertice in the mesh
       div_sum = 0._r8
       do ineigh = 1, mesh%vtx(icell)%nnb
          index_edge  = mesh%vtx(icell)%ed(ineigh)
          div_sum     = div_sum+up_flux_at_edge(index_edge)*&
                                mesh%plg(icell)%nr(ineigh)*&
                                mesh%edp(index_edge)%leng*rearth
       end do
       scalar_td_at_cell(icell) = scalar_at_cell(icell)-dtime*div_sum/((rearth**2)*mesh%plg(icell)%areag)
    end do
!
! [2] compute R values at cell
!
    do icell = 1, mesh%nv_halo(1)!cell_number
       p_inflow  = zero 
       p_outflow = zero
       qmax      = max(scalar_at_cell(icell),scalar_td_at_cell(icell))
       qmin      = min(scalar_at_cell(icell),scalar_td_at_cell(icell))

       do ineigh = 1, mesh%vtx(icell)%nnb
          p_inflow  = p_inflow+min(zero ,mesh%plg(icell)%nr(ineigh)*df_flux_at_edge(mesh%vtx(icell)%ed(ineigh)))
          p_outflow = p_outflow+max(zero,mesh%plg(icell)%nr(ineigh)*df_flux_at_edge(mesh%vtx(icell)%ed(ineigh)))
          qmax      = max(qmax, scalar_at_cell(mesh%vtx(icell)%nb(ineigh)), scalar_td_at_cell(mesh%vtx(icell)%nb(ineigh)))
          qmin      = min(qmin, scalar_at_cell(mesh%vtx(icell)%nb(ineigh)), scalar_td_at_cell(mesh%vtx(icell)%nb(ineigh)))
       end do

       maf_inflow  = (qmax-scalar_td_at_cell(icell))*mesh%plg(icell)%areag*(rearth**2)/dtime
       maf_outflow = (scalar_td_at_cell(icell)-qmin)*mesh%plg(icell)%areag*(rearth**2)/dtime
         
       if(p_inflow.lt.zero)then
          r_inflow_at_cell(icell)  = min(one,maf_inflow/abs(p_inflow))
       else
          r_inflow_at_cell(icell)  = zero
       end if

       if(p_outflow.gt.dbleps)then
          r_outflow_at_cell(icell) = min(one,maf_outflow/abs(p_outflow))
       else
          r_outflow_at_cell(icell) = zero
       end if

    end do
!
![3] compute limited flux at edge
!
    do ie = 1, edge_number
       v0                = mesh%edt(ie)%v(1)
       v1                = mesh%edt(ie)%v(2)
       v0v1              = mesh%vtx(v1)%p-mesh%vtx(v0)%p
       flag              = sign(one,dot_product(v0v1,mesh%edp(ie)%nr))
       if(flag*df_flux_at_edge(ie).ge.0._r8)then ! flow is 0->1
          coef           = min(r_inflow_at_cell(v1),r_outflow_at_cell(v0))
       else  ! flow is 1->0
          coef           = min(r_inflow_at_cell(v0),r_outflow_at_cell(v1))
       end if
       lm_flux_at_edge(ie) = up_flux_at_edge(ie)+coef*(ho_flux_at_edge(ie)-up_flux_at_edge(ie))
    end do

    deallocate(r_inflow_at_cell)
    deallocate(r_outflow_at_cell)
    deallocate(scalar_td_at_cell)
    deallocate(df_flux_at_edge)

    return
  end subroutine flux_limiter_fct_hx

!================================================
!            FCT limiter on TR cells
!================================================

  subroutine flux_limiter_fct_tr(mesh,scalar_at_cell ,&
                                      ho_flux_at_edge,&
                                      up_flux_at_edge,&
                                      lm_flux_at_edge,&
                                      dtime)
!
! io
!
    type(global_domain),     intent(in)    :: mesh
    real(r8), allocatable  , intent(in)    :: scalar_at_cell(:)     ! scalar at time n
    real(r8), allocatable  , intent(in)    :: ho_flux_at_edge(:)    ! high-order flux at edge, nr positive
    real(r8), allocatable  , intent(in)    :: up_flux_at_edge(:)    ! low-order flux from upwind, nr positive
    real(r8), allocatable  , intent(inout) :: lm_flux_at_edge(:)    ! limited flux
    real(r8),                intent(in)    :: dtime
!
! local
!
    real(r8), allocatable               :: df_flux_at_edge(:)    ! ho-up flux
    real(r8), allocatable               :: scalar_td_at_cell(:)  ! TD solution
    real(r8), allocatable               :: r_inflow_at_cell(:)   ! R
    real(r8), allocatable               :: r_outflow_at_cell(:)  ! R
    real(r8), dimension(1:3)            :: v0v1
    real(r8)                            :: div_sum
    real(r8)                            :: p_inflow
    real(r8)                            :: p_outflow
    real(r8)                            :: maf_inflow
    real(r8)                            :: maf_outflow
    real(r8)                            :: qmax 
    real(r8)                            :: qmin
    real(r8)                            :: flag
    real(r8)                            :: coef
    integer(i4)                         :: cell_number
    integer(i4)                         :: edge_number
    integer(i4)                         :: ineigh
    integer(i4)                         :: index_edge
    integer(i4)                         :: icell
    integer(i4)                         :: ie
    integer(i4)                         :: v0
    integer(i4)                         :: v1
    integer(i4)                         :: v_nb 

    cell_number = mesh%nt
    edge_number = mesh%ne

    allocate(scalar_td_at_cell(1:cell_number))
    allocate( r_inflow_at_cell(1:cell_number))
    allocate(r_outflow_at_cell(1:cell_number))
    allocate(  df_flux_at_edge(1:edge_number))

    do ie = 1, edge_number
       df_flux_at_edge(ie) = (ho_flux_at_edge(ie)-up_flux_at_edge(ie))*rearth*mesh%edt(ie)%leng
    end do
!
![1] preupdate TD solution at cell
!
    do icell = 1, cell_number     ! for each vertice in the mesh
       div_sum = 0._r8
       do ineigh = 1, 3
          index_edge  = mesh%tri(icell)%ed(ineigh)
          div_sum     = div_sum+up_flux_at_edge(index_edge)*&
                                mesh%tri(icell)%nr(ineigh)*&
                                mesh%edt(index_edge)%leng*rearth
       end do
       scalar_td_at_cell(icell) = scalar_at_cell(icell)-dtime*div_sum/((rearth**2)*mesh%tri(icell)%areag)
    end do
!
![2] compute R values at cell
!
    do icell = 1, cell_number
       p_inflow  = 0._r8
       p_outflow = 0._r8
       qmax      = max(scalar_at_cell(icell), scalar_td_at_cell(icell))
       qmin      = min(scalar_at_cell(icell), scalar_td_at_cell(icell))

       do ineigh = 1, 3 
          p_inflow  = p_inflow+min(0._r8 ,mesh%tri(icell)%nr(ineigh)*df_flux_at_edge(mesh%tri(icell)%ed(ineigh)))
          p_outflow = p_outflow+max(0._r8,mesh%tri(icell)%nr(ineigh)*df_flux_at_edge(mesh%tri(icell)%ed(ineigh)))
          qmax      = max(qmax, scalar_at_cell(mesh%tri(icell)%nb(ineigh)), scalar_td_at_cell(mesh%tri(icell)%nb(ineigh)))
          qmin      = min(qmin, scalar_at_cell(mesh%tri(icell)%nb(ineigh)), scalar_td_at_cell(mesh%tri(icell)%nb(ineigh)))
       end do

       maf_inflow  = (qmax-scalar_td_at_cell(icell))*mesh%tri(icell)%areag*(rearth**2)/dtime
       maf_outflow = (scalar_td_at_cell(icell)-qmin)*mesh%tri(icell)%areag*(rearth**2)/dtime
         
       if(p_inflow.lt.0._r8)then
          r_inflow_at_cell(icell)  = min(1._r8,maf_inflow/abs(p_inflow))
       else
          r_inflow_at_cell(icell)  = 0._r8
       end if

       if(p_outflow.gt.0._r8)then
          r_outflow_at_cell(icell) = min(1._r8,maf_outflow/abs(p_outflow))
       else
          r_outflow_at_cell(icell) = 0._r8
       end if

    end do
!
![3] compute limited flux at edge
!
    do ie = 1, edge_number
       v0                = mesh%edp(ie)%v(1)
       v1                = mesh%edp(ie)%v(2)
       v0v1              = mesh%tri(v1)%c%p-mesh%tri(v0)%c%p
       flag              = sign(one,dot_product(v0v1,mesh%edt(ie)%nr))
       if(flag*df_flux_at_edge(ie).ge.0)then ! flow is 0->1
          coef           = min(r_inflow_at_cell(v1),r_outflow_at_cell(v0))
       else  ! flow is 1->0
          coef           = min(r_inflow_at_cell(v0),r_outflow_at_cell(v1))
       end if
       lm_flux_at_edge(ie) = up_flux_at_edge(ie)+coef*(ho_flux_at_edge(ie)-up_flux_at_edge(ie))
    end do

    deallocate(r_inflow_at_cell)
    deallocate(r_outflow_at_cell)
    deallocate(scalar_td_at_cell)
    deallocate(df_flux_at_edge)

    return
  end subroutine flux_limiter_fct_tr

  end module grist_limiter_module
