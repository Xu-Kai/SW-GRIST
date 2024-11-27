
!----------------------------------------------------------------------------
! Created on 2016
! Version 1.0
! Description: High-level flux operators used by solvers
! Revision history:
!----------------------------------------------------------------------------

 module grist_flux_operators

   use grist_constants       , only: gravity, i4, r8, rearth
   use grist_nml_module      , only: two_polyfit,one_polyfit,half_polyfit,&
                                     stencil_exchange_flag,lbdry_flag,stencil_width
   use grist_domain_types    , only: global_domain
   use grist_math_module     , only: norm
   use grist_base_flux_module, only: flux_upwind_nowind, flux_wrf34, flux_wrf3, flux_wrf56, flux_wrf5, flux_ppm

   use grist_ffsl_flux_module, only: calc_hx_normal_flux_mr07, &
                                     calc_hx_normal_flux_sm10, &
                                     calc_tr_normal_flux_mr07, &
                                     calc_tr_normal_flux_sm10

   use grist_recon_module    , only: calc_der_2nd_at_hx,&
                                     calc_double_der_2nd_at_hx,&
                                     calc_der_2nd_at_tr,&
                                     calc_der_4th_at_hx,&
                                     calc_der_4th_at_tr
   use grist_tspas_flux_module, only: calc_hx_normal_flux_tspas

  implicit none
  private

  public ::  calc_normal_flux_hx_edge, & ! contains q
             calc_normal_flux_tr_edge    ! contains q

  real(r8), parameter  :: one=1._r8

  contains

 subroutine calc_normal_flux_hx_edge(mesh,scalar_normal_velocity_at_edge  , &
                                          scalar_tangent_velocity_at_edge, & ! for ffsl
                                          scalar_height_at_prime_cell    , &
                                          scalar_normal_flux_at_edge     , &
                                          adv_flag                       , &
                                          dtime                          , &
                                          scalar_dptz_at_prime_cell      , &
                                          scalar_geop_at_prime_cell      , &
                                          called_by_ptadv)
! io
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge(:)
   real(r8), allocatable, intent(in)    :: scalar_tangent_velocity_at_edge(:)
   real(r8), allocatable, intent(in)    :: scalar_height_at_prime_cell(:)
   real(r8), allocatable, intent(inout) :: scalar_normal_flux_at_edge(:)
   integer(i4),           intent(in)    :: adv_flag
   real(r8)   ,           intent(in)    :: dtime
   real(r8), allocatable, optional, intent(in) :: scalar_dptz_at_prime_cell(:)
   real(r8), allocatable, optional, intent(in) :: scalar_geop_at_prime_cell(:)  ! actually height in meter
   logical,               optional, intent(in) :: called_by_ptadv
! local
   real(r8),dimension(1:3)              :: v0v1 !0->1
   real(r8)                             :: beta
   real(r8)                             :: wind_edge
   real(r8)                             :: der0, der1, der2, der3
   real(r8)                             :: der0_4th, der1_4th
   real(r8)                             :: qface1, qface2
   real(r8), allocatable                :: der_2nd(:,:)
   real(r8)                             :: local_der_2nd_center
   real(r8)                             :: local_der_2nd(10)
   integer(i4)                          :: v0, v1, v_upwind
   integer(i4)                          :: ie, inb, kk
   integer(i4)                          :: cell_index,edge_index, isum, flag, mark, kkk, tmp_ne, tmp_ne_halo

   select case(adv_flag)

   case (1)  ! 1st-order upwind

       do ie = 1, mesh%ne
          v0        = mesh%edt(ie)%v(1)       ! 1st end of edge
          v1        = mesh%edt(ie)%v(2)       ! 2nd end of edge
          !v0v1      = mesh%vtx(v1)%p-mesh%vtx(v0)%p
          !flag      = sign(one,dot_product(v0v1,mesh%edp(ie)%nr))
          !wind_edge = flag*scalar_normal_velocity_at_edge(ie)
          wind_edge = mesh%edt(ie)%edpNr_edtTg*scalar_normal_velocity_at_edge(ie)

          call flux_upwind_nowind(wind_edge,scalar_height_at_prime_cell(v0),&
                                            scalar_height_at_prime_cell(v1),&
                                            scalar_normal_flux_at_edge(ie))
       end do

   case(2)  ! 2nd-order center difference 

       do ie = 1, mesh%ne
          wind_edge                        = scalar_normal_velocity_at_edge(ie)
          v0                               = mesh%edt(ie)%v(1)
          v1                               = mesh%edt(ie)%v(2)
          scalar_normal_flux_at_edge(ie)   = 0.5*(scalar_height_at_prime_cell(v0)+&
                                                  scalar_height_at_prime_cell(v1))
       end do

   case(3,4,5,6,7)  ! WRF3/4, damping depends on beta value

    if(adv_flag.eq.3) beta=1._r8    ! 3rd order
    if(adv_flag.eq.4) beta=0        ! 4th order
    if(adv_flag.eq.5) beta=0.25     ! hybrid 1
    if(adv_flag.eq.6) beta=0.5      ! hybrid 2
    if(adv_flag.eq.7) beta=0.75     ! hybrid 3

    do ie = 1, mesh%ne

          !v0        = mesh%edt(ie)%v(1)       ! 1st end of edge
          !v1        = mesh%edt(ie)%v(2)       ! 2nd end of edge
!=========!=======================================
!      Mak!e sure wind follows C0C1 positive
!=========!=======================================
          !!v0v1      = mesh%vtx(v1)%p-mesh%vtx(v0)%p
          !!flag      = sign(one,dot_product(v0v1,mesh%edp(ie)%nr))
          !!wind_edge = flag*scalar_normal_velocity_at_edge(ie)
          !wind_edge = mesh%edt(ie)%edpNr_edtTg*scalar_normal_velocity_at_edge(ie)
      !        
          !call calc_der_2nd_at_hx(mesh,scalar_height_at_prime_cell,v0,ie,0,der0)
          !call calc_der_2nd_at_hx(mesh,scalar_height_at_prime_cell,v1,ie,1,der1)
          !call flux_wrf34(wind_edge, scalar_height_at_prime_cell(v0),&
          !                           scalar_height_at_prime_cell(v1),&
          !                der0,der1, scalar_normal_flux_at_edge(ie) ,&
          !                mesh%edt(ie)%leng,beta)

          v0        = mesh%edt_v(1, ie)       ! 1st end of edge
          v1        = mesh%edt_v(2, ie)       ! 2nd end of edge
!================================================
!      Make sure wind follows C0C1 positive
!================================================
          !v0v1      = mesh%vtx(v1)%p-mesh%vtx(v0)%p
          !flag      = sign(one,dot_product(v0v1,mesh%edp(ie)%nr))
          !wind_edge = flag*scalar_normal_velocity_at_edge(ie)
          wind_edge = mesh%edt_edpNr_edtTg(ie)*scalar_normal_velocity_at_edge(ie)

          call calc_der_2nd_at_hx(mesh,scalar_height_at_prime_cell,v0,ie,0,der0)
          call calc_der_2nd_at_hx(mesh,scalar_height_at_prime_cell,v1,ie,1,der1)
          call flux_wrf34(wind_edge, scalar_height_at_prime_cell(v0),&
                                     scalar_height_at_prime_cell(v1),&
                          der0,der1, scalar_normal_flux_at_edge(ie) ,&
                          mesh%edt_leng(ie),beta)

    end do

   case(33)  ! pure 3ord 

    do ie = 1, mesh%ne

          !v0        = mesh%edt(ie)%v(1)       ! 1st end of edge
          !v1        = mesh%edt(ie)%v(2)       ! 2nd end of edge
!=========!=======================================
!      Mak!e sure wind follows C0C1 positive
!=========!=======================================
          !!v0v1      = mesh%vtx(v1)%p-mesh%vtx(v0)%p
          !!flag      = sign(one,dot_product(v0v1,mesh%edp(ie)%nr))
          !!wind_edge = flag*scalar_normal_velocity_at_edge(ie)
          !wind_edge = mesh%edt(ie)%edpNr_edtTg*scalar_normal_velocity_at_edge(ie)
          !mark      = int(-0.5*sign(1._r8,wind_edge)+1.5)  ! 1 or 2
          !v_upwind  = mesh%edt(ie)%v(mark)    ! global index
  
          !call calc_der_2nd_at_hx(mesh,scalar_height_at_prime_cell,v_upwind,ie,mark-1,der0)

          !call flux_wrf3(scalar_height_at_prime_cell(v0),&
          !               scalar_height_at_prime_cell(v1),&
          !               der0, scalar_normal_flux_at_edge(ie) ,&
          !               mesh%edt(ie)%leng)

          v0        = mesh%edt_v(1, ie)       ! 1st end of edge
          v1        = mesh%edt_v(2, ie)       ! 2nd end of edge
!================================================
!      Make sure wind follows C0C1 positive
!================================================
          !v0v1      = mesh%vtx(v1)%p-mesh%vtx(v0)%p
          !flag      = sign(one,dot_product(v0v1,mesh%edp(ie)%nr))
          !wind_edge = flag*scalar_normal_velocity_at_edge(ie)
          wind_edge = mesh%edt_edpNr_edtTg(ie)*scalar_normal_velocity_at_edge(ie)
          mark      = int(-0.5*sign(1._r8,wind_edge)+1.5)  ! 1 or 2
          v_upwind  = mesh%edt_v(mark, ie)    ! global index
  
          call calc_der_2nd_at_hx(mesh,scalar_height_at_prime_cell,v_upwind,ie,mark-1,der0)

          call flux_wrf3(scalar_height_at_prime_cell(v0),&
                         scalar_height_at_prime_cell(v1),&
                         der0, scalar_normal_flux_at_edge(ie) ,&
                         mesh%edt_leng(ie))
    end do
   case(8,9,10,11,12)  ! WRF, 5th, 6th

    if(adv_flag.eq.8)  beta=1        ! 5rd order
    if(adv_flag.eq.9)  beta=0        ! 6th order
    if(adv_flag.eq.10) beta=0.25     ! hybrid
    if(adv_flag.eq.11) beta=0.5      ! hybrid
    if(adv_flag.eq.12) beta=0.75     ! hybrid

    if (.not. lbdry_flag) then 
       tmp_ne = mesh%ne
    else
       !tmp_ne = mesh%ne_bdry(stencil_width)
       tmp_ne = mesh%ne_bdry(2)
    end if

    do ie = 1, tmp_ne

          v0        = mesh%edt(ie)%v(1)       ! 1st end of edge
          v1        = mesh%edt(ie)%v(2)       ! 2nd end of edge
          !v0v1      = mesh%vtx(v1)%p-mesh%vtx(v0)%p
          !flag      = sign(one,dot_product(v0v1,mesh%edp(ie)%nr))
          !wind_edge = flag*scalar_normal_velocity_at_edge(ie)
          wind_edge = mesh%edt(ie)%edpNr_edtTg*scalar_normal_velocity_at_edge(ie)

          IF(two_polyfit)then ! two polyfit and hybrid equal to machine round off
             call calc_der_4th_at_hx(mesh,scalar_height_at_prime_cell,v0,ie,0,der0_4th,der0)
             call calc_der_4th_at_hx(mesh,scalar_height_at_prime_cell,v1,ie,1,der1_4th,der1)
             call calc_der_2nd_at_hx(mesh,scalar_height_at_prime_cell,v0,ie,0,der0)
             call calc_der_2nd_at_hx(mesh,scalar_height_at_prime_cell,v1,ie,1,der1)
             call flux_wrf56(wind_edge,scalar_height_at_prime_cell(v0),&
                                       scalar_height_at_prime_cell(v1),&
                             der0,der1,der0_4th,der1_4th,scalar_normal_flux_at_edge(ie),&
                             mesh%edt(ie)%leng,beta)
          ELSEIF(one_polyfit)then
             call calc_der_4th_at_hx(mesh,scalar_height_at_prime_cell,v0,ie,0,der0_4th,der0)
             call calc_der_4th_at_hx(mesh,scalar_height_at_prime_cell,v1,ie,1,der1_4th,der1)
             call flux_wrf56(wind_edge,scalar_height_at_prime_cell(v0),&
                                       scalar_height_at_prime_cell(v1),&
                             der0,der1,der0_4th,der1_4th,scalar_normal_flux_at_edge(ie),&
                             mesh%edt(ie)%leng,beta)
          ELSEIF(half_polyfit)then
             if(sign(1._r8,wind_edge).gt.0.)then ! v0 is upwind
                call calc_der_4th_at_hx(mesh,scalar_height_at_prime_cell,v0,ie,0,der0_4th,der0)
                call calc_der_2nd_at_hx(mesh,scalar_height_at_prime_cell,v1,ie,1,der1)
             else
                call calc_der_4th_at_hx(mesh,scalar_height_at_prime_cell,v1,ie,1,der0_4th,der0)
                call calc_der_2nd_at_hx(mesh,scalar_height_at_prime_cell,v0,ie,0,der1)
             end if

                call flux_wrf5(wind_edge,scalar_height_at_prime_cell(v0),&
                                         scalar_height_at_prime_cell(v1),&
                            der0,der1,der0_4th,scalar_normal_flux_at_edge(ie),&
                            mesh%edt(ie)%leng,beta)
          ELSE ! hybrid for O5
             if(sign(1._r8,wind_edge).gt.0.)then ! v0 is upwind
                call calc_der_4th_at_hx(mesh,scalar_height_at_prime_cell,v0,ie,0,der0_4th,der0)
             else
                call calc_der_4th_at_hx(mesh,scalar_height_at_prime_cell,v1,ie,1,der0_4th,der0)
             end if
                call calc_der_2nd_at_hx(mesh,scalar_height_at_prime_cell,v0,ie,0,der0)
                call calc_der_2nd_at_hx(mesh,scalar_height_at_prime_cell,v1,ie,1,der1)

                call flux_wrf5(wind_edge,scalar_height_at_prime_cell(v0),&
                                         scalar_height_at_prime_cell(v1),&
                            der0,der1,der0_4th,scalar_normal_flux_at_edge(ie),&
                            mesh%edt(ie)%leng,beta)
          END IF

    end do

    if (lbdry_flag) then 
      do ie = tmp_ne+1, mesh%ne

            v0        = mesh%edt(ie)%v(1)       ! 1st end of edge
            v1        = mesh%edt(ie)%v(2)       ! 2nd end of edge
!================================================
!      Make sure wind follows C0C1 positive
!================================================
            !v0v1      = mesh%vtx(v1)%p-mesh%vtx(v0)%p
            !flag      = sign(one,dot_product(v0v1,mesh%edp(ie)%nr))
            !wind_edge = flag*scalar_normal_velocity_at_edge(ie)
            wind_edge = mesh%edt(ie)%edpNr_edtTg*scalar_normal_velocity_at_edge(ie)

            call calc_der_2nd_at_hx(mesh,scalar_height_at_prime_cell,v0,ie,0,der0)
            call calc_der_2nd_at_hx(mesh,scalar_height_at_prime_cell,v1,ie,1,der1)
            call flux_wrf34(wind_edge, scalar_height_at_prime_cell(v0),&
                                       scalar_height_at_prime_cell(v1),&
                            der0,der1, scalar_normal_flux_at_edge(ie) ,&
                            mesh%edt(ie)%leng,beta)

      end do
    end if

   case(58)  ! a 5th order flux operator based on twice quadratic integration

    if(stencil_exchange_flag)then
      tmp_ne_halo = mesh%ne_halo(1)
    else
      tmp_ne_halo = mesh%ne_halo(2)
    endif
    allocate(der_2nd(tmp_ne_halo,2))

    do ie = 1, tmp_ne_halo
          !v0        = mesh%edt(ie)%v(1)       ! 1st end of edge
          !v1        = mesh%edt(ie)%v(2)       ! 2nd end of edge
          v0        = mesh%edt_v(1, ie)       ! 1st end of edge
          v1        = mesh%edt_v(2, ie)       ! 2nd end of edge
! calc 2nd-order derivative at two ends of edge ie, and store
          call calc_der_2nd_at_hx(mesh,scalar_height_at_prime_cell,v0,ie,0,der_2nd(ie,1))
          call calc_der_2nd_at_hx(mesh,scalar_height_at_prime_cell,v1,ie,1,der_2nd(ie,2))
    end do
!
    if (.not. lbdry_flag) then 
       tmp_ne = mesh%ne
    else
       !tmp_ne = mesh%ne_bdry(stencil_width)
       tmp_ne = mesh%ne_bdry(2)
    end if

    do ie = 1, tmp_ne

          !v0        = mesh%edt(ie)%v(1)       ! 1st end of edge
          !v1        = mesh%edt(ie)%v(2)       ! 2nd end of edge
          !!v0v1      = mesh%vtx(v1)%p-mesh%vtx(v0)%p
          !!flag      = sign(one,dot_product(v0v1,mesh%edp(ie)%nr))
          !!wind_edge = flag*scalar_normal_velocity_at_edge(ie)
          !wind_edge = mesh%edt(ie)%edpNr_edtTg*scalar_normal_velocity_at_edge(ie)
          !mark      = int(-0.5*sign(1._r8,wind_edge)+1.5)  ! 1 or 2
          !v_upwind  = mesh%edt(ie)%v(mark)    ! global index
          !local_der_2nd_center = der_2nd(ie,mark)
! for v0(v!1) , we select 2nd-order derivatives connecting V0V1, and VxV0(X=1,nnb)
! this seq!uence is same as that when constructing the B matrix
          !do inb = 1, mesh%vtx(v_upwind)%nnb
          !   edge_index  = mesh%edt(ie)%my_edge_on_edge(mark,inb)
          !   cell_index  = mesh%edt(ie)%ur_cell_on_edge(mark,inb)
          !   local_der_2nd(inb) = der_2nd(edge_index,cell_index)
          !end do

          !call calc_double_der_2nd_at_hx(mesh,local_der_2nd_center,local_der_2nd,&
          !                                    mesh%vtx(v_upwind)%nnb, &
          !                                    v_upwind,ie,mark-1,der0_4th)
          !call flux_wrf5(wind_edge,scalar_height_at_prime_cell(v0)   , &
          !                         scalar_height_at_prime_cell(v1)   , &
          !                   der_2nd(ie,1),der_2nd(ie,2),der0_4th    , &
          !                   scalar_normal_flux_at_edge(ie), mesh%edt(ie)%leng,beta)

          v0        = mesh%edt_v(1, ie)       ! 1st end of edge
          v1        = mesh%edt_v(2, ie)       ! 2nd end of edge
          !v0v1      = mesh%vtx(v1)%p-mesh%vtx(v0)%p
          !flag      = sign(one,dot_product(v0v1,mesh%edp(ie)%nr))
          !wind_edge = flag*scalar_normal_velocity_at_edge(ie)
          wind_edge = mesh%edt_edpNr_edtTg(ie)*scalar_normal_velocity_at_edge(ie)
          mark      = int(-0.5*sign(1._r8,wind_edge)+1.5)  ! 1 or 2
          v_upwind  = mesh%edt_v(mark, ie)    ! global index
          local_der_2nd_center = der_2nd(ie,mark)
! for v0(v1) , we select 2nd-order derivatives connecting V0V1, and VxV0(X=1,nnb)
! this sequence is same as that when constructing the B matrix
          do inb = 1, mesh%vtx_nnb(v_upwind)
             edge_index  = mesh%edt_my_edge_on_edge(mark,inb,ie)
             cell_index  = mesh%edt_ur_cell_on_edge(mark,inb,ie)
             local_der_2nd(inb) = der_2nd(edge_index,cell_index)
          end do

          call calc_double_der_2nd_at_hx(mesh,local_der_2nd_center,local_der_2nd,&
                                              mesh%vtx_nnb(v_upwind), &
                                              v_upwind,ie,mark-1,der0_4th)

          call flux_wrf5(wind_edge,scalar_height_at_prime_cell(v0)   , &
                                   scalar_height_at_prime_cell(v1)   , &
                             der_2nd(ie,1),der_2nd(ie,2),der0_4th    , &
                             scalar_normal_flux_at_edge(ie), mesh%edt_leng(ie),beta)
    end do
    if (lbdry_flag) then 
      do ie = tmp_ne+1, mesh%ne

          !v0        = mesh%edt(ie)%v(1)       ! 1st end of edge
          !v1        = mesh%edt(ie)%v(2)       ! 2nd end of edge
!=========!=======================================
!      Mak!e sure wind follows C0C1 positive
!=========!=======================================
          !!v0v1      = mesh%vtx(v1)%p-mesh%vtx(v0)%p
          !!flag      = sign(one,dot_product(v0v1,mesh%edp(ie)%nr))
          !!wind_edge = flag*scalar_normal_velocity_at_edge(ie)
          !wind_edge = mesh%edt(ie)%edpNr_edtTg*scalar_normal_velocity_at_edge(ie)

          !!call calc_der_2nd_at_hx(mesh,scalar_height_at_prime_cell,v0,ie,0,der0)
          !!call calc_der_2nd_at_hx(mesh,scalar_height_at_prime_cell,v1,ie,1,der1)
          !call flux_wrf34(wind_edge, scalar_height_at_prime_cell(v0),&
          !                           scalar_height_at_prime_cell(v1),&
          !                der_2nd(ie,1),der_2nd(ie,1), scalar_normal_flux_at_edge(ie) ,&
          !                mesh%edt(ie)%leng,beta)

          v0        = mesh%edt_v(1, ie)       ! 1st end of edge
          v1        = mesh%edt_v(2, ie)       ! 2nd end of edge
!================================================
!      Make sure wind follows C0C1 positive
!================================================
          !v0v1      = mesh%vtx(v1)%p-mesh%vtx(v0)%p
          !flag      = sign(one,dot_product(v0v1,mesh%edp(ie)%nr))
          !wind_edge = flag*scalar_normal_velocity_at_edge(ie)
          wind_edge = mesh%edt_edpNr_edtTg(ie)*scalar_normal_velocity_at_edge(ie)

          !call calc_der_2nd_at_hx(mesh,scalar_height_at_prime_cell,v0,ie,0,der0)
          !call calc_der_2nd_at_hx(mesh,scalar_height_at_prime_cell,v1,ie,1,der1)
          call flux_wrf34(wind_edge, scalar_height_at_prime_cell(v0),&
                                     scalar_height_at_prime_cell(v1),&
                          der_2nd(ie,1),der_2nd(ie,1), scalar_normal_flux_at_edge(ie) ,&
                          mesh%edt_leng(ie),beta)
      end do
    end if

    deallocate(der_2nd)

   case(27)  ! MS13, ULA
      
      call calc_hx_normal_flux_mr07(mesh,scalar_normal_velocity_at_edge ,&
                                         scalar_tangent_velocity_at_edge,&
                                         scalar_height_at_prime_cell    ,&
                                         scalar_normal_flux_at_edge     ,&
                                         adv_flag                       ,&
                                         dtime)

   case(28)  ! MS13, UQA-1
      
      call calc_hx_normal_flux_sm10(mesh,scalar_normal_velocity_at_edge ,&
                                         scalar_tangent_velocity_at_edge,&
                                         scalar_height_at_prime_cell    ,&
                                         scalar_normal_flux_at_edge     ,&
                                         adv_flag                       ,&
                                         dtime)

   case(99) ! TSPAS
      call calc_hx_normal_flux_tspas(mesh,scalar_normal_velocity_at_edge,&
                                          scalar_height_at_prime_cell   ,&
                                          scalar_normal_flux_at_edge    ,&
                                          dtime                        )

   case default
        print*,"you must select an advection scheme/flux divergence operator"
   end select

   return
  end subroutine calc_normal_flux_hx_edge

!================================================
!             only contain q
!================================================

  subroutine calc_normal_flux_tr_edge(mesh,scalar_pv_at_dual_cell         ,& ! dual-cell value 
                                           scalar_normal_velocity_at_edge ,& ! normal wind at tr edge
                                           scalar_tangent_velocity_at_edge,& ! tangent wind at tr edge, for ffsl
                                           scalar_pv_at_edge              ,& ! edge-based value mapped from dual cell
                                           called_by_tr                   ,& ! logical
                                           pv_flag                        ,&
                                           dtime       )
! io
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_pv_at_dual_cell(:)
   real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge(:)
   real(r8), allocatable, intent(in)    :: scalar_tangent_velocity_at_edge(:)
   real(r8), allocatable, intent(inout) :: scalar_pv_at_edge(:)
   logical,               intent(in)    :: called_by_tr
   integer(i4),           intent(in)    :: pv_flag
   real(r8),              intent(in)    :: dtime
! local
   real(r8),dimension(1:3)              :: v0v1            ! direction 0->1
   real(r8)                             :: beta
   real(r8)                             :: der0
   real(r8)                             :: der1
   real(r8)                             :: der0_4th
   real(r8)                             :: der1_4th
   real(r8)                             :: flag
   real(r8)                             :: wind_edge
   real(r8)                             :: qface1
   real(r8)                             :: qface2
   integer(i4)                          :: ie
   integer(i4)                          :: v0
   integer(i4)                          :: v1, mark, v_upwind

!================================================
!       potential vorticity at hx's edge
!================================================

      select case(pv_flag)

      case (1)  ! 1st-order upwind
        if(.not. called_by_tr) then
          do ie = 1, mesh%ne
             v0        = mesh%edp(ie)%v(1)           ! 1st end of edge
             v1        = mesh%edp(ie)%v(2)           ! 2nd end of edge
             !v0v1      = mesh%tri(v1)%c%p-mesh%tri(v0)%c%p

             !flag      = sign(one,dot_product(v0v1,mesh%edp(ie)%tg))
   !
   ! called by SWE solver on HX for PV
   !
             wind_edge = mesh%edp(ie)%edp12_edp_tg*scalar_normal_velocity_at_edge(ie)
             !wind_edge = flag*scalar_normal_velocity_at_edge(ie)
   
             call flux_upwind_nowind(wind_edge,scalar_pv_at_dual_cell(v0),&
                                               scalar_pv_at_dual_cell(v1),&
                                               scalar_pv_at_edge(ie))
   !
   ! substract wind_edge
   !
             scalar_pv_at_edge(ie)  = scalar_pv_at_edge(ie)
          end do
        else
          do ie = 1, mesh%ne
             v0        = mesh%edp(ie)%v(1)           ! 1st end of edge
             v1        = mesh%edp(ie)%v(2)           ! 2nd end of edge
             !v0v1      = mesh%tri(v1)%c%p-mesh%tri(v0)%c%p

             !flag      = sign(one,dot_product(v0v1,mesh%edt(ie)%nr))
   !
   ! called by ADV solver on TR
   !
             wind_edge = mesh%edp(ie)%edp12_edt_nr*scalar_normal_velocity_at_edge(ie)
             !wind_edge = flag*scalar_normal_velocity_at_edge(ie)
   
             call flux_upwind_nowind(wind_edge,scalar_pv_at_dual_cell(v0),&
                                               scalar_pv_at_dual_cell(v1),&
                                               scalar_pv_at_edge(ie))
   !
   ! substract wind_edge
   !
             scalar_pv_at_edge(ie)  = scalar_pv_at_edge(ie)
          end do
        end if
      case(2) ! standard 2nd in R10
        do ie = 1, mesh%ne
             v0 = mesh%edp(ie)%v(1)
             v1 = mesh%edp(ie)%v(2)
             scalar_pv_at_edge(ie) = 0.5*(scalar_pv_at_dual_cell(v0)+&
                                            scalar_pv_at_dual_cell(v1))
        end do
      case(3,4,5,6,7)  ! WRF, 3rd/4th

        if(pv_flag.eq.3) beta=1        ! 3rd order
        if(pv_flag.eq.4) beta=0        ! 4th order
        if(pv_flag.eq.5) beta=0.25     ! hybrid 1
        if(pv_flag.eq.6) beta=0.5      ! hybrid 2
        if(pv_flag.eq.7) beta=0.75     ! hybrid 3

        if(.not. called_by_tr) then
          !do ie = 1, mesh%ne
          !   v0        = mesh%edp(ie)%v(1)       ! 1st end of edge
          !   v1        = mesh%edp(ie)%v(2)       ! 2nd end of edge
          !   !v0v1      = mesh%tri(v1)%c%p-mesh%tri(v0)%c%p
          !   !flag      = sign(one,dot_product(v0v1,mesh%edp(ie)%tg))
  !
  ! called! by SWE solver on HX for PV
  !
          !   wind_edge = mesh%edp(ie)%edp12_edp_tg*scalar_normal_velocity_at_edge(ie)
          !   !wind_edge = flag*scalar_normal_velocity_at_edge(ie)
          !   mark      = int(-0.5*sign(1._r8,wind_edge)+1.5)  ! 1 or 2
          !   v_upwind  = mesh%edp(ie)%v(mark)    ! global index
          !   call calc_der_2nd_at_tr(mesh,scalar_pv_at_dual_cell,v_upwind,ie,mark-1,der0)
          !   call flux_wrf3(scalar_pv_at_dual_cell(v0),&
          !                  scalar_pv_at_dual_cell(v1),&
          !                  der0,scalar_pv_at_edge(ie),         &
          !                  mesh%edp(ie)%leng)
          !end do
          do ie = 1, mesh%ne
             v0        = mesh%edp_v(1, ie)       ! 1st end of edge
             v1        = mesh%edp_v(2, ie)       ! 2nd end of edge
             !v0v1      = mesh%tri(v1)%c%p-mesh%tri(v0)%c%p
             !flag      = sign(one,dot_product(v0v1,mesh%edp(ie)%tg))
  !
  ! called by SWE solver on HX for PV
  !
             wind_edge = mesh%edp12_edp_tg(ie)*scalar_normal_velocity_at_edge(ie)
             !wind_edge = flag*scalar_normal_velocity_at_edge(ie)
             mark      = int(-0.5*sign(1._r8,wind_edge)+1.5)  ! 1 or 2
             v_upwind  = mesh%edp_v(mark, ie)    ! global index
             call calc_der_2nd_at_tr(mesh,scalar_pv_at_dual_cell,v_upwind,ie,mark-1,der0)
             call flux_wrf3(scalar_pv_at_dual_cell(v0),&
                            scalar_pv_at_dual_cell(v1),&
                            der0,scalar_pv_at_edge(ie),         &
                            mesh%edp_leng(ie))
          end do
        else
          !do ie = 1, mesh%ne
          !   v0        = mesh%edp(ie)%v(1)       ! 1st end of edge
          !   v1        = mesh%edp(ie)%v(2)       ! 2nd end of edge
          !   !v0v1      = mesh%tri(v1)%c%p-mesh%tri(v0)%c%p
          !   !flag      = sign(one,dot_product(v0v1,mesh%edt(ie)%nr))
  !
  ! called! by ADV solver on TR
  !
          !   wind_edge = mesh%edp(ie)%edp12_edt_nr*scalar_normal_velocity_at_edge(ie)
          !   !wind_edge = flag*scalar_normal_velocity_at_edge(ie)
          !   mark      = int(-0.5*sign(1._r8,wind_edge)+1.5)  ! 1 or 2
          !   v_upwind  = mesh%edp(ie)%v(mark)    ! global index
          !   call calc_der_2nd_at_tr(mesh,scalar_pv_at_dual_cell,v_upwind,ie,mark-1,der0)
          !   call flux_wrf3(scalar_pv_at_dual_cell(v0),&
          !                  scalar_pv_at_dual_cell(v1),&
          !                  der0,scalar_pv_at_edge(ie),         &
          !                  mesh%edp(ie)%leng)
          !end do
          do ie = 1, mesh%ne
             v0        = mesh%edp_v(1, ie)       ! 1st end of edge
             v1        = mesh%edp_v(2, ie)       ! 2nd end of edge
             !v0v1      = mesh%tri(v1)%c%p-mesh%tri(v0)%c%p
             !flag      = sign(one,dot_product(v0v1,mesh%edt(ie)%nr))
  !
  ! called by ADV solver on TR
  !
             wind_edge = mesh%edp12_edt_nr(ie)*scalar_normal_velocity_at_edge(ie)
             !wind_edge = flag*scalar_normal_velocity_at_edge(ie)
             mark      = int(-0.5*sign(1._r8,wind_edge)+1.5)  ! 1 or 2
             v_upwind  = mesh%edp_v(mark, ie)    ! global index
             call calc_der_2nd_at_tr(mesh,scalar_pv_at_dual_cell,v_upwind,ie,mark-1,der0)
             call flux_wrf3(scalar_pv_at_dual_cell(v0),&
                            scalar_pv_at_dual_cell(v1),&
                            der0,scalar_pv_at_edge(ie),         &
                            mesh%edp_leng(ie))
          end do
        end if

      case(8,9,10,11,12)  ! WRF, 5rd/6th

        if(pv_flag.eq.8)  beta=1        ! 5rd order
        if(pv_flag.eq.9)  beta=0        ! 6th order
        if(pv_flag.eq.10) beta=0.25     ! hybrid
        if(pv_flag.eq.11) beta=0.5      ! hybrid
        if(pv_flag.eq.12) beta=0.75     ! hybrid

        if(.not. called_by_tr) then
          do ie = 1, mesh%ne

             v0        = mesh%edp(ie)%v(1)  ! 1st end of edge
             v1        = mesh%edp(ie)%v(2)  ! 2nd end of edge
             !v0v1      = mesh%tri(v1)%c%p-mesh%tri(v0)%c%p

             !flag      = sign(one,dot_product(v0v1,mesh%edp(ie)%tg))
  !
  ! called by SWE solver on HX for PV
  !
             wind_edge = mesh%edp(ie)%edp12_edp_tg*scalar_normal_velocity_at_edge(ie)
             !wind_edge = flag*scalar_normal_velocity_at_edge(ie)
  
             call calc_der_2nd_at_tr(mesh,scalar_pv_at_dual_cell,v0,ie,0,der0)
             call calc_der_2nd_at_tr(mesh,scalar_pv_at_dual_cell,v1,ie,1,der1)
             call calc_der_4th_at_tr(mesh,scalar_pv_at_dual_cell,v0,ie,0,der0_4th)
             call calc_der_4th_at_tr(mesh,scalar_pv_at_dual_cell,v1,ie,1,der1_4th)
  
  !================================================
  !               Compute edge-based value(flux)
  !================================================
  
             call flux_wrf56(wind_edge,scalar_pv_at_dual_cell(v0),&
                                       scalar_pv_at_dual_cell(v1),&
                                       der0,der1,der0_4th, der1_4th,&
                             scalar_pv_at_edge(ie), mesh%edp(ie)%leng,beta)
          end do
        else
          do ie = 1, mesh%ne

             v0        = mesh%edp(ie)%v(1)  ! 1st end of edge
             v1        = mesh%edp(ie)%v(2)  ! 2nd end of edge
             !v0v1      = mesh%tri(v1)%c%p-mesh%tri(v0)%c%p

             !flag      = sign(one,dot_product(v0v1,mesh%edt(ie)%nr))
  !
  ! called by ADV solver on TR
  !
             wind_edge = mesh%edp(ie)%edp12_edt_nr*scalar_normal_velocity_at_edge(ie)
             !wind_edge = flag*scalar_normal_velocity_at_edge(ie)
  
             call calc_der_2nd_at_tr(mesh,scalar_pv_at_dual_cell,v0,ie,0,der0)
             call calc_der_2nd_at_tr(mesh,scalar_pv_at_dual_cell,v1,ie,1,der1)
             call calc_der_4th_at_tr(mesh,scalar_pv_at_dual_cell,v0,ie,0,der0_4th)
             call calc_der_4th_at_tr(mesh,scalar_pv_at_dual_cell,v1,ie,1,der1_4th)
  
  !================================================
  !               Compute edge-based value(flux)
  !================================================
  
             call flux_wrf56(wind_edge,scalar_pv_at_dual_cell(v0),&
                                       scalar_pv_at_dual_cell(v1),&
                                       der0,der1,der0_4th, der1_4th,&
                             scalar_pv_at_edge(ie), mesh%edp(ie)%leng,beta)
          end do
        end if
      case (27)  ! MS13, ULA
        call calc_tr_normal_flux_mr07(mesh,scalar_normal_velocity_at_edge ,&
                                           scalar_tangent_velocity_at_edge,&
                                           scalar_pv_at_dual_cell         ,&
                                           scalar_pv_at_edge              ,&
                                           called_by_tr                   ,&
                                           dtime)
      case (28)  ! MS13, UQA-I
        call calc_tr_normal_flux_sm10(mesh,scalar_normal_velocity_at_edge ,&
                                           scalar_tangent_velocity_at_edge,&
                                           scalar_pv_at_dual_cell         ,&
                                           scalar_pv_at_edge              ,&
                                           called_by_tr                   ,&
                                           dtime)
      case default
         print*,"PLEASE SET PV_ORDER IN NAMELIST!!! GRIST STOPS..."
         stop
      end select

     return
  end subroutine calc_normal_flux_tr_edge

  end module grist_flux_operators
