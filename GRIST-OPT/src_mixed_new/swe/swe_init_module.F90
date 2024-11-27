
!----------------------------------------------------------------------------
! Created on 2016
! Version 1.0
! Description: Initialize SWE model
! Revision history:
!----------------------------------------------------------------------------

 module swe_init_module

  use grist_constants,       only: gravity,  pi, i4, r8, rearth,day2sec,omega
  use grist_domain_types,    only: global_domain
  use grist_data_types,      only: scalar_1d_field, vector_type
  use grist_math_module,     only: arcdistll,  convert_vector_sph2cart, &
                                   sph2cart, cart2sph, tri_kite_areas
  use grist_nml_module,      only: testcase, initialfield, modon_radius, scalar_type,&
                                   use_streamf
  use grist_data_types,      only: exchange_field_list_2d, scalar_2d_field
  use grist_config_partition,only: exchange_data_2d_add, exchange_data_2d

  implicit none

  private

  public   :: swe_initialize,       &
              swe_initialize_tr,    &
              overwrite_wind,       &
              overwrite_wind_streamf

  contains

!================================================
!                init SWE
!================================================

  subroutine swe_initialize(mesh,scalar_normal_velocity_at_edge,&
                                 scalar_height_at_prime_cell   ,&
                                 scalar_topo_at_prime_cell     ,&
                                 scalar_u_at_edge              ,&
                                 scalar_v_at_edge)

! io
   type(global_domain),     intent(inout) :: mesh
   type(scalar_1d_field),   intent(inout) :: scalar_normal_velocity_at_edge
   type(scalar_1d_field),   intent(inout) :: scalar_height_at_prime_cell
   type(scalar_1d_field),   intent(inout) :: scalar_topo_at_prime_cell
   type(scalar_1d_field),   intent(inout) :: scalar_u_at_edge      ! U wind at initial time
   type(scalar_1d_field),   intent(inout) :: scalar_v_at_edge      ! V ..............
! local
   type(vector_type)             :: vector_velocity_at_edge
   type(scalar_1d_field)               :: stream_function_dual_cell
   real(r8)                            :: lon,lat, time
   real(r8)                            :: sub_h_value,sub_t_value, sub_area, tmp_h_value, tmp_t_value
   integer(i4)                         :: it, ie, iv, inb
   integer(i4)                         :: v0, v1 
   real(r8)                            :: flag,v0v1(3)
   type(scalar_2d_field)               :: tmp_tri_kite_area
   type(exchange_field_list_2d),pointer:: field_head_2d

   time = 0._r8 ! only for initial conditions

!================================================
!         set scalar initial conditions
!================================================

   if(trim(scalar_type).eq.'point')then
      do iv = 1, mesh%nv
         scalar_topo_at_prime_cell%f(iv)   =   ini_topo(mesh%vtx(iv)%lon, mesh%vtx(iv)%lat)
         scalar_height_at_prime_cell%f(iv) = ini_scalar(mesh%vtx(iv)%lon, mesh%vtx(iv)%lat, time)-scalar_topo_at_prime_cell%f(iv)
      end do
   end if

   if(trim(scalar_type).eq.'cell')then
      do iv = 1, mesh%nv
         sub_t_value = 0._r8
         sub_h_value = 0._r8
         sub_area    = 0._r8 
         do inb = 1, mesh%vtx(iv)%nnb
            do ie  = 1, 3
               call cart2sph(mesh%plg(iv)%sub_triangle_midp_3d(inb,ie,1),&
                             mesh%plg(iv)%sub_triangle_midp_3d(inb,ie,2),&
                             mesh%plg(iv)%sub_triangle_midp_3d(inb,ie,3),&
                             lon, lat)
               tmp_t_value =   ini_topo(lon, lat)
               tmp_h_value = ini_scalar(lon, lat, time)
               sub_t_value = sub_t_value+tmp_t_value*mesh%plg(iv)%sub_triangle_area(inb)/3._r8
               sub_h_value = sub_h_value+tmp_h_value*mesh%plg(iv)%sub_triangle_area(inb)/3._r8
            end do
            sub_area  = sub_area+mesh%plg(iv)%sub_triangle_area(inb)
         end do
         scalar_topo_at_prime_cell%f(iv)     = sub_t_value/sub_area
         scalar_height_at_prime_cell%f(iv)   = sub_h_value/sub_area
      end do
   end if

!================================================
! Full velocity vector at intersecting point
! between tr x hx
!================================================

   vector_velocity_at_edge%pos = 6

   if(.not.allocated(vector_velocity_at_edge%p))then
      allocate(vector_velocity_at_edge%p(1:mesh%ne))
   end if

!================================================
!                normal velocity 
!================================================

   call ini_vector_velocity(mesh,vector_velocity_at_edge,&
                                 time                   ,&
                                 scalar_u_at_edge       ,&
                                 scalar_v_at_edge)

   do ie=1, mesh%ne
      scalar_normal_velocity_at_edge%f(ie) = dot_product(vector_velocity_at_edge%p(ie)%v, mesh%edp(ie)%nr)
   end do

   if(use_streamf)then ! overwrite normal wind using stream function
      if(.not.allocated(stream_function_dual_cell%f)) allocate(stream_function_dual_cell%f(mesh%nt_full))
      do it = 1, mesh%nt_full
         stream_function_dual_cell%f(it)  = ini_stream_function(mesh%tri(it)%c%lon,mesh%tri(it)%c%lat,0._r8)
      end do

      do ie = 1, mesh%ne
         v0   = mesh%edp(ie)%v(1) ! tri point
         v1   = mesh%edp(ie)%v(2) ! tri point
         v0v1 = mesh%tri(v1)%c%p-mesh%tri(v0)%c%p
         flag = sign(1._r8, dot_product(v0v1,mesh%edp(ie)%tg))
         scalar_normal_velocity_at_edge%f(ie) = flag*(stream_function_dual_cell%f(v1)-stream_function_dual_cell%f(v0))/(rearth*mesh%edp(ie)%leng)
      end do
      deallocate(stream_function_dual_cell%f)
   end if

!================================================
!           calculate intersecting area
!================================================

   !call tri_kite_areas(mesh)
   field_head_2d=>null()
   
   mesh%nt = mesh%nt_compute
   call tri_kite_areas(mesh)
   mesh%nt = mesh%nt_full
   
   allocate(tmp_tri_kite_area%f(3,mesh%nt_full))
   do it=1,mesh%nt_compute
     tmp_tri_kite_area%f(1:3,it) = mesh%tri(it)%kite_area(1:3)
   end do
   
   call exchange_data_2d_add(mesh,field_head_2d,tmp_tri_kite_area)
   call exchange_data_2d(mesh%local_block,field_head_2d)

   do it=mesh%nt_compute+1,mesh%nt_full
     mesh%tri(it)%kite_area(1:3) = tmp_tri_kite_area%f(1:3,it)
   end do
   deallocate(tmp_tri_kite_area%f)
   
   deallocate(vector_velocity_at_edge%p)

   return
  end subroutine swe_initialize

!================================================
!              init SWE on TR
!================================================

  subroutine swe_initialize_tr(mesh,scalar_normal_velocity_at_tr_edge,&
                                    scalar_height_at_dual_cell, &
                                    scalar_u_at_edge,     &
                                    scalar_v_at_edge      ) 

! io
   type(global_domain),     intent(inout) :: mesh
   type(scalar_1d_field),   intent(inout) :: scalar_normal_velocity_at_tr_edge
   type(scalar_1d_field),   intent(inout) :: scalar_height_at_dual_cell
   type(scalar_1d_field),   intent(inout) :: scalar_u_at_edge
   type(scalar_1d_field),   intent(inout) :: scalar_v_at_edge
! local
   type(vector_type)                   :: vector_velocity_at_edge
   real(r8)                            :: time
   integer(i4)                         :: it
   integer(i4)                         :: ie

   time = 0._r8 ! only for initial conditions

!================================================
!         set scalar initial conditions
!================================================

   do it=1,mesh%nt
      scalar_height_at_dual_cell%f(it) = ini_scalar(mesh%tri(it)%c%lon,mesh%tri(it)%c%lat)
   end do

!================================================
! Full velocity vector at intersecting point 
! between tr x hx
!================================================

   vector_velocity_at_edge%pos = 2

   if(.not.allocated(vector_velocity_at_edge%p))then
      allocate(vector_velocity_at_edge%p(1:mesh%ne))
   end if

!================================================
!                normal velocity 
!================================================


   call ini_vector_velocity(mesh,vector_velocity_at_edge,&
                                 time                   ,&
                                 scalar_u_at_edge       ,&
                                 scalar_v_at_edge)

   do ie=1, mesh%ne
      scalar_normal_velocity_at_tr_edge%f(ie) = dot_product(vector_velocity_at_edge%p(ie)%v, mesh%edt(ie)%nr)
   end do

   deallocate(vector_velocity_at_edge%p)

   return
  end subroutine swe_initialize_tr

!================================================
!   overwrite normal wind for adv on hx and tr 
!================================================

  subroutine overwrite_wind(mesh, normal_wind, tangent_wind, itimestep, timestep,is_tr)
! io
   type(global_domain),  intent(in)    :: mesh
   type(scalar_1d_field),intent(inout) :: normal_wind
   type(scalar_1d_field),intent(inout) :: tangent_wind
   integer(i4)       ,   intent(in)    :: itimestep
   real(r8),             intent(in)    :: timestep
   logical,              intent(in)    :: is_tr
! local
   type(vector_type)       :: vector_velocity_at_edge
   type(scalar_1d_field)   :: tmp_u_at_edge ! not used, just for storing
   type(scalar_1d_field)   :: tmp_v_at_edge ! not used, just for storing
   real(r8)                :: time
   integer(i4)             :: ie


   allocate(vector_velocity_at_edge%p(1:mesh%ne_full))
   allocate(tmp_u_at_edge%f(1:mesh%ne_full))
   allocate(tmp_v_at_edge%f(1:mesh%ne_full))

   vector_velocity_at_edge%pos           = 6 ! hx's edge
   if(is_tr) vector_velocity_at_edge%pos = 2 ! tr's edge

   time = real(itimestep,r8)*timestep-timestep ! u(n-1)

   call ini_vector_velocity(mesh,vector_velocity_at_edge, time, tmp_u_at_edge, tmp_v_at_edge)

   do ie=1, mesh%ne
      normal_wind%f(ie)           = dot_product(vector_velocity_at_edge%p(ie)%v, mesh%edp(ie)%nr)  ! for plg
      tangent_wind%f(ie)          = dot_product(vector_velocity_at_edge%p(ie)%v, mesh%edp(ie)%tg)  ! for plg
      if(is_tr) normal_wind%f(ie) = dot_product(vector_velocity_at_edge%p(ie)%v, mesh%edt(ie)%nr)  ! for tri
      if(is_tr) tangent_wind%f(ie)= dot_product(vector_velocity_at_edge%p(ie)%v, mesh%edt(ie)%tg)  ! for tri
   end do

   deallocate(vector_velocity_at_edge%p)
   deallocate(tmp_u_at_edge%f)
   deallocate(tmp_v_at_edge%f)
   
   return
  end subroutine overwrite_wind

  subroutine overwrite_wind_streamf(mesh, normal_wind, tangent_wind, itimestep, timestep,is_tr)
! io
   type(global_domain),  intent(in)    :: mesh
   type(scalar_1d_field),intent(inout) :: normal_wind
   type(scalar_1d_field),intent(inout) :: tangent_wind
   integer(i4)       ,   intent(in)    :: itimestep
   real(r8),             intent(in)    :: timestep
   logical,              intent(in)    :: is_tr
! local
   type(scalar_1d_field)   :: stream_function_prime_cell
   type(scalar_1d_field)   :: stream_function_dual_cell
   real(r8)                :: time,v0v1(3),v2v3(3)
   real(r8)                :: flag
   integer(i4)             :: v0, v1, v2, v3
   integer(i4)             :: it, ie, iv

   time = real(itimestep,r8)*timestep-timestep ! u(n-1)

   allocate(stream_function_prime_cell%f(mesh%nv_full))
   allocate(stream_function_dual_cell%f(mesh%nt_full))

   do iv = 1, mesh%nv_full
      stream_function_prime_cell%f(iv) = ini_stream_function(mesh%vtx(iv)%lon,mesh%vtx(iv)%lat,time)
   end do 

   do it = 1, mesh%nt_full
      stream_function_dual_cell%f(it)  = ini_stream_function(mesh%tri(it)%c%lon,mesh%tri(it)%c%lat,time)
   end do

   do ie = 1, mesh%ne

      v0   = mesh%edp(ie)%v(1) ! tri point
      v1   = mesh%edp(ie)%v(2) ! tri point
      v2   = mesh%edt(ie)%v(1) ! vtx point
      v3   = mesh%edt(ie)%v(2) ! vtx point
  
      v0v1 = mesh%tri(v1)%c%p-mesh%tri(v0)%c%p
      v2v3 = mesh%vtx(v3)%p-mesh%vtx(v2)%p
     
      flag = sign(1._r8, dot_product(v0v1,mesh%edp(ie)%tg))
      normal_wind%f(ie)  = flag*(stream_function_dual_cell%f(v1)-stream_function_dual_cell%f(v0))/(rearth*mesh%edp(ie)%leng)

      flag = sign(1._r8, dot_product(v2v3,mesh%edt(ie)%tg))
      tangent_wind%f(ie) = -flag*(stream_function_prime_cell%f(v3)-stream_function_prime_cell%f(v2))/(rearth*mesh%edt(ie)%leng)

   end do

   deallocate(stream_function_prime_cell%f)
   deallocate(stream_function_dual_cell%f)
   
   return
  end subroutine overwrite_wind_streamf

!================================================
! Calculates velocity field on nodes at a time 
! step. vel must have the position already set
!================================================

  subroutine ini_vector_velocity(mesh,vector_velocity,time, &
                                      scalar_u_at_edge,     &
                                      scalar_v_at_edge)
! io
    type(global_domain),    intent(in)    :: mesh
    type(vector_type),      intent(inout) :: vector_velocity
    real(r8),               intent(in)    :: time
    type(scalar_1d_field),  intent(inout) :: scalar_u_at_edge     ! U wind at initial time
    type(scalar_1d_field),  intent(inout) :: scalar_v_at_edge     ! V ..............
! local
    integer(i4)           :: iv
    integer(i4)           :: it
    integer(i4)           :: ie
    integer(i4)           :: kk
    real (r8)             :: utmp
    real (r8)             :: vtmp

    select case(vector_velocity%pos)
    case(0) !Vectors on HX's center
       if(.not. allocated(vector_velocity%p))then
          vector_velocity%n=mesh%nv
          allocate(vector_velocity%p(1:mesh%nv_full))
       end if
       do iv=1, mesh%nv
          utmp = u_init(mesh%vtx(iv)%lon,mesh%vtx(iv)%lat, time)
          vtmp = v_init(mesh%vtx(iv)%lon,mesh%vtx(iv)%lat, time)
          call convert_vector_sph2cart(utmp,vtmp,mesh%vtx(iv)%p,vector_velocity%p(iv)%v)
       end do
    case(1) !Vectors on Triangle Circumcenters
       if(.not. allocated(vector_velocity%p))then
          vector_velocity%n=mesh%nt
          allocate(vector_velocity%p(1:mesh%nt_full))
       end if
       do it=1, mesh%nt
          utmp = u_init(mesh%tri(it)%c%lon,mesh%tri(it)%c%lat, time)
          vtmp = v_init(mesh%tri(it)%c%lon,mesh%tri(it)%c%lat, time)
          call convert_vector_sph2cart(utmp,vtmp,mesh%tri(it)%c%p,vector_velocity%p(it)%v)
       end do
    case(2) !Vectors on TR's edge (mid-point)

       if(.not. allocated(vector_velocity%p))then
          vector_velocity%n=mesh%ne
          allocate(vector_velocity%p(1:mesh%ne_full))
       end if

       do ie=1, mesh%ne
          utmp = u_init(mesh%edt(ie)%c%lon,mesh%edt(ie)%c%lat, time)
          vtmp = v_init(mesh%edt(ie)%c%lon,mesh%edt(ie)%c%lat, time)
          scalar_u_at_edge%f(ie) = utmp
          scalar_v_at_edge%f(ie) = vtmp
          call convert_vector_sph2cart(utmp,vtmp,mesh%edt(ie)%c%p,vector_velocity%p(ie)%v)
       end do

    case(6) !Vectors on HX edges

       if(.not. allocated(vector_velocity%p))then
          vector_velocity%n=mesh%ne
          allocate(vector_velocity%p(1:mesh%ne_full))
       end if

       do ie=1, mesh%ne
          utmp = u_init(mesh%edp(ie)%c%lon,mesh%edp(ie)%c%lat, time)
          vtmp = v_init(mesh%edp(ie)%c%lon,mesh%edp(ie)%c%lat, time)
          scalar_u_at_edge%f(ie) = utmp
          scalar_v_at_edge%f(ie) = vtmp
          call convert_vector_sph2cart(utmp,vtmp,mesh%edp(ie)%c%p,vector_velocity%p(ie)%v)
       end do
    case default
       print*, "calc_vel error: Please set a correct position for the velocity :", vector_velocity%pos
       stop
    end select

   return
  end subroutine ini_vector_velocity

!================================================
!                 Functions
!================================================

!================================================
! A pointwise initial of scalar field
! Lat in [-pi/2,pi/2], Lon in [-pi,pi]
!================================================
  real(r8) function ini_scalar(lon, lat, time)
! io
    real(r8), intent(in) :: lon
    real(r8), intent(in) :: lat
    real(r8), intent(in), optional :: time
! local
    !Center of bell, slot, hill
    real(r8)   :: lat1, lon1, lat2, lon2
    real(r8)   :: r1, r2, h1, h2
    real(r8)   :: b_value, c_value
    real(r8)   :: x0, y0, z0
    real(r8)   :: x1, y1, z1
    real(r8)   :: x2, y2, z2
    real(r8)   :: hmax
    real(r8)   :: base_radius
    real(r8)   :: ae, dis, u0
    ! For RH wave
    real(r8)   :: h0 
    real(r8)   :: wave_r
    real(r8)   :: omega_0
    real(r8)   :: kkk
    real(r8)   :: part_a, part_b, part_c
    !Exponential/Gaussian parameters
    real(r8)   :: alpha
    ! Shamir and Paldor, 2016
    real(r8)   :: phase_speed, uTilde, vTilde, hTilde

! initial scalar fields
    select case(initialfield)
    case(0) ! Constant
       ini_scalar=1._r8
!================================================
!           Test Cases in W1992
!================================================
    case(1) ! One Cosine Bell
       base_radius = 1._r8/3._r8
       dis         = arcdistll(lon, lat, pi/2., 0._r8)
       if(dis.lt.base_radius)then
          ini_scalar = 500._r8*(1.0_r8+cos(pi*dis/base_radius))
       else
          ini_scalar = 0.0_r8
       endif
    case(99)! One Gaussian hill
       dis         = arcdistll(lon, lat, pi/2., 0._r8)
       ini_scalar  = 1000._r8*exp(-3._r8*pi*pi*dis*dis)
    case(2) ! Global Steady State Nonlinear Zonal Geostrophic Flow 
       alpha      = 0._r8
       u0         = 2._r8*pi*rearth/(12._r8*day2sec)
       ini_scalar = (2.94e4-(rearth*omega*u0+0.5_r8*(u0**2))*&
                           (-cos(lon)*cos(lat)*sin(alpha)+sin(lat)*cos(alpha))**2)/gravity
    case(5) ! Mountain induced Rossby wave
       alpha      = 0._r8
       u0         = 20._r8
       ini_scalar = (gravity*5960._r8-(rearth*omega*u0+0.5_r8*(u0**2))*&
                           (-cos(lon)*cos(lat)*sin(alpha)+sin(lat)*cos(alpha))**2)/gravity
    case(6) ! Rossby-Haurwitz Wave

       h0     = 8000._r8
       kkk    = 7.848e-6_r8
       omega_0= 7.848e-6_r8
       wave_r = 4._r8

       part_a = 0.5_r8*omega_0*(2*omega+omega_0)*(cos(lat)**2._r8)
       part_a = part_a+0.25_r8*(kkk**2._r8)*(cos(lat)**(2._r8*wave_r))*((wave_r+1._r8)*cos(lat)**2._r8+&
               (2._r8*wave_r**2._r8-wave_r-2._r8)-2._r8*(wave_r**2._r8)*(cos(lat)**(-2._r8)))

       part_b = (2._r8*(omega+omega_0)*kkk*(cos(lat)**wave_r))/&
                ((wave_r+1)*(wave_r+2._r8))*((wave_r**2._r8+2._r8*wave_r+2._r8)-((wave_r+1)**2._r8)*(cos(lat)**2))
       part_c = 0.25_r8*(kkk**2._r8)*(cos(lat)**(wave_r*2._r8))*((wave_r+1._r8)*cos(lat)**2._r8-(wave_r+2._r8))

       ini_scalar = gravity*h0
       ini_scalar = ini_scalar+(rearth**2)*part_a
       ini_scalar = ini_scalar+(rearth**2)*part_b*cos(wave_r*lon)
       ini_scalar = ini_scalar+(rearth**2)*part_c*cos(2*wave_r*lon)
       ini_scalar = ini_scalar/gravity

!================================================
!           Test Cases in NL 2010
!================================================

    case(11)

       lon1        = pi
       lat1        = pi/3._r8
       lon2        = pi
       lat2        = -pi/3._r8
       hmax        = 500._r8
       base_radius = rearth/3._r8
       b_value     = 0._r8
       c_value     = 1._r8

       r1 = rearth*arcdistll(lon, lat, lon1, lat1)
       r2 = rearth*arcdistll(lon, lat, lon2, lat2)
       h1 = hmax*(1._r8+cos(pi*r1/base_radius))
       h2 = hmax*(1._r8+cos(pi*r2/base_radius))

       if(r1.lt.base_radius)then
          ini_scalar = b_value+c_value*h1
       else if(r2.lt.base_radius)then
          ini_scalar = b_value+c_value*h2
       else
          ini_scalar = b_value
       end if

    case(12,14)
       lon1        =-5_r8*pi/6._r8
       lat1        = 0_r8
       lon2        = 5_r8*pi/6._r8
       lat2        = 0_r8
       hmax        = 0.5_r8
       base_radius = 0.5_r8
       b_value     = 0._r8
       c_value     = 1._r8

       r1 = arcdistll(lon, lat, lon1, lat1)
       r2 = arcdistll(lon, lat, lon2, lat2)
       h1 = hmax*(1._r8+cos(pi*r1/base_radius))
       h2 = hmax*(1._r8+cos(pi*r2/base_radius))

       if(r1.lt.base_radius)then
          ini_scalar = b_value+c_value*h1
       else if(r2.lt.base_radius)then
          ini_scalar = b_value+c_value*h2
       else
          ini_scalar = b_value
       end if

    case(24) ! 24 is 14 with b c changed 

       lon1        =-5_r8*pi/6._r8
       lat1        = 0_r8
       lon2        = 5_r8*pi/6._r8
       lat2        = 0_r8
       hmax        = 0.5_r8
       base_radius = 0.5_r8
       b_value     = 0.1_r8
       c_value     = 0.9_r8

       r1 = arcdistll(lon, lat, lon1, lat1)
       r2 = arcdistll(lon, lat, lon2, lat2)
       h1 = hmax*(1._r8+cos(pi*r1/base_radius))
       h2 = hmax*(1._r8+cos(pi*r2/base_radius))

       if(r1.lt.base_radius)then
          ini_scalar = b_value+c_value*h1
       else if(r2.lt.base_radius)then
          ini_scalar = b_value+c_value*h2
       else
          ini_scalar = b_value
       end if

     case(34) ! 34 is correlated with 24 as in L2012

       lon1        =-5_r8*pi/6._r8
       lat1        = 0_r8
       lon2        = 5_r8*pi/6._r8
       lat2        = 0_r8
       hmax        = 0.5_r8
       base_radius = 0.5_r8
       b_value     = 0.1_r8
       c_value     = 0.9_r8

       r1 = arcdistll(lon, lat, lon1, lat1)
       r2 = arcdistll(lon, lat, lon2, lat2)
       h1 = hmax*(1._r8+cos(pi*r1/base_radius))
       h2 = hmax*(1._r8+cos(pi*r2/base_radius))

       if(r1.lt.base_radius)then
          ini_scalar = b_value+c_value*h1
       else if(r2.lt.base_radius)then
          ini_scalar = b_value+c_value*h2
       else
          ini_scalar = b_value
       end if
       ini_scalar = -0.8_r8*ini_scalar*ini_scalar+0.9_r8

    case(13)
       lon1        =-3_r8*pi/4._r8
       lat1        = 0_r8 
       lon2        = 3_r8*pi/4._r8
       lat2        = 0_r8
       hmax        = 0.5_r8
       base_radius = 0.5_r8
       b_value     = 0.1_r8
       c_value     = 1._r8

       r1 = arcdistll(lon, lat, lon1, lat1)
       r2 = arcdistll(lon, lat, lon2, lat2)
       h1 = hmax*(1._r8+cos(pi*r1/base_radius))
       h2 = hmax*(1._r8+cos(pi*r2/base_radius))

       if(r1.lt.base_radius)then
          ini_scalar = b_value+c_value*h1
       else if(r2.lt.base_radius)then
          ini_scalar = b_value+c_value*h2
       else
          ini_scalar = b_value
       end if
    
    case(15) ! One Slotted Cylinder

       lon1        = pi/2_r8
       lat1        = 0_r8

       base_radius = 0.5_r8
       c_value     = 1_r8
       b_value     = 0.1_r8
       r1          = arcdistll(lon, lat, lon1, lat1)

       if(r1.le.base_radius.and.abs(lon-lon1).ge.(base_radius/6_r8))then
         ini_scalar = c_value
       else if(r1.le.base_radius.and.abs(lon-lon1).lt.(base_radius/6_r8).and.&
               (lat-lat1).lt.(-5_r8*base_radius/12_r8))then
         ini_scalar = c_value
       else
         ini_scalar = b_value
       end if

    case(25) ! Two Slotted Cylinders

       lon1        =-5_r8*pi/6._r8
       lat1        = 0_r8
       lon2        = 5_r8*pi/6._r8
       lat2        = 0_r8

       base_radius = 0.5_r8
       c_value     = 1._r8
       b_value     = 0.1_r8
       r1          = arcdistll(lon, lat, lon1, lat1)
       r2          = arcdistll(lon, lat, lon2, lat2)

       if(r1.le.base_radius.and.abs(lon-lon1).ge.(base_radius/6_r8))then
         ini_scalar = c_value
       else if(r2.le.base_radius.and.abs(lon-lon2).ge.(base_radius/6_r8))then
         ini_scalar = c_value
       else if(r1.le.base_radius.and.abs(lon-lon1).lt.(base_radius/6_r8).and.&
               (lat-lat1).lt.(-5_r8*base_radius/12_r8))then
         ini_scalar = c_value
       else if(r2.le.base_radius.and.abs(lon-lon2).lt.(base_radius/6_r8).and.&
               (lat-lat2).gt.(5_r8*base_radius/12_r8))then
         ini_scalar = c_value
       else
         ini_scalar = b_value
       end if

    case(20,30) ! Gassian Hills

       lon1        = -5_r8*pi/6._r8
       lat1        = 0._r8
       lon2        = 5._r8*pi/6._r8
       lat2        = 0._r8

       hmax        = 1._r8
       b_value     = 5._r8

       x0   = cos(lat)*cos(lon)
       y0   = cos(lat)*sin(lon)
       z0   = sin(lat)

       x1   = cos(lat1)*cos(lon1)
       y1   = cos(lat1)*sin(lon1)
       z1   = sin(lat1)

       x2   = cos(lat2)*cos(lon2)
       y2   = cos(lat2)*sin(lon2)
       z2   = sin(lat2)

       h1   = hmax*exp(-b_value*((x0-x1)**2+(y0-y1)**2+(z0-z1)**2))
       h2   = hmax*exp(-b_value*((x0-x2)**2+(y0-y2)**2+(z0-z2)**2))
       ini_scalar = h1+h2

    case(22)  ! two Gaussian hills as from case99
      
       lon1        = -5_r8*pi/6._r8
       lat1        = 0._r8
       lon2        = 5._r8*pi/6._r8
       lat2        = 0._r8

       dis         = arcdistll(lon, lat, lon1,lat1)
       h1          = 1._r8*exp(-5._r8*pi*pi*dis*dis)
       dis         = arcdistll(lon, lat, lon2,lat2)
       h2          = 1._r8*exp(-5._r8*pi*pi*dis*dis)
       ini_scalar  = h1+h2

    case(90) ! Colliding Modons
       ini_scalar = 10000._r8
    case(51) ! SP16, Rossby
       call getPhaseSpeed(phase_speed,0)
       call getAmplitudes(lat,uTilde,vTilde,hTilde,0)
       ini_scalar = 5.0e4_r8/gravity+hTilde*cos(10._r8*lon)
    case(52) ! SP16, WIG
       call getPhaseSpeed(phase_speed,-1)
       call getAmplitudes(lat,uTilde,vTilde,hTilde,-1)
       ini_scalar = 5.0e4_r8/gravity+hTilde*cos(10._r8*lon)
    case(53) ! SP16, EIG
       call getPhaseSpeed(phase_speed,1)
       call getAmplitudes(lat,uTilde,vTilde,hTilde,1)
       ini_scalar = 5.0e4_r8/gravity+hTilde*cos(10._r8*lon)
    case default
       ini_scalar=0.
    end select

    return
  end function ini_scalar

  real(r8) function ini_topo(lon, lat)

!================================================
!  F - initial conditions for scalar fields
!   Lat in [-pi/2,pi/2], Lon in [-pi,pi]
!================================================
! io
    real(r8), intent(in) :: lon
    real(r8), intent(in) :: lat
! local
    real(r8) :: rrr
    real(r8) :: lat1, lon1, lat2, lon2
    real(r8) :: lonp, latp
    real(r8) :: u0

       lon1=0.
       lat1=0.
       lon2=0.
       lat2=0.
    !Initial scalar fields
    select case(initialfield)
    case(5) ! Mountain induced Rossby wave
      rrr = min((pi/9._r8)**2,((lon-pi/2._r8)**2+(lat-pi/6._r8)**2))
      rrr = sqrt(rrr)
      ini_topo=2000._r8*(1._r8-rrr*9._r8/pi)
    case(33) ! PTB's small depth 
      rrr = min((pi/9._r8)**2,((lon-pi/2._r8)**2+(lat-pi/6._r8)**2))
      rrr = sqrt(rrr)
      ini_topo=2000._r8*(1._r8-rrr*9._r8/pi)
    case default
      ini_topo=0.
    end select

    return
  end function ini_topo

  function velocity(p, time)
!================================================
!  V - velocity in Cartesian coordinates
!  p in Cartesian coords
!================================================
    real (r8), intent(in) :: p(1:3)
    real (r8):: velocity(1:3)

    !Period, actual time
    real (r8), intent(in)::  time

    !Auxiliar variables
    real (r8):: lon
    real (r8):: lat
    real (r8):: utmp
    real (r8):: vtmp

    call cart2sph(p(1), p(2), p(3), lon, lat)
    utmp=u_init(lon, lat, time)
    vtmp=v_init(lon, lat, time)
    call convert_vector_sph2cart(utmp, vtmp, p, velocity)
    return
  end function velocity

!================================================
!  U - velocity in West-East direction
!   for a given testcase
!
!   Lat in [-pi/2,pi/2], Lon in [-pi,pi]
!================================================
  real(r8) function u_init(lon, lat, time)
! io
    real(r8), intent(in) :: lon
    real(r8), intent(in) :: lat
    real(r8), intent(in) :: time
! local
    real(r8)   :: omega_m
    real(r8)   :: omega_0
    !real(r8)   :: a
    real(r8)   :: duration
    !Period, actual time
    !Auxiliar variables
    real(r8)   :: k
    real(r8)   :: lonp
    real(r8)   :: ae 
    real(r8)   :: alpha
    real(r8)   :: u0
    real(r8)   :: kkk
    real(r8)   :: wave_r
    integer(i4):: wave_m
!
! colliding modons
!
    real(r8)   :: center_a_lat, center_b_lat
    real(r8)   :: center_a_lon, center_b_lon
    real(r8)   :: modon_a, modon_b
    real(r8)   :: dist_a,  dist_b
! Shamir and Paldor, 2016
    real(r8)   :: phase_speed, uTilde, vTilde, hTilde

    u_init = 0._r8
    ! Velocity for each testcase
    select case(testcase)

!================================================
!      Test Cases in Williamson et al. (1992)
!================================================
    case(1,15) ! Solid body rotation
       alpha  = 0._r8
       u0     = 2*pi*rearth/(12._r8*day2sec)
       u_init = u0*(cos(lat)*cos(alpha)+sin(lat)*cos(lon)*sin(alpha))
    case(2)    ! Global Steady State Nonlinear Zonal Geostrophic Flow
       alpha  = 0._r8 
       u0     = 2._r8*pi*rearth/(12._r8*day2sec)
       u_init = u0*(cos(lat)*cos(alpha)+cos(lon)*sin(lat)*sin(alpha))
    case(5)    ! Mountain induced Rossby wave
       alpha  = 0._r8 
       u0     = 20._r8 
       u_init = u0*(cos(lat)*cos(alpha)+cos(lon)*sin(lat)*sin(alpha))
    case(6)    ! Rossby-Haurwitz Wave
       omega_0= 7.848e-6_r8
       kkk    = 7.848e-6_r8
       wave_r = 4._r8
       u_init = rearth*omega_0*cos(lat)+&
                rearth*kkk*((cos(lat))**(wave_r-1))*&
               (wave_r*(sin(lat)**2)-(cos(lat)**2))*&
               cos(wave_r*lon)
!================================================
!             Test Cases in NL(2010)
!================================================
    case(11)    ! TC1
       duration = 12._r8*day2sec
       k        = 10._r8*rearth/duration
       u_init   = k*(sin((lon)/2._r8)**2)*(sin(2._r8*lat))*(cos(pi*time/duration))
    case(12)    ! TC2
       duration = 12._r8*day2sec
       k        = 10._r8*rearth/duration
       u_init   = k*(sin((lon))**2)*(sin(2_r8*lat))*(cos(pi*time/duration))
    case(13)    ! TC3
       duration = 12._r8*day2sec
       k        = 5._r8*rearth/duration
       u_init   = -k*(sin((lon)*0.5_r8)**2)*(sin(2.*lat))*(cos(lat)**2)*(cos(pi*time/duration))
    case(14,20) ! TC4
       duration = 12._r8*day2sec
       k        = 10._r8*rearth/duration
       lonp     = lon-(2._r8*pi*time/duration)
       u_init   = k*((sin(lonp))**2)*(sin(2_r8*lat))*(cos(pi*time/duration))+(2._r8*pi*rearth*cos(lat))/duration
    case(30)    ! case1+case4 as in sg11?
       duration = 12._r8*day2sec
       k        = 10._r8*rearth/duration
       lonp     = lon-(2._r8*pi*time/duration)
       u_init   = k*((sin(lonp*0.5_r8))**2)*(sin(2_r8*lat))*(cos(pi*time/duration))+(2_r8*pi*rearth*cos(lat))/duration
    case(90)    ! Colliding Modons

       center_a_lat = 0._r8
       center_b_lat = 0._r8
       center_a_lon = pi*0.5_r8
       center_b_lon = pi*1.5_r8
       dist_a       = rearth*arcdistll(lon,lat,center_a_lon,center_a_lat)
       dist_b       = rearth*arcdistll(lon,lat,center_b_lon,center_b_lat)
       modon_a      = 40._r8*exp(-(dist_a/modon_radius)**2)
       modon_b      = 40._r8*exp(-(dist_b/modon_radius)**2)
       u_init       = modon_a-modon_b

    case(51) ! SP16, Rossby
       call getPhaseSpeed(phase_speed,0)
       call getAmplitudes(lat,uTilde,vTilde,hTilde,0)
       u_init = uTilde*cos(10._r8*lon)
    case(52) ! SP16, WIG
       call getPhaseSpeed(phase_speed,-1)
       call getAmplitudes(lat,uTilde,vTilde,hTilde,-1)
       u_init = uTilde*cos(10._r8*lon)
    case(53) ! SP16, EIG
       call getPhaseSpeed(phase_speed,1)
       call getAmplitudes(lat,uTilde,vTilde,hTilde,1)
       u_init = uTilde*cos(10._r8*lon)
    case default
       print*, "Error in transport tests: Unknown vector field ", testcase
       stop
    end select

    return
  end function u_init

!================================================
!  V - velocity in South-North direction
!   for a given testcase
!   Lat in [-pi/2,pi/2], Lon in [-pi,pi]
!================================================

  real(r8) function v_init(lon, lat, time)
! io
    real(r8), intent(in) :: lon
    real(r8), intent(in) :: lat
    real(r8), intent(in) :: time
! local
    integer(i4):: wave_m
    real(r8)   :: omega_m
    real(r8)   :: omega_0
    real(r8)   :: k
    real(r8)   :: lonp
    real(r8)   :: ae
    real(r8)   :: alpha
    real(r8)   :: u0
    real(r8)   :: kkk
    real(r8)   :: wave_r
    real(r8)   :: duration 
! Shamir and Paldor, 2016
    real(r8)   :: phase_speed, uTilde, vTilde, hTilde


    ! Velocity for each testcase
    v_init = 0._r8
    select case(testcase)
!================================================
!      Test Cases in Williamson et al. (1992)
!================================================

    case(1,15) ! Solid body rotation
       alpha  = 0._r8
       u0     = 2*pi*rearth/(12._r8*day2sec)
       v_init = -u0*sin(lon)*sin(alpha)
    case(2)    ! Global Steady State Nonlinear Zonal Geostrophic Flow
       alpha  = 0._r8
       u0     = 2._r8*pi*rearth/(12*day2sec)
       v_init = -u0*sin(lon)*sin(alpha)
    case(5)    ! Mountain induced Rossby Wave
       alpha  = 0._r8
       u0     = 20._r8
       v_init = -u0*sin(lon)*sin(alpha)
    case(6)    ! Rossby-Haurwitz Wave
       kkk    = 7.848e-6_r8
       wave_r = 4._r8
       v_init = -rearth*kkk*wave_r*((cos(lat))**(wave_r-1))*&
                 sin(lat)*sin(wave_r*lon)

!================================================
!            Test Cases in NL (2010)
!================================================

    case(11)    ! TC1
       duration = 12._r8*day2sec
       k        = 10._r8*rearth/duration
       v_init   = 0.5_r8*k*(sin(lon))*(cos(lat))*(cos(pi*time/(12._r8*day2sec)))
    case(12)    ! TC2
       duration = 12._r8*day2sec
       k        = 10._r8*rearth/duration
       v_init   = k*(sin(2*(lon)))*(cos(lat))*(cos(pi*time/duration)) 
    case(13)    ! TC3
       duration = 12._r8*day2sec
       k        = 5._r8*rearth/duration
       v_init   = 0.5_r8*k*(sin((lon)))*(cos(lat)**3)*(cos(pi*time/duration))
    case(14,20) ! TC4
       duration = 12._r8*day2sec
       k        = 10._r8*rearth/duration
       lonp     = lon-(2._r8*pi*time/duration)
       v_init   = k*(sin(2*lonp))*(cos(lat))*(cos(pi*time/duration))
    case(30)    ! case1+case4 as in sg11?
       duration = 12._r8*day2sec
       k        = 10._r8*rearth/duration
       lonp     = lon-(2._r8*pi*time/duration)
       v_init   = 0.5_r8*k*(sin(lonp))*(cos(lat))*(cos(pi*time/duration))
    case(90)    ! Colliding Modons
       v_init   = 0._r8
    case(51) ! SP16, Rossby
       call getPhaseSpeed(phase_speed,0)
       call getAmplitudes(lat,uTilde,vTilde,hTilde,0)
       v_init = vTilde*cos(10._r8*lon- 0.5_r8*pi)
    case(52) ! SP16, WIG
       call getPhaseSpeed(phase_speed,-1)
       call getAmplitudes(lat,uTilde,vTilde,hTilde,-1)
       v_init = vTilde*cos(10._r8*lon- 0.5_r8*pi)
    case(53) ! SP16, EIG
       call getPhaseSpeed(phase_speed,1)
       call getAmplitudes(lat,uTilde,vTilde,hTilde,1)
       v_init = vTilde*cos(10._r8*lon- 0.5_r8*pi)
    case default
       print*, "Error in transport tests: Unknown vector field ", testcase
       stop
    end select
    return
  end function v_init

  function ini_stream_function(lon, lat, time)

!================================================
!  Initial Conditions for scalar fields
!  Lat in [-pi/2,pi/2], Lon in [-pi,pi]
!================================================
! io
    real(r8), intent(in) :: lon
    real(r8), intent(in) :: lat
    real(r8), intent(in), optional :: time
    real(r8)             :: ini_stream_function
! local
    real(r8)             :: duration, k, lonp
    real(r8)             :: ooo, kkk, rrr
    

    select case(testcase)
    case(1,15) ! solid body rotation flow
        ini_stream_function = -2._r8*pi*rearth*rearth*sin(lat)/(12._r8*day2sec)
    case(6)    ! RH wave
        ooo     = 7.848e-6_r8
        kkk     = 7.848e-6_r8
        rrr     = 4._r8
        ini_stream_function = (rearth**2)*(-ooo*sin(lat)+kkk*(cos(lat)**rrr)*sin(lat)*cos(rrr*lon))
    case(14,20)! used by deformatinoal flow of two cosine bells and gaussian hills 
        duration = 12._r8*day2sec
        k        = 10._r8*rearth/duration
        lonp     = lon-(2._r8*pi*time/duration)
        ini_stream_function = rearth*(k*((sin(lonp))**2)*(cos(lat)**2)*(cos(pi*time/duration))-&
                              (2._r8*pi*rearth*sin(lat))/duration)
    case default
        ini_stream_function = 0._r8
    end select

    return
  end function ini_stream_function

!========================================================================
! Below subroutines are used for building test case in 
! Shamir and Paldor 2016
! phase speed
!========================================================================

  subroutine getPhaseSpeed(C,waveFlag)

    implicit none

    integer,          intent(in)    :: waveFlag    ! -1=WIG, 0=Rossby, 1=EIG
    real(kind=r8),    intent(inout) :: C           ! phase speed (rad/sec)

    real(kind=r8)                   :: Cj(1:3)     ! see Eq. 7  in text
    real(kind=r8)                   :: Delta0      ! see Eq. 8a in text
    complex(kind=r8)                :: Deltaj      ! see eq. 8b in text
    real(kind=r8)                   :: Delta4      ! see Eq. 8c in text
    real(kind=r8)                   :: En          ! see Eq. 9  in text

    complex(kind=r8), parameter     :: i1 = cmplx(0.0_r8,1.0_r8,r8)
    complex(kind=r8), parameter     :: r1 = cmplx(1.0_r8,0.0_r8,r8)
    real(kind=r8),    parameter     :: r2 = 0.5_r8
    real(kind=r8),    parameter     :: r3 = 1.0_r8/3.0_r8
    integer,          parameter     :: n  = 5                        ! chosen mode number
    integer,          parameter     :: k  = 10                       ! chosen wave number

    integer                         :: j
    real(kind=r8),    parameter     :: H0 = 5.0e4_r8/gravity                 ! Layer's mean depth (m)
    real(kind=r8),    parameter  :: sigma = 0.5_r8 + ( 0.25_r8 + k**2 ) ** 0.5_r8           ! see Eq.10 in text


    En     =  gravity*H0 / rearth**2 * ( n + sigma )**2
    Delta0 =  03.0_r8 * k**2 * En
    Delta4 = -54.0_r8 * k**4 * gravity*H0 * omega / rearth**2

    do j=1,3
       Deltaj = ( r1 * ( Delta4**2 - 4.0_r8 * Delta0**3 ) )**r2
       Deltaj = ( r2 * ( Delta4    +          Deltaj    ) )**r3
       Deltaj = Deltaj * exp(2.0_r8*pi * i1 * j * r3)

       Cj(j)  = real( -r3 / k**2 * ( Deltaj + Delta0 / Deltaj ) )
    end do

    select case (waveFlag)
      case(0)
        C = -minval(abs(Cj))
      case(1)
        C = maxval(Cj)
      case(-1)
        C = minval(Cj)
    end select

  end subroutine getPhaseSpeed

!========================================================================
! eigenfunction
!========================================================================
  subroutine getPsi(lat,psi,dpsi)

    implicit none

    real(kind=r8), intent(in)    :: lat            ! latitude (rad)
    real(kind=r8), intent(inout) :: psi            ! see Eq. 14 in text
    real(kind=r8), intent(inout) :: dpsi           ! meridional derivative of psi

    integer,          parameter  :: n     = 5                        ! chosen mode number
    integer,          parameter  :: k     = 10                       ! chosen wave number
    real(kind=r8), parameter     :: amp = 1.0e-8_r8   ! arbitrary amplitude for linear waves
    real(kind=r8)                :: C5     ! see Eq. 19 in text
    real(kind=r8)                :: C5p    ! meridional derivative of C5
    real(kind=r8), parameter     :: sigma = 0.5_r8 + ( 0.25_r8 + k**2 ) ** 0.5_r8           ! see Eq.10 in text

    real(kind=r8), parameter     :: a3 = sigma * ( sigma + 1 ) * ( sigma + 2 )
    real(kind=r8), parameter     :: a4 = a3    * ( sigma + 3 )
    real(kind=r8), parameter     :: a5 = a4    * ( sigma + 4 )

    C5    = ( 04.0_r8 * a5  * sin(lat)**4  - &
              20.0_r8 * a4  * sin(lat)**2  + &
              15.0_r8 * a3) * sin(lat)     / 15.0_r8

    C5p   = ( 04.0_r8 * a5  * sin(lat)**4  - &
              12.0_r8 * a4  * sin(lat)**2  + &
              03.0_r8 * a3) * cos(lat)     / 3.0_r8

    psi   = amp * cos(lat) ** (sigma) * C5

    dpsi  = amp * cos(lat) ** (sigma) * ( - sigma * tan(lat) * C5 + C5p )

  end subroutine getPsi

!========================================================================
! latitude-dependent amplitudes
!========================================================================
  subroutine getAmplitudes(lat,uTilde,vTilde,hTilde,waveFlag)

    implicit none

    integer,       intent(in)    :: waveFlag       ! -1=WIG, 0=Rossby, 1=EIG
    real(kind=r8), intent(in)    :: lat            ! latitude (rad)
    real(kind=r8), intent(inout) :: uTilde         ! see Eq. 18 in text
    real(kind=r8), intent(inout) :: vTilde         ! see Eqs. 12a and 15b in text
    real(kind=r8), intent(inout) :: hTilde         ! see Eqs. 12b and 15a in text

    real(kind=r8)                :: Kp     ! see Eq. 13 in text
    real(kind=r8)                :: Km     ! see Eq. 13 in text
    real(kind=r8)                :: W1     ! see Eq. 16 in text
    real(kind=r8)                :: W2     ! see Eq. 17 in text
    real(kind=r8)                :: C                 ! phase speed (rad/sec)
    real(kind=r8)                :: psi    ! see Eq. 13 in text
    real(kind=r8)                :: dpsi   ! meridional derivative of psi

    real(kind=r8),    parameter  :: r2    = 0.5_r8
    real(kind=r8),    parameter  :: o2    = 2.0_r8 * omega
    integer,          parameter  :: n     = 5                        ! chosen mode number
    integer,          parameter  :: k     = 10                       ! chosen wave number
    real(kind=r8),    parameter  :: sigma = 0.5_r8 + ( 0.25_r8 + k**2 ) ** 0.5_r8           ! see Eq.10 in text
    real(kind=r8),    parameter  :: H0 = 5.0e4_r8/gravity                 ! Layer's mean depth (m)

    call getPhaseSpeed(C,waveFlag)
    call getPsi(lat,psi,dpsi)

    select case (waveFlag)
    case(0)
       Kp = (gravity*H0 + rearth**2 * C**2 * cos(lat)**2) / (C*cos(lat))
       Km = (gravity*H0 - rearth**2 * C**2 * cos(lat)**2) / (C*cos(lat))

       vTilde = ( o2    * abs(Km)  / cos(lat)**2  )**r2 * psi
       hTilde = ( o2    * abs(Km)  * rearth**2 * H0**2 )**r2 / Km * &
                ( dpsi  + tan(lat) * ( r2 * Kp/Km - o2 / C)  * psi )
       uTilde = ( o2    * sin(lat) / C ) * vTilde + &
                ( gravity / rearth / cos(lat) / C ) * hTilde

    case(-1,1)
       W1 = ( (k*C)**2 - o2**2 * sin(lat)**2 ) / ( C  * cos(lat) )
       W2 = ( (k*C)**2 - o2**2 * sin(lat)**2 ) / ( o2 * cos(lat)**2 )

       vTilde = ( o2    * abs(W1)  * gravity*H0 / cos(lat)**2 )**r2 / W1 * &
                ( dpsi  + tan(lat) * ( r2 - o2 / W2 + o2 / C ) * psi )
       hTilde = ( o2    * abs(W1)  * rearth**2 * H0 / gravity )**r2 * psi
       uTilde = ( o2    * sin(lat) / C ) * vTilde + &
                ( gravity / rearth / cos(lat) / C ) * hTilde
    end select

    vTilde = k * vTilde

  end subroutine getAmplitudes

  end module swe_init_module
