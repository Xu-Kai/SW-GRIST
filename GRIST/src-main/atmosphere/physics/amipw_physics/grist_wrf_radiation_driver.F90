
!----------------------------------------------------------------------------
! Created on 2019
! Author: Yi Zhang
! Version 1.0
! Description: This driver mainly follows the naming convention of WRF, but 
!              the dim conversion requires a transpose from GRIST based vars, 
!              currently only support RRTMG-SW/LW as default
! Revision history:
!                 1. RRTMG SW&LW are default now
!                 2. RRTM LW is ok, CAM-SW is ok
!--------------------------------------------------------------------------

module grist_wrf_radiation_driver

  use grist_hpe_constants,     only: eta_full, eta_full_a, eta_face
  use grist_mpi
  use grist_constants,         only: deg2rad, i4, r8, r_d=>rdry, r_v=>rvap, cp, cv, gravity, one, zero
  use grist_wrf_data_structure,only: pstate_wrf, ptend_wrf, psurf_wrf, restart, &
                                     warm_rain
  use grist_wrf_data_structure,only: p_qv, p_qc, p_qr, p_qi, p_qs, p_qg, param_first_scalar
  use grist_nml_module,        only: start_ymd, start_tod, doAquaPlanet
  use grist_wrfphys_nml_module,only: wrfphys_cf_scheme, wrfphys_rasw_scheme, wrfphys_ralw_scheme, step_ra, wphys_has_req
  use grist_time_manager,      only: get_curr_calday
  use grist_handle_error,      only: endrun
!
! WRF module RA
!
  use module_wrfmodel_constants,only: xls, xlv, xlf, rhowater, rhosnow, ep_2, svp1, svp2, svp3, svpt0, stbolt 
  use module_ra_sw_v2,          only: swinit, swrad          ! mm5
  use module_ra_rrtm_v2,        only: rrtminit, rrtmlwrad    ! rrtm lw
  use module_ra_cam,            only: camradinit, camrad     ! camrad
#ifdef RRTMG_V341
  use module_ra_rrtmg_sw_v341,  only: rrtmg_swrad, rrtmg_swinit
  use module_ra_rrtmg_lw_v341,  only: rrtmg_lwrad, rrtmg_lwinit
#endif
#ifdef RRTMG_V381
  use module_ra_rrtmg_sw_v381,  only: rrtmg_swrad, rrtmg_swinit
  use module_ra_rrtmg_lw_v381,  only: rrtmg_lwrad, rrtmg_lwinit
#endif
  use cam3_cloud_fraction,      only: cam3_cldfrc_init, cam3_cldfrc

  implicit  none

  private
  public  ::  grist_wrf_radiation_init, &
              grist_wrf_radiation_run,  &
              grist_wrf_radiation_final

    integer(i4)      :: icloud   = 1
    real(r8)         :: dpd      = 360._r8/365.2422_r8
    integer(i4)      :: stepra
    integer(i4)      :: stepabs
    logical, public  :: call_radiation
!
! for CAM rad
! 
    integer(i4), parameter    :: cam_abs_dim1 = 4
    integer(i4)               :: cam_abs_dim2
    integer(i4), parameter    :: levsiz       = 59
    integer(i4), parameter    :: num_months   = 12
    integer(i4), parameter    :: paerlev      = 29
    integer(i4), parameter    :: naer_c       = 13
    integer(i4), parameter    :: n_cldadv     = 3 ! for CAM, this is for all tracers, but only set to 3 now
    real(r8)                  :: pin(levsiz)
    real(r8)                  :: m_hybi(paerlev)
    real(r8), allocatable     :: m_psp(:, :)
    real(r8), allocatable     :: m_psn(:, :)
    real(r8), allocatable     :: ozmixm   (:,:,:,:)
    real(r8), allocatable     :: aerosolcp(:,:,:,:)
    real(r8), allocatable     :: aerosolcn(:,:,:,:)
    real(r8), allocatable     :: abstot_3d(:,:,:,:)
    real(r8), allocatable     :: absnxt_3d(:,:,:,:)
    real(r8), allocatable     :: emstot_3d(:,:,:)

contains

    subroutine grist_wrf_radiation_init(ncell, nLevel)
      integer(i4),  intent(in) :: ncell, nLevel
      real(r8) :: hypm(nLevel)

      hypm(1:nLevel-1) = eta_full(nLevel-1:1:-1)*1e5_r8

      call_radiation = .true. ! init radiation as .true.

      if(.not.allocated(ptend_wrf%rthraten))    allocate(ptend_wrf%rthraten   (1:ncell, 1:nLevel, 1:1));ptend_wrf%rthraten=zero
      if(.not.allocated(ptend_wrf%rthraten_lw)) allocate(ptend_wrf%rthraten_lw(1:ncell, 1:nLevel, 1:1));ptend_wrf%rthraten_lw=zero
      if(.not.allocated(ptend_wrf%rthraten_sw)) allocate(ptend_wrf%rthraten_sw(1:ncell, 1:nLevel, 1:1));ptend_wrf%rthraten_sw=zero

!----------------------------------------------------------------------
!            LWRAD 
!----------------------------------------------------------------------

      select case(trim(wrfphys_ralw_scheme))
      case('RRTM')
      call rrtminit(ptend_wrf%rthraten(   1:ncell, 1:nLevel, 1:1),   &
                    ptend_wrf%rthraten_lw(1:ncell, 1:nLevel, 1:1),   &
                    pstate_wrf%cldfra(    1:ncell, 1:nLevel, 1:1),   &
                    restart,                 &
                    1, ncell, 1, 1, 1, nLevel, &
                    1, ncell, 1, 1, 1, nLevel, &
                    1, ncell, 1, 1, 1, nLevel  )
      case('CAM')
        ! do nothing until sw init
      case('RRTMGV341','RRTMGV381')
        call rrtmg_lwinit(eta_full_a(1)*1e5, .true., &
                          1, ncell, 1, 1, 1, nLevel, &
                          1, ncell, 1, 1, 1, nLevel, &
                          1, ncell, 1, 1, 1, nLevel-1  ) 
      case default
        if(mpi_rank().eq.0) print*,"you must select a radiation LW scheme in WRFphys: wrfphys_ralw_scheme"
        call endrun
      end select

!----------------------------------------------------------------------
!            SWRAD 
!----------------------------------------------------------------------

      select case(trim(wrfphys_rasw_scheme))
      case('MM5')
      call swinit(ptend_wrf%rthraten(   1:ncell, 1:nLevel, 1:1),   &
                  ptend_wrf%rthraten_sw(1:ncell, 1:nLevel, 1:1),   &
                  restart,                 &
                  1, ncell, 1, 1, 1, nLevel, &
                  1, ncell, 1, 1, 1, nLevel, &
                  1, ncell, 1, 1, 1, nLevel  )
      case('CAM')

      if(.not.allocated(m_psp))     allocate(m_psp(    1:ncell,          1:1))
      if(.not.allocated(m_psn))     allocate(m_psn(    1:ncell,          1:1))
      if(.not.allocated(ozmixm))    allocate(ozmixm(   1:ncell, levsiz,  1:1, num_months ))
      if(.not.allocated(aerosolcp)) allocate(aerosolcp(1:ncell, paerlev, 1:1, naer_c ))
      if(.not.allocated(aerosolcn)) allocate(aerosolcn(1:ncell, paerlev, 1:1, naer_c ))
      if(.not.allocated(abstot_3d)) allocate(abstot_3d(1:ncell, 1:nLevel,1:cam_abs_dim2 ,1:1 ));abstot_3d=zero
      if(.not.allocated(absnxt_3d)) allocate(absnxt_3d(1:ncell, 1:nLevel,1:cam_abs_dim1 ,1:1 ));absnxt_3d=zero
      if(.not.allocated(emstot_3d)) allocate(emstot_3d(1:ncell, 1:nLevel,                1:1 ));emstot_3d=zero
      cam_abs_dim2 = nLevel

      call camradinit(r_d,r_v,cp,gravity,stbolt,ep_2,hypm,eta_face(1)*1e5_r8, hypm, &
                      ozmixm,    & ! real out and interpolated o3 mix ratio
                      pin,       & ! out, read from ozone_plev.formatted
                      levsiz,    & ! 59
                      pstate_wrf%xlat(1:ncell, 1:1), & ! in
                      num_months,& ! nmonth, in
                      m_psp,     & ! out
                      m_psn,     & ! out
                      m_hybi,    & ! out
                      aerosolcp, & ! out
                      aerosolcn, & ! out
                      paerlev,   & ! 29
                      naer_c,&     ! 13
                      1, ncell, 1, 1, 1, nLevel,&
                      1, ncell, 1, 1, 1, nLevel,&
                      1, ncell, 1, 1, 1, nLevel-1)
       aerosolcp = 1e-12
       aerosolcn = 1e-12
       ozmixm    = 1e-12

      case('RRTMGV341','RRTMGV381')
       call rrtmg_swinit(.true. ,                    &
                          1, ncell, 1, 1, 1, nLevel, &
                          1, ncell, 1, 1, 1, nLevel, &
                          1, ncell, 1, 1, 1, nLevel-1) 
      case default
        if(mpi_rank().eq.0) print*,"you must select a radiation LW scheme in WRFphys: wrfphys_rasw_scheme"
        call endrun
      end select
!
! cam3 cloud fraction para init
!
      call cam3_cldfrc_init(nLevel-1, eta_full(1:nLevel-1)*1e5_r8)

      return
    end subroutine grist_wrf_radiation_init

    subroutine grist_wrf_radiation_final
      if(allocated(ptend_wrf%rthraten))    deallocate(ptend_wrf%rthraten   )
      if(allocated(ptend_wrf%rthraten_lw)) deallocate(ptend_wrf%rthraten_lw)
      if(allocated(ptend_wrf%rthraten_sw)) deallocate(ptend_wrf%rthraten_sw)

      if(allocated(m_psp))     deallocate(m_psp)
      if(allocated(m_psn))     deallocate(m_psn)
      if(allocated(ozmixm))    deallocate(ozmixm)
      if(allocated(aerosolcp)) deallocate(aerosolcp)
      if(allocated(aerosolcn)) deallocate(aerosolcn)
      if(allocated(abstot_3d)) deallocate(abstot_3d)
      if(allocated(absnxt_3d)) deallocate(absnxt_3d)
      if(allocated(emstot_3d)) deallocate(emstot_3d)

     return
    end subroutine grist_wrf_radiation_final

    subroutine grist_wrf_radiation_run(ncell,nLevel,nspecies,itimestep,dtime,coszrs)

!======================================================================
! grid structure in physics part of wrf
!----------------------------------------------------------------------
! the horizontal velocities used in the physics are unstaggered
! relative to temperature/moisture variables. all predicted
! variables are carried at half levels except w, which is at full
! levels. some arrays with names (*8w) are at w (full) levels.
!
!----------------------------------------------------------------------
! in wrf, kms (smallest number) is the bottom level and kme (largest
! number) is the top level.  in your scheme, if 1 is at the top level,
! then you have to reverse the order in the k direction.
!
!         kme      -   half level (no data at this level)
!         kme    ----- full level
!         kme-1    -   half level
!         kme-1  ----- full level
!         .
!         .
!         .
!         kms+2    -   half level
!         kms+2  ----- full level
!         kms+1    -   half level
!         kms+1  ----- full level
!         kms      -   half level
!         kms    ----- full level
!
!======================================================================
! grid structure in physics part of wrf
! 
!-------------------------------------
! the horizontal velocities used in the physics are unstaggered 
! relative to temperature/moisture variables. all predicted 
! variables are carried at half levels except w, which is at full 
! levels. some arrays with names (*8w) are at w (full) levels.
!
!==================================================================
! definitions
!-----------
! theta      potential temperature (k)
! qv         water vapor mixing ratio (kg/kg)
! qc         cloud water mixing ratio (kg/kg)
! qr         rain water mixing ratio (kg/kg)
! qi         cloud ice mixing ratio (kg/kg)
! qs         snow mixing ratio (kg/kg)
!-----------------------------------------------------------------
!-- rthraten      theta tendency 
!                 due to radiation (k/s)
!-- rthratenlw    theta tendency 
!                 due to long wave radiation (k/s)
!-- rthratensw    theta temperature tendency 
!                 due to short wave radiation (k/s)
!-- dt            time step (s)
!-- itimestep     number of time steps
!-- glw           downward long wave flux at ground surface (w/m^2)
!-- gsw           net short wave flux at ground surface (w/m^2)
!-- swdown        downward short wave flux at ground surface (w/m^2)
!-- xlat          latitude, south is negative (degree)
!-- xlong         longitude, west is negative (degree)
!-- albedo                albedo (between 0 and 1)
!-- cldfra        cloud fraction (between 0 and 1)
!-- emiss         surface emissivity (between 0 and 1)
!-- rho_phy       density (kg/m^3)
!-- rr            dry air density (kg/m^3)
!-- moist         moisture array (4d - last index is species) (kg/kg)
!-- nspecies       number of moisture species
!-- p8w           pressure at full levels (pa)
!-- p_phy         pressure (pa)
!-- pb            base-state pressure (pa)
!-- pi_phy        exner function (dimensionless)
!-- dz8w          dz between full levels (m)
!-- t_phy         temperature (k)
!-- t8w           temperature at full levels (k)
!-- gmt           greenwich mean time hour of model start (hour)
!-- julday        the initial day (julian day)
!-- config_flags  boundary condition flag
!-- radt          time for calling radiation (min)
!-- degrad        conversion factor for 
!                 degrees to radians (pi/180.) (rad/deg)
!-- dpd           degrees per day for earth's 
!                 orbital position (deg/day)
!-- r_d           gas constant for dry air (j/kg/k)
!-- cp            heat capacity at constant pressure for dry air (j/kg/k)
!-- g             acceleration due to gravity (m/s^2)
!-- rvovrd        r_v divided by r_d (dimensionless)
!-- xtime         time since simulation start (min)
!-- declin        solar declination angle (rad)
!-- solcon        solar constant (w/m^2)
!-- p_qv          species index for water vapor
!-- p_qc          species index for cloud water
!-- p_qr          species index for rain water
!-- p_qi          species index for cloud ice
!-- p_qs          species index for snow
!-- p_qg          species index for graupel
!-- ids           start index for i in domain
!-- ide           end index for i in domain
!-- jds           start index for j in domain
!-- jde           end index for j in domain
!-- kds           start index for k in domain
!-- kde           end index for k in domain
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-- i_start       start indices for i in tile
!-- i_end         end indices for i in tile
!-- j_start       start indices for j in tile
!-- j_end         end indices for j in tile
!-- kts           start index for k in tile
!-- kte           end index for k in tile
!-- num_tiles     number of tiles
!
!==================================================================

! io
   integer(i4), intent(in)    :: ncell
   integer(i4), intent(in)    :: nLevel
   integer(i4), intent(in)    :: nspecies
   integer(i4), intent(in)    :: itimestep
   real(r8),    intent(in)    :: dtime
   real(r8),    intent(in)    :: coszrs(ncell) 
!
   integer(i4)  :: julday
   real(r8)     :: gmt
   real(r8)     :: radt
   real(r8)     :: calday
   real(r8)     :: julian
! local  var
!   real(r8)     :: swdown(1:ncell, 1:1)
   real(r8)     :: xtime, declin, solcon
   real(r8)     :: obecl, sinob, sxlong, arg, decdeg, djul, rjul, eccfac
   real(r8), dimension( 1:ncell, 1:nLevel, 1:1) :: f_ice_phy, f_rain_phy
   logical      :: doabsems
   real(r8)     :: swvisdir(1:ncell,1), swvisdif(1:ncell,1), swnirdir(1:ncell,1), swnirdif(1:ncell,1)
   real(r8)     :: xcoszen(ncell,1)

   xcoszen(1:ncell,1) = coszrs(1:ncell) 
!
! check time
!
   !call get_curr_calday(start_ymd, start_tod, 0, dtime, calday) 
    call get_curr_calday(start_ymd, start_tod, itimestep, dtime, calday) ! Regression, we use rrtmg now, the use of 'julday' is not
                                                                         ! essential; other phys may still need
   julday = calday

   stepra  = step_ra
   stepabs = stepra
   radt    = dtime*stepra/60._r8
!
! cloud fraction call every step
!
   select case(trim(wrfphys_cf_scheme))
   case('BINARY')
   call cal_cldfra(pstate_wrf%cldfra(1:ncell,1:nLevel,1:1)     , &  ! out
                   pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qc), &  ! in
                   pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qi), &  ! in
                   p_qi,  &                   ! in                        
                   p_qc,  &                   ! in
                   param_first_scalar,     &  ! in
                   1,ncell, 1,1, 1,nLevel, &
                   1,ncell, 1,1, 1,nLevel, &
                   1,ncell, 1,1, 1,nLevel  )

   case('RANDALL')
    call cal_cldfra2(cldfra=pstate_wrf%cldfra(1:ncell,1:nLevel,1:1)     , & ! out
                     qv    =pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qv), & ! in
                     qc    =pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qc), & ! in
                     qi    =pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qi), & ! in
                     qs    =pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qs), & ! in
                     f_qv=.true., f_qc=.true., f_qi=.true., f_qs=.true., &  
                     t_phy =pstate_wrf%t_phy( 1:ncell,1:nLevel,1:1),  & ! in 
                     p_phy =pstate_wrf%p_phy( 1:ncell,1:nLevel,1:1),  & ! in
                     ids=1,ide=ncell, jds=1,jde=1, kds=1,kde=nLevel,  &
                     ims=1,ime=ncell, jms=1,jme=1, kms=1,kme=nLevel,  &
                     its=1,ite=ncell, jts=1,jte=1, kts=1,kte=nLevel-1)
   case('CAM3')
    !call cam3_cldfrc(ncell, nLevel-1, 1, pstate_wrf%relhum(1:ncell,1:nLevel-1,1))
    call cam3_cldfrc(ncell, nLevel-1, 1)
    
   case default
     if(mpi_rank().eq.0) print*,"you must select a cloud fraction scheme in WRFphys: wrfphys_cf_scheme"
     call endrun()
   end select

   call_radiation = .false.

   if (itimestep.eq.0 .or. itimestep .eq. 1 .or.  mod(itimestep,stepra)  .eq. 0) call_radiation = .true.
   if (                    itimestep .ne. 1 .and. mod(itimestep,stepabs) .eq. 0) doabsems       = .true.

   IF(call_radiation)then
       !psurf_wrf%gsw          = zero
       !psurf_wrf%glw          = zero
       !psurf_wrf%swdown       = zero
       ptend_wrf%rthraten     = zero
       ptend_wrf%rthraten_lw  = zero
       ptend_wrf%rthraten_sw  = zero
!
! calculate constant for short wave radiation
!
       gmt = zero

       call radconst(xtime, &     ! out
                     declin,&     ! out
                     solcon,&     ! out
                     gmt,   &     ! in
                     calday,&     ! in
                     julian,&     ! out
                     deg2rad,   & ! in, constant
                     dpd,       & ! in
                     itimestep, & ! in
                     dtime )      ! in

!----------------------------------------------------------------------
!            LWRAD 
!----------------------------------------------------------------------

  select case(trim(wrfphys_ralw_scheme))
    case('RRTM')
       call rrtmlwrad(ptend_wrf%rthraten_lw(1:ncell,1:nLevel, 1:1),      &  ! inout, 1st-calc-inside
                      psurf_wrf%glw(     1:ncell,          1:1),         &  ! inout
                      psurf_wrf%emiss(   1:ncell,          1:1),         &  ! in
                      pstate_wrf%moist(  1:ncell,1:nLevel, 1:1, p_qv),   &  ! in
                      pstate_wrf%moist(  1:ncell,1:nLevel, 1:1, p_qc),   &  ! in
                      pstate_wrf%moist(  1:ncell,1:nLevel, 1:1, p_qr),   &  ! in
                      pstate_wrf%moist(  1:ncell,1:nLevel, 1:1, p_qi),   &  ! in
                      pstate_wrf%moist(  1:ncell,1:nLevel, 1:1, p_qs),   &  ! in
                      pstate_wrf%moist(  1:ncell,1:nLevel, 1:1, p_qg),   &  ! in
                      pstate_wrf%p8w(    1:ncell,1:nLevel, 1:1),         &  ! in, pressure at full levels (pa)
                      pstate_wrf%p_phy(  1:ncell,1:nLevel, 1:1),         &  ! in, pressure (pa)
                      pstate_wrf%pi_phy( 1:ncell,1:nLevel, 1:1),         &  ! in, exner function (dimensionless)
                      pstate_wrf%dz8w(   1:ncell,1:nLevel, 1:1),         &  ! in, dz between full levels (m)
                      pstate_wrf%t_phy(  1:ncell,1:nLevel, 1:1),         &  ! in, temperature (k)
                      pstate_wrf%t8w(    1:ncell,1:nLevel, 1:1),         &  ! in, temperature at full levels (k)
                      pstate_wrf%rhom(   1:ncell,1:nLevel, 1:1),         &  ! in, density (kg/m^3)
                      pstate_wrf%cldfra( 1:ncell,1:nLevel, 1:1),         &  ! in, cloud fraction (binary)
                      r_d,       &  ! constant
                      gravity,   &  ! constant
                      p_qv,      &  ! constant
                      p_qc,      &  ! constant
                      p_qr,      &  ! constant
                      p_qi,      &  ! constant
                      p_qs,      &  ! constant
                      p_qg,      &  ! constant
                      param_first_scalar,     & ! constant
                      icloud,                 & ! constant
                      warm_rain,              & ! flag
                      1,ncell, 1,1, 1,nLevel, &     
                      1,ncell, 1,1, 1,nLevel, &
                      1,ncell, 1,1, 1,nLevel-1)
    case ('CAM')
     ! do nothing until swrad, but note that  camlw is problematic now

    case('RRTMGV341')
! because v341 and v381 rrtmg has many common modules, compiling them
! is mutuallly exclusive, v341 is an old version for regression,
! test v381 now
#ifdef RRTMG_V341
      call rrtmg_lwrad(rthratenlw = ptend_wrf%rthraten_lw(1:ncell,1:nLevel, 1:1), & ! inout
                       lwupt      = pstate_wrf%lwupt (1:ncell,1:1), & ! diagnosed
                       lwuptc     = pstate_wrf%lwuptc(1:ncell,1:1), & ! diagnosed
                       lwdnt      = pstate_wrf%lwdnt (1:ncell,1:1), & ! diagnosed
                       lwdntc     = pstate_wrf%lwdntc(1:ncell,1:1), & ! diagnosed
                       lwupb      = pstate_wrf%lwupb (1:ncell,1:1), & ! diagnosed
                       lwupbc     = pstate_wrf%lwupbc(1:ncell,1:1), & ! diagnosed
                       lwdnb      = pstate_wrf%lwdnb (1:ncell,1:1), & ! diagnosed
                       lwdnbc     = pstate_wrf%lwdnbc(1:ncell,1:1), & ! diagnosed
                       glw        = psurf_wrf%glw    (1:ncell,1:1), & ! inout
                       olr        = pstate_wrf%olr   (1:ncell,1:1), & ! diagnosed
                       lwcf       = pstate_wrf%lwcf  (1:ncell,1:1), & ! diagnosed
                       emiss      = psurf_wrf%emiss(  1:ncell,          1:1), & ! in
                       p8w        = pstate_wrf%p8w(   1:ncell,1:nLevel, 1:1), & ! in
                       p3d        = pstate_wrf%p_phy( 1:ncell,1:nLevel, 1:1), &
                       pi3d       = pstate_wrf%pi_phy(1:ncell,1:nLevel, 1:1), &
                       dz8w       = pstate_wrf%dz8w(  1:ncell,1:nLevel, 1:1), &
                       tsk        = psurf_wrf%tsk(    1:ncell,            1), &
                       t3d        = pstate_wrf%t_phy( 1:ncell,1:nLevel, 1:1), &
                       t8w        = pstate_wrf%t8w(   1:ncell,1:nLevel, 1:1), &
                       rho3d      = pstate_wrf%rhom(  1:ncell,1:nLevel, 1:1), &
                       r          = r_d,       &
                       g          = gravity,   &
                       icloud     = icloud,    &
                       warm_rain  = warm_rain, &
                       cldfra3d   = pstate_wrf%cldfra( 1:ncell,1:nLevel, 1:1), &
                       f_ice_phy  = f_ice_phy, &
                       f_rain_phy = f_rain_phy,&
                       xland      = psurf_wrf%xland (  1:ncell,          1:1), &
                       xice       = psurf_wrf%xice  (  1:ncell,          1:1), &
                       snow       = psurf_wrf%snow  (  1:ncell,          1:1), &
                       qv3d       = pstate_wrf%moist(  1:ncell,1:nLevel, 1:1, p_qv), & 
                       qc3d       = pstate_wrf%moist(  1:ncell,1:nLevel, 1:1, p_qc), & 
                       qr3d       = pstate_wrf%moist(  1:ncell,1:nLevel, 1:1, p_qr), &
                       qi3d       = pstate_wrf%moist(  1:ncell,1:nLevel, 1:1, p_qi), & 
                       qs3d       = pstate_wrf%moist(  1:ncell,1:nLevel, 1:1, p_qs), & 
                       qg3d       = pstate_wrf%moist(  1:ncell,1:nLevel, 1:1, p_qg), &
                       f_qv       = .true., &
                       f_qc       = .true., &
                       f_qr       = .true., &
                       f_qi       = .true., &
                       f_qs       = .true., &
                       f_qg       = .true., &
                       ! Do not consider aerosol and WRF-CHEM now
                       !tauaerlw1,tauaerlw2,tauaerlw3,tauaerlw4,    & ! czhao 
                       !tauaerlw5,tauaerlw6,tauaerlw7,tauaerlw8,    & ! czhao 
                       !tauaerlw9,tauaerlw10,tauaerlw11,tauaerlw12, & ! czhao 
                       !tauaerlw13,tauaerlw14,tauaerlw15,tauaerlw16,& ! czhao 
                       !aer_ra_feedback,                            & !czhao
                       !progn,                                      & !czhao
                       !qndrop3d,f_qndrop,                          & !czhao
                       ids=1,ide=ncell, jds=1,jde=1, kds=1,kde=nLevel, & 
                       ims=1,ime=ncell, jms=1,jme=1, kms=1,kme=nLevel, &
                       its=1,ite=ncell, jts=1,jte=1, kts=1,kte=nLevel-1 )
                       ! Do not diagnose these now
                       !lwupflx, lwupflxc, lwdnflx, lwdnflxc
#endif

    case('RRTMGV381')
#ifdef RRTMG_V381
        call RRTMG_LWRAD(&
                       rthratenlw = ptend_wrf%rthraten_lw(1:ncell,1:nLevel, 1:1),   &
                       lwupt      = pstate_wrf%lwupt (1:ncell,1:1), &
                       lwuptc     = pstate_wrf%lwuptc(1:ncell,1:1), &
                       lwdnt      = pstate_wrf%lwdnt (1:ncell,1:1), &
                       lwdntc     = pstate_wrf%lwdntc(1:ncell,1:1), &
                       lwupb      = pstate_wrf%lwupb (1:ncell,1:1), &
                       lwupbc     = pstate_wrf%lwupbc(1:ncell,1:1), &
                       lwdnb      = pstate_wrf%lwdnb (1:ncell,1:1), &
                       lwdnbc     = pstate_wrf%lwdnbc(1:ncell,1:1), &
!                      lwupflx, lwupflxc, lwdnflx, lwdnflxc,      &
                       glw        =  psurf_wrf%glw   (1:ncell,1:1), & ! inout, 
                       olr        = pstate_wrf%olr   (1:ncell,1:1), &
                       lwcf       = pstate_wrf%lwcf  (1:ncell,1:1), &
                       emiss      =  psurf_wrf%emiss (1:ncell,          1:1), &
                       p8w        = pstate_wrf%p8w(   1:ncell,1:nLevel, 1:1), &
                       p3d        = pstate_wrf%p_phy( 1:ncell,1:nLevel, 1:1), &
                       pi3d       = pstate_wrf%pi_phy(1:ncell,1:nLevel, 1:1), &
                       dz8w       = pstate_wrf%dz8w(  1:ncell,1:nLevel, 1:1), &
                       tsk        = psurf_wrf%tsk(    1:ncell,          1:1), &
                       t3d        = pstate_wrf%t_phy( 1:ncell,1:nLevel, 1:1), &
                       t8w        = pstate_wrf%t8w(   1:ncell,1:nLevel, 1:1), &
                       rho3d      = pstate_wrf%rhom(  1:ncell,1:nLevel, 1:1), &
                       r          = r_d, g = gravity,                         &
                       icloud     = icloud, warm_rain = warm_rain,            &
                       cldfra3d   = pstate_wrf%cldfra( 1:ncell,1:nLevel, 1:1),&
                       ! comment these 2 because is_cammgmp_used = .false.
                       !lradius,iradius,                           & 
                       is_cammgmp_used = .false.,                           & 
                       !f_ice_phy, f_rain_phy,                     &
                       xland = psurf_wrf%xland (1:ncell,    1:1), &
                       xice  = psurf_wrf%xice  (1:ncell,    1:1), &
                       snow  = psurf_wrf%snow  (1:ncell,    1:1), &
                       qv3d  = pstate_wrf%moist(  1:ncell,1:nLevel, 1:1, p_qv), &
                       qc3d  = pstate_wrf%moist(  1:ncell,1:nLevel, 1:1, p_qc), &
                       qr3d  = pstate_wrf%moist(  1:ncell,1:nLevel, 1:1, p_qr), &
                       qi3d  = pstate_wrf%moist(  1:ncell,1:nLevel, 1:1, p_qi), &
                       qs3d  = pstate_wrf%moist(  1:ncell,1:nLevel, 1:1, p_qs), &
                       qg3d  = pstate_wrf%moist(  1:ncell,1:nLevel, 1:1, p_qg), &
                       ! varying O3, not used now, may use later
                       !o3input, 
                       !o33d,                             &
                       f_qv = .true., f_qc= .true., f_qr= .true., f_qi= .true., f_qs= .true., f_qg= .true., &
                       re_cloud = pstate_wrf%re_cloud(1:ncell,1:nLevel, 1:1), &
                       re_ice   = pstate_wrf%re_ice  (1:ncell,1:nLevel, 1:1), &
                       re_snow  = pstate_wrf%re_snow (1:ncell,1:nLevel, 1:1), &  ! G. Thompson
!
! if wphys_has_req .gt.0, but re_cloud/ice/snow was not calculated (i.e.,0), rad code
! will still compute a value for it, but not reicalc and relcalc
!
                       has_reqc = wphys_has_req, has_reqi = wphys_has_req, has_reqs = wphys_has_req, &  ! G. Thompson
                       !optional chem not used now
                       !tauaerlw1,tauaerlw2,tauaerlw3,tauaerlw4,   & ! czhao 
                       !tauaerlw5,tauaerlw6,tauaerlw7,tauaerlw8,   & ! czhao 
                       !tauaerlw9,tauaerlw10,tauaerlw11,tauaerlw12,   & ! czhao 
                       !tauaerlw13,tauaerlw14,tauaerlw15,tauaerlw16,   & ! czhao 
                       !aer_ra_feedback,                           & !czhao
!jdfcz                 progn,prescribe,                           & !czhao
                       !progn,                                     & !czhao
                       !qndrop3d,f_qndrop,                         & !czhao
!ccc added for time varying gases. yizhang: not used now, may use later
                       yr = 2000,julian = julian,                      &
!ccc
                       mp_physics = 999,  & ! not used by grist now, dum input
                       ids=1,ide=ncell,jds=1,jde=1, kds=1,kde=nLevel,  & 
                       ims=1,ime=ncell,jms=1,jme=1, kms=1,kme=nLevel,  &
                       its=1,ite=ncell,jts=1,jte=1, kts=1,kte=nLevel-1)
                       ! yizhang: these are flux at each layer, commented now
                       !lwupflx, lwupflxc, lwdnflx, lwdnflxc  )
#endif
    case default
      if(mpi_rank().eq.0) print*,"you must select a radiation LW scheme in WRFphys: wrfphys_ralw_scheme"
      call endrun()
    end select

!----------------------------------------------------------------------
!            SWRAD
!----------------------------------------------------------------------

  select case (trim(wrfphys_rasw_scheme))
  case('MM5')
       call swrad(dtime,                                         & ! in
                  ptend_wrf%rthraten_sw(1:ncell,1:nLevel,1:1)  , & ! inout, 1st-calc-inside
                  psurf_wrf%gsw(     1:ncell,         1:1)     , & ! inout, 1st-calc-inside
                  pstate_wrf%xlat(   1:ncell,         1:1)     , & ! in
                  pstate_wrf%xlong(  1:ncell,         1:1)     , & ! in
                  psurf_wrf%albedo(  1:ncell,1:nLevel    )     , & ! in
                  pstate_wrf%rhom (  1:ncell,1:nLevel,1:1)     , & ! in
                  pstate_wrf%t_phy(  1:ncell,1:nLevel,1:1)     , & ! in
                  pstate_wrf%moist(  1:ncell,1:nLevel,1:1,p_qv), & ! in
                  pstate_wrf%moist(  1:ncell,1:nLevel,1:1,p_qc), & ! in
                  pstate_wrf%moist(  1:ncell,1:nLevel,1:1,p_qr), & ! in
                  pstate_wrf%moist(  1:ncell,1:nLevel,1:1,p_qi), & ! in
                  pstate_wrf%moist(  1:ncell,1:nLevel,1:1,p_qs), & ! in
                  pstate_wrf%moist(  1:ncell,1:nLevel,1:1,p_qg), & ! in
                  pstate_wrf%p_phy(  1:ncell,1:nLevel,1:1)     , & ! in
                  pstate_wrf%pi_phy( 1:ncell,1:nLevel,1:1)     , & ! in
                  pstate_wrf%dz8w(   1:ncell,1:nLevel,1:1)     , & ! in
                  gmt,      & ! in, constant
                  r_d,      & ! in, constant
                  cp,       & ! in, constant
                  gravity,  & ! in, constant
                  julday,   & ! in, constant
                  xtime,    & ! in, constant
                  declin,   & ! in, constant
                  solcon,   & ! in, constant
                  p_qv,     & ! in, constant
                  p_qc,     & ! in, constant
                  p_qr,     & ! in, constant
                  p_qi,     & ! in, constant
                  p_qs,     & ! in, constant
                  p_qg,     & ! in, constant
                  param_first_scalar,   & ! in, constant
                  radt,                 & ! in, constant
                  icloud,               & ! in, constant
                  deg2rad,              & ! in, constant
                  warm_rain,            & ! in, constant
                  1,ncell, 1,1, 1,nLevel, &    
                  1,ncell, 1,1, 1,nLevel, & 
                  1,ncell, 1,1, 1,nLevel-1)

    case('CAM')
!
! CAMRAD sw can run, lw is problemtic now, use rrtm
!
       call camrad(rthratenlw=ptend_wrf%rthraten_lw(1:ncell,1:nLevel,1:1), &
                   rthratensw=ptend_wrf%rthraten_sw(1:ncell,1:nLevel,1:1), &
                   dolw=.false., dosw=.true.,               &
                   ! beg diagnostics
                   swupt  = pstate_wrf%swupt (1:ncell,1:1), &
                   swuptc = pstate_wrf%swuptc(1:ncell,1:1), &
                   swdnt  = pstate_wrf%swdnt (1:ncell,1:1), &
                   swdntc = pstate_wrf%swdntc(1:ncell,1:1), &
                   lwupt  = pstate_wrf%lwupt (1:ncell,1:1), &
                   lwuptc = pstate_wrf%lwuptc(1:ncell,1:1), &
                   lwdnt  = pstate_wrf%lwdnt (1:ncell,1:1), &
                   lwdntc = pstate_wrf%lwdntc(1:ncell,1:1), &
                   swupb  = pstate_wrf%swupb (1:ncell,1:1), &
                   swupbc = pstate_wrf%swupbc(1:ncell,1:1), &
                   swdnb  = pstate_wrf%swdnb (1:ncell,1:1), &
                   swdnbc = pstate_wrf%swdnbc(1:ncell,1:1), &
                   lwupb  = pstate_wrf%lwupb (1:ncell,1:1), &
                   lwupbc = pstate_wrf%lwupbc(1:ncell,1:1), &
                   lwdnb  = pstate_wrf%lwdnb (1:ncell,1:1), &
                   lwdnbc = pstate_wrf%lwdnbc(1:ncell,1:1), &
                   swcf   = pstate_wrf%swcf  (1:ncell,1:1), &
                   lwcf   = pstate_wrf%lwcf  (1:ncell,1:1), &
                   olr    = pstate_wrf%olr   (1:ncell,1:1), &
                   cemiss = pstate_wrf%cemiss (1:ncell, 1:nLevel,1:1),&
                   taucldc= pstate_wrf%taucldc(1:ncell, 1:nLevel,1:1),&
                   taucldi= pstate_wrf%taucldi(1:ncell, 1:nLevel,1:1),&
                   coszr  = pstate_wrf%coszr  (1:ncell,          1:1),&
                   ! end diagnostics
                   gsw    = psurf_wrf%gsw   (1:ncell,1:1), & ! inout, 1st out
                   glw    = psurf_wrf%glw   (1:ncell,1:1), & ! inout, 1st out
                   xlat   = pstate_wrf%xlat (1:ncell,1:1), & ! in
                   xlong  = pstate_wrf%xlong(1:ncell,1:1), & ! in
                   albedo = psurf_wrf%albedo(1:ncell,1:1), & ! in
                   t_phy  = pstate_wrf%t_phy(1:ncell,1:nLevel,1:1) ,& ! in
                   tsk    = psurf_wrf%tsk (  1:ncell,1)            ,& ! in
                   emiss  = psurf_wrf%emiss( 1:ncell,         1:1) ,     & ! in
                   qv3d   = pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qv), & ! in
                   qc3d   = pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qc), & ! in
                   qr3d   = pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qr), & ! in
                   qi3d   = pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qi), & ! in
                   qs3d   = pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qs), & ! in
                   qg3d   = pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qg), & ! in
                   alswvisdir = psurf_wrf%asdir(1:ncell,1),              &
                   alswvisdif = psurf_wrf%asdif(1:ncell,1),              &
                   alswnirdir = psurf_wrf%aldir(1:ncell,1),              &
                   alswnirdif = psurf_wrf%aldif(1:ncell,1),              &
                   swvisdir   = swvisdir(1:ncell,1),                     & ! out, not used yet
                   swvisdif   = swvisdif(1:ncell,1),                     & ! out, not used yet
                   swnirdir   = swnirdir(1:ncell,1),                     & ! out, not used yet
                   swnirdif   = swnirdif(1:ncell,1),                     & ! out, not used yet
                   sf_surface_physics = 8,                               & !
                   f_qv=.true.,f_qc=.true.,f_qr=.false.,f_qi=.true.,f_qs=.true.,f_qg=.false., &
                   f_ice_phy=f_ice_phy, f_rain_phy=f_rain_phy, & ! dum, not used
                   p_phy     = pstate_wrf%p_phy(   1:ncell,1:nLevel,1:1), & ! in
                   p8w       = pstate_wrf%p8w  (   1:ncell,1:nLevel,1:1), & ! in
                   z         = pstate_wrf%zzz  (   1:ncell,1:nLevel,1:1), & ! in
                   pi_phy    = pstate_wrf%pi_phy(  1:ncell,1:nLevel,1:1), & ! in
                   rho_phy   = pstate_wrf%rhom  (  1:ncell,1:nLevel,1:1), & ! in
                   dz8w      = pstate_wrf%dz8w  (  1:ncell,1:nLevel,1:1), & ! in
                   cldfra    = pstate_wrf%cldfra(  1:ncell,1:nLevel,1:1), & ! in
                   xland     = psurf_wrf%xland  (  1:ncell,         1:1), & ! in
                   xice      = psurf_wrf%xice   (  1:ncell,         1:1), & ! in
                   snow      = psurf_wrf%snow   (  1:ncell,         1:1), & ! in
                   ozmixm    = ozmixm(1:ncell, 1:levsiz,  1:1, 1:num_months ), & ! in
                   pin0      = pin,        &  ! in
                   levsiz    = levsiz,     &  ! in
                   num_months= num_months, &  ! in
                   m_psp     = m_psp(1:ncell,1:1), & ! in
                   m_psn     = m_psn(1:ncell,1:1), & ! in
                   aerosolcp = aerosolcp(1:ncell, 1:paerlev, 1:1, 1:naer_c ), & ! in
                   aerosolcn = aerosolcn(1:ncell, 1:paerlev, 1:1, 1:naer_c ), & ! in
                   m_hybi0   = m_hybi, & ! in
                   cam_abs_dim1=cam_abs_dim1, cam_abs_dim2=cam_abs_dim2, paerlev=paerlev,naer_c=naer_c, & ! in, dimension 
                   gmt=gmt,julday=julday,julian=julian,yr=2000,dt=dtime,xtime=xtime,declin=declin,solcon=solcon,& ! (quasi) constants
                   radt      = radt, degrad = deg2rad, n_cldadv=n_cldadv,   &
                   abstot_3d = abstot_3d,   & ! inout, mainly out 
                   absnxt_3d = absnxt_3d,   & ! inout
                   emstot_3d = emstot_3d,   & ! inout
                   doabsems  = doabsems,    &
                   ids=1,ide = ncell, jds=1,jde=1, kds=1,kde=nLevel,  &
                   ims=1,ime = ncell, jms=1,jme=1, kms=1,kme=nLevel,  &
                   its=1,ite = ncell, jts=1,jte=1, kts=1,kte=nLevel-1   )

        psurf_wrf%swdown(1:ncell,1:1)  = pstate_wrf%swdnb (1:ncell,1:1)

    case('RRTMGV341')
#ifdef RRTMG_V341
      call rrtmg_swrad(coszrs, rthratensw  = ptend_wrf%rthraten_sw(1:ncell,1:nLevel,1:1), &
                        swupt       = pstate_wrf%swupt (1:ncell,1:1), & ! diagnose
                        swuptc      = pstate_wrf%swuptc(1:ncell,1:1), & ! diagnose
                        swdnt       = pstate_wrf%swdnt (1:ncell,1:1), & ! diagnose
                        swdntc      = pstate_wrf%swdntc(1:ncell,1:1), & ! diagnose
                        swupb       = pstate_wrf%swupb (1:ncell,1:1), & ! diagnose
                        swupbc      = pstate_wrf%swupbc(1:ncell,1:1), & ! diagnose
                        swdnb       = pstate_wrf%swdnb (1:ncell,1:1), & ! diagnose
                        swdnbc      = pstate_wrf%swdnbc(1:ncell,1:1), & ! diagnose
                        swcf        = pstate_wrf%swcf  (1:ncell,1:1), & ! diagnose
                        gsw         = psurf_wrf%gsw    (1:ncell,1:1), &
                        xtime       = xtime,                          & ! local
                        gmt         = gmt,                            & ! local
                        xlat        = pstate_wrf%xlat (1:ncell,1:1),  &
                        xlong       = pstate_wrf%xlong(1:ncell,1:1),  &
                        radt        = radt,                           &
                        degrad      = deg2rad,                        &
                        declin      = declin,                         &
                        coszr       = pstate_wrf%coszr  (1:ncell,1:1),& ! inout 
                        julday      = julday,                         &
                        solcon      = solcon,                         &
                        albedo      = psurf_wrf%albedo(1:ncell,1:1),  &
                        t3d         = pstate_wrf%t_phy( 1:ncell,1:nLevel, 1:1), &
                        t8w         = pstate_wrf%t8w(   1:ncell,1:nLevel, 1:1), &
                        tsk         = psurf_wrf%tsk (   1:ncell,1),             &
                        p3d         = pstate_wrf%p_phy( 1:ncell,1:nLevel, 1:1), &
                        p8w         = pstate_wrf%p8w(   1:ncell,1:nLevel, 1:1), &
                        pi3d        = pstate_wrf%pi_phy(1:ncell,1:nLevel, 1:1), &
                        rho3d       = pstate_wrf%rhom(  1:ncell,1:nLevel, 1:1), &
                        dz8w        = pstate_wrf%dz8w(  1:ncell,1:nLevel, 1:1), &
                        cldfra3d    = pstate_wrf%cldfra(  1:ncell,1:nLevel,1:1),&
                        r           = r_d,          &
                        g           = gravity,      &
                        icloud      = icloud,       &
                        warm_rain   = warm_rain,    &
                        f_ice_phy   = f_ice_phy,    & ! dum, not used
                        f_rain_phy  = f_rain_phy,   & ! dum, not used
                        xland       = psurf_wrf%xland(  1:ncell,          1:1),    &
                        xice        = psurf_wrf%xice (  1:ncell,          1:1),    &
                        snow        = psurf_wrf%snow (  1:ncell,          1:1),    &
                        qv3d        = pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qv), &
                        qc3d        = pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qc), &
                        qr3d        = pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qr), &
                        qi3d        = pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qi), &
                        qs3d        = pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qs), &
                        qg3d        = pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qg), &
                        ! consider as noahmp as weell
                        alswvisdir  = psurf_wrf%asdir(1:ncell,1)                 , &
                        alswvisdif  = psurf_wrf%asdif(1:ncell,1)                 , &
                        alswnirdir  = psurf_wrf%aldir(1:ncell,1)                 , &
                        alswnirdif  = psurf_wrf%aldif(1:ncell,1)                 , &
                        ! do not consider ssib output now
                        !swvisdir, swvisdif,                        &  !Zhenxin ssib swr comp (06/20/2011)
                        !swnirdir, swnirdif,                        &  !Zhenxin ssib swi comp (06/20/2011)
                        sf_surface_physics= 8,&  ! thou not ssib, still use those four albedo vars
                        f_qv        = .true., &
                        f_qc        = .true., &
                        f_qr        = .true., &
                        f_qi        = .true., &
                        f_qs        = .true., &
                        f_qg        = .true., &
                        ! do not consider aerosol and wrf-chem now
                        !tauaer300,tauaer400,tauaer600,tauaer999,   & ! czhao 
                        !gaer300,gaer400,gaer600,gaer999,           & ! czhao 
                        !waer300,waer400,waer600,waer999,           & ! czhao 
                        !aer_ra_feedback,                           &
                        !progn,                                     &
                        !qndrop3d,f_qndrop,                         & !czhao
                        ids=1,ide=ncell, jds=1,jde=1, kds=1,kde=nLevel,  & 
                        ims=1,ime=ncell, jms=1,jme=1, kms=1,kme=nLevel,  &
                        its=1,ite=ncell, jts=1,jte=1, kts=1,kte=nLevel-1 )
                        ! do not diagnose these now
                        !swupflx, swupflxc, swdnflx, swdnflxc       &

        psurf_wrf%swdown(1:ncell,1:1)  = pstate_wrf%swdnb (1:ncell,1:1)
#endif

    case('RRTMGV381')
#ifdef RRTMG_V381
      call RRTMG_SWRAD(rthratensw = ptend_wrf%rthraten_sw(1:ncell,1:nLevel,1:1), &
                       swupt      = pstate_wrf%swupt (1:ncell,1:1), & ! diagnose
                       swuptc     = pstate_wrf%swuptc(1:ncell,1:1), & ! diagnose
                       swdnt      = pstate_wrf%swdnt (1:ncell,1:1), & ! diagnose
                       swdntc     = pstate_wrf%swdntc(1:ncell,1:1), & ! diagnose
                       swupb      = pstate_wrf%swupb (1:ncell,1:1), & ! diagnose
                       swupbc     = pstate_wrf%swupbc(1:ncell,1:1), & ! diagnose
                       swdnb      = pstate_wrf%swdnb (1:ncell,1:1), & ! diagnose
                       swdnbc     = pstate_wrf%swdnbc(1:ncell,1:1), & ! diagnose
!                      swupflx, swupflxc, swdnflx, swdnflxc,        &
                       swcf       = pstate_wrf%swcf  (1:ncell,1:1), & ! diagnose
                       gsw        = psurf_wrf%gsw    (1:ncell,1:1), & ! diagnose
                       xtime      = xtime, gmt  = gmt,              & ! use still
                       xlat       = pstate_wrf%xlat (1:ncell,1:1),  &
                       xlong      = pstate_wrf%xlong(1:ncell,1:1),  &
                       radt       = radt,                           &
                       degrad     = deg2rad,                        &
                       declin     = declin,                         &
                       coszr      = pstate_wrf%coszr  (1:ncell,1:1),& ! inout, not used
                       julday     = julday,      & ! not used
                       solcon     = solcon,      & ! input
                       albedo     = psurf_wrf%albedo(1:ncell,1:1),            &
                       t3d        = pstate_wrf%t_phy( 1:ncell,1:nLevel, 1:1), &
                       t8w        = pstate_wrf%t8w  ( 1:ncell,1:nLevel, 1:1), &
                       tsk        = psurf_wrf%tsk   ( 1:ncell,          1:1), &
                       p3d        = pstate_wrf%p_phy( 1:ncell,1:nLevel, 1:1), &
                       p8w        = pstate_wrf%p8w  ( 1:ncell,1:nLevel, 1:1), &
                       pi3d       = pstate_wrf%pi_phy(1:ncell,1:nLevel, 1:1), &
                       rho3d      = pstate_wrf%rhom ( 1:ncell,1:nLevel, 1:1), &
                       dz8w       = pstate_wrf%dz8w ( 1:ncell,1:nLevel, 1:1), &
                       cldfra3d   = pstate_wrf%cldfra(1:ncell,1:nLevel, 1:1), &
                       !lradius, iradius,          &
                       is_cammgmp_used = .false.,  &
                       r          = r_d, g = gravity,                          &
                       re_cloud   = pstate_wrf%re_cloud(1:ncell,1:nLevel,1:1), &
                       re_ice     = pstate_wrf%re_ice  (1:ncell,1:nLevel,1:1), &
                       re_snow    = pstate_wrf%re_snow (1:ncell,1:nLevel,1:1), &
                       has_reqc   = wphys_has_req, has_reqi = wphys_has_req, has_reqs = wphys_has_req, &
                       icloud     = icloud, warm_rain  = warm_rain,            &
                       !f_ice_phy , f_rain_phy,                     &
                       xland      = psurf_wrf%xland( 1:ncell,         1:1),    & 
                       xice       = psurf_wrf%xice(  1:ncell,         1:1),    & 
                       snow       = psurf_wrf%snow(  1:ncell,         1:1),    &
                       qv3d       = pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qv), &
                       qc3d       = pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qc), &
                       qr3d       = pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qr), &
                       qi3d       = pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qi), &
                       qs3d       = pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qs), &
                       qg3d       = pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qg), &
                       ! comment now
                       !o3input, o33d,                            &
                       !aer_opt    = 0, &! thou optional, required, or seg fault, due to internal inconsistency 
                       !aerod, 
                       no_src = 1,                   & ! not active because aer_opt and aerod are not input, dummy now
                       alswvisdir = psurf_wrf%asdir(1:ncell,1),   & 
                       alswvisdif = psurf_wrf%asdif(1:ncell,1),   &  !Zhenxin ssib alb comp (06/20/2011)
                       alswnirdir = psurf_wrf%aldir(1:ncell,1),   &
                       alswnirdif = psurf_wrf%aldif(1:ncell,1),   &  !Zhenxin ssib alb comp (06/20/2011)
                       ! do not consider ssib output
                       !swvisdir, swvisdif,                        &  !Zhenxin ssib swr comp (06/20/2011)
                       !swnirdir, swnirdif,                        &  !Zhenxin ssib swi comp (06/20/2011)
                       sf_surface_physics = 8,                     &  !Zhenxin
                       f_qv = .true., f_qc= .true., f_qr= .true.,  &
                       f_qi = .true., f_qs= .true., f_qg= .true.,  & ! true only for WSM6, when use other MP, should check, seem ok as well 
                       !tauaer300,tauaer400,tauaer600,tauaer999,   & ! czhao
                       !gaer300,gaer400,gaer600,gaer999,           & ! czhao
                       !waer300,waer400,waer600,waer999,           & ! czhao
                       !aer_ra_feedback,                           &
!jdfcz                 progn,prescribe,                            &
                       !progn,                                     &
                       !qndrop3d,f_qndrop,                         & ! czhao
                       mp_physics = 999,                           & ! dummy
                       ids=1,ide=ncell, jds=1,jde=1, kds=1,kde=nLevel,                 & 
                       ims=1,ime=ncell, jms=1,jme=1, kms=1,kme=nLevel,                 &
                       its=1,ite=ncell, jts=1,jte=1, kts=1,kte=nLevel-1,               &
                       xcoszen = xcoszen(1:ncell,1:1) ) !,julian = julday     )
                       !swupflx, swupflxc, swdnflx, swdnflxc,      &
                       !tauaer3d_sw,ssaaer3d_sw,asyaer3d_sw,       & ! jararias 2013/11
                       !swddir, swddni, swddif,                    & ! jararias 2013/08

      psurf_wrf%swdown(1:ncell,1:1)  = pstate_wrf%swdnb (1:ncell,1:1)
#endif

    case default
      if(mpi_rank().eq.0) print*,"you must select a radiation SW scheme in WRFphys: wrfphys_rasw_scheme"
      call endrun
    end select

! final renewed tendency
! must be here, or affect restart
    ptend_wrf%rthraten(1:ncell,1:nLevel,1:1) = ptend_wrf%rthraten_lw(1:ncell,1:nLevel,1:1)+ptend_wrf%rthraten_sw(1:ncell,1:nLevel,1:1)

    END IF

    return
   end subroutine grist_wrf_radiation_run

!---------------------------------------------------------------------
!bop
! !iroutine: radconst - compute radiation terms
! !interfac:
   subroutine radconst(xtime,declin,solcon,gmt,julday,julian, &
                       degrad,dpd,step,dt                            )
   use grist_zenith,                     only: orb_decl
!---------------------------------------------------------------------
!   use module_wrf_error
   implicit none
!---------------------------------------------------------------------

! !arguments:
   real, intent(in   )      ::       julday
   integer, intent(in   )   ::       step
   real, intent(in   )      ::       gmt,dt,degrad,dpd
   real, intent(out  )      ::       xtime,declin,solcon,julian
   real                     ::       obecl,sinob,sxlong,arg,  &
                                     decdeg,djul,rjul,eccfac
!
! !description:
! compute terms used in radiation physics
!eop

! for short wave radiation

   declin=0.
   solcon=0.

!-----obecl : obliquity = 23.5 degree.
        
   obecl=23.5*degrad
   sinob=sin(obecl)
   xtime=float(step)*dt/60. ! model time since start in mitutes
        
!-----calculate longitude of the sun from vernal equinox:
        
   !julian=float(julday-1)+(xtime/60.+gmt)/24.
   julian=julday  ! zhangyi changes for grist

   if(julian.ge.80.)sxlong=dpd*(julian-80.)
   if(julian.lt.80.)sxlong=dpd*(julian+285.)
   sxlong=sxlong*degrad
   arg=sinob*sin(sxlong)
   declin=asin(arg)
   decdeg=declin/degrad
!----solar constant eccentricity factor (paltridge and platt 1976)
   djul=julian*360./365.
   rjul=djul*degrad
   eccfac=1.000110+0.034221*cos(rjul)+0.001280*sin(rjul)+0.000719*  &
          cos(2*rjul)+0.000077*sin(2*rjul)
   solcon=1370.*eccfac
#ifdef AMIPW_CLIMATE
   solcon=1365.*eccfac
#endif
!
! For AP, be consistent with those in grist_zenith/CAM5 Physics
! otherewise use default for regresssing old wrfphys run
!
   if(doAquaPlanet)then
      call orb_decl(julian,declin,eccfac,doAquaPlanet)
      solcon = 1366._r8*eccfac
   end if
#ifdef AMIPW_CLIMATE
! overwrite previous vars, to be consistent with other parts (dtp, lsm)
   call orb_decl(julian,declin,eccfac,doAquaPlanet)
   solcon = 1365._r8*eccfac
#endif

!   write(6,10) decdeg, solcon
!10 format(1x,'*** solar declination angle = ',f6.2,' degrees.',     &
!        ' solar constant = ',f8.2,' w/m**2 ***')
!   call wrf_debug (50, wrf_err_message)

   end subroutine radconst

!---------------------------------------------------------------------
!bop
! !iroutine: cal_cldfra - compute cloud fraction
! !interface:
   subroutine cal_cldfra(cldfra,qc,qi,p_qi,p_qc,                     &
          param_first_scalar,                                        &
          ids,ide, jds,jde, kds,kde,                                 &
          ims,ime, jms,jme, kms,kme,                                 &
          its,ite, jts,jte, kts,kte                                  )
!---------------------------------------------------------------------
   implicit none
!---------------------------------------------------------------------
   integer,  intent(in   )   ::           ids,ide, jds,jde, kds,kde, &
                                          ims,ime, jms,jme, kms,kme, &
                                          its,ite, jts,jte, kts,kte

   integer,  intent(in   )   ::           p_qi,p_qc,param_first_scalar
!
   real, dimension( ims:ime, kms:kme, jms:jme ), intent(out  ) ::    &
                                                             cldfra

   real, dimension( ims:ime, kms:kme, jms:jme ), intent(in   ) ::    &
                                                                 qi, &
                                                                 qc

   real thresh
   integer:: i,j,k
! !description:
! compute cloud fraction from input ice and cloud water fields
! if provided.
!
! whether qi or qc is active or not is determined from the indices of
! the fields into the 4d scalar arrays in wrf. these indices are 
! p_qi and p_qc, respectively, and they are passed in to the routine
! to enable testing to see if qi and qc represent active fields in
! the moisture 4d scalar array carried by wrf.
! 
! if a field is active its index will have a value greater than or
! equal to param_first_scalar, which is also an input argument to 
! this routine.
!eop
!---------------------------------------------------------------------
     thresh=1.0e-6

     if ( p_qi .ge. param_first_scalar .and. p_qc .ge. param_first_scalar ) then
        do j = jts,jte
        do k = kts,kte
        do i = its,ite
           if ( qc(i,k,j)+qi(i,k,j) .gt. thresh) then
              cldfra(i,k,j)=1.
           else
              cldfra(i,k,j)=0.
           endif
        enddo
        enddo
        enddo
     else if ( p_qc .ge. param_first_scalar ) then
        do j = jts,jte
        do k = kts,kte
        do i = its,ite
           if ( qc(i,k,j) .gt. thresh) then
              cldfra(i,k,j)=1.
           else
              cldfra(i,k,j)=0.
           endif
        enddo
        enddo
        enddo
     else 
        do j = jts,jte
        do k = kts,kte
        do i = its,ite
           cldfra(i,k,j)=0.
        enddo
        enddo
        enddo
     endif
    return
   end subroutine cal_cldfra


!bop
! !iroutine: cal_cldfra2 - compute cloud fraction
! !interface:
! cal_cldfra_xr - compute cloud fraction.
! code adapted from that in module_ra_gfdleta.f in wrf_v2.0.3 by james done
!!
!!---  cloud fraction parameterization follows randall, 1994
!!     (see hong et al., 1998)
!!     (modified by ferrier, feb '02)
!

   subroutine cal_cldfra2(cldfra, qv, qc, qi, qs,                &
                          f_qv, f_qc, f_qi, f_qs, t_phy, p_phy,  &
                          f_ice_phy,f_rain_phy,                  &
                          ids,ide, jds,jde, kds,kde,             &
                          ims,ime, jms,jme, kms,kme,             &
                          its,ite, jts,jte, kts,kte              )
!---------------------------------------------------------------------
   implicit none
!---------------------------------------------------------------------
   integer,  intent(in   )   ::           ids,ide, jds,jde, kds,kde, &
                                          ims,ime, jms,jme, kms,kme, &
                                          its,ite, jts,jte, kts,kte

!
   real, dimension( ims:ime, kms:kme, jms:jme ), intent(out  ) ::    &
                                                             cldfra

   real, dimension( ims:ime, kms:kme, jms:jme ), intent(in   ) ::    &
                                                                 qv, &
                                                                 qi, &
                                                                 qc, &
                                                                 qs, &
                                                              t_phy, &
                                                              p_phy
!                                                              p_phy, &
!                                                          f_ice_phy, &
!                                                         f_rain_phy

   real, dimension( ims:ime, kms:kme, jms:jme ),                     &
         optional,                                                   &
         intent(in   ) ::                                            &
                                                          f_ice_phy, &
                                                         f_rain_phy

   logical,optional,intent(in) :: f_qc,f_qi,f_qv,f_qs

!  real thresh
   integer:: i,j,k
   real    :: rhum, tc, esw, esi, weight, qvsw, qvsi, qvs_weight, qimid, qwmid, qcld, denom, arg, subsat

   real    ,parameter :: alpha0=100., gamma=0.49, qcldmin=1.e-12,    &
                                        pexp=0.25, rhgrid=1.0
   real    , parameter ::  svp1=0.61078
   real    , parameter ::  svp2=17.2693882
   real    , parameter ::  svpi2=21.8745584
   real    , parameter ::  svp3=35.86
   real    , parameter ::  svpi3=7.66
   real    , parameter ::  svpt0=273.15
   real    , parameter ::  r_d = 287.
   real    , parameter ::  r_v = 461.6
   real    , parameter ::  ep_2=r_d/r_v
! !description:
! compute cloud fraction from input ice and cloud water fields
! if provided.
!
! whether qi or qc is active or not is determined from the indices of
! the fields into the 4d scalar arrays in wrf. these indices are 
! p_qi and p_qc, respectively, and they are passed in to the routine
! to enable testing to see if qi and qc represent active fields in
! the moisture 4d scalar array carried by wrf.
! 
! if a field is active its index will have a value greater than or
! equal to param_first_scalar, which is also an input argument to 
! this routine.
!eop


!-----------------------------------------------------------------------
!---  compute grid-scale cloud cover for radiation
!
!     (modified by ferrier, feb '02)
!
!---  cloud fraction parameterization follows randall, 1994
!     (see hong et al., 1998)
!-----------------------------------------------------------------------
! note: ep_2=287./461.6 rd/rv
! note: r_d=287.

! alternative calculation for critical rh for grid saturation
!     rhgrid=0.90+.08*((100.-dx)/95.)**.5

! calculate saturation mixing ratio weighted according to the fractions of
! water and ice.
! following:
! murray, f.w. 1966. ``on the computation of saturation vapor pressure''  j. appl. meteor.  6 p.204
!    es (in mb) = 6.1078 . exp[ a . (t-273.16)/ (t-b) ]
!
!       over ice        over water
! a =   21.8745584      17.2693882
! b =   7.66            35.86

!---------------------------------------------------------------------

    do j = jts,jte
    do k = kts,kte
    do i = its,ite
      tc         = t_phy(i,k,j) - svpt0
      esw     = 1000.0 * svp1 * exp( svp2  * tc / ( t_phy(i,k,j) - svp3  ) )
      esi     = 1000.0 * svp1 * exp( svpi2 * tc / ( t_phy(i,k,j) - svpi3 ) )
      qvsw = ep_2 * esw / ( p_phy(i,k,j) - esw )
      qvsi = ep_2 * esi / ( p_phy(i,k,j) - esi )

      if ( present(f_qi) .and. present(f_qc) .and. present(f_qs) ) then

! mji - for mp options 2, 4, 6, 7, 8, etc. (qc = liquid, qi = ice, qs = snow)
         if ( f_qi .and. f_qc .and. f_qs) then
            qcld = qi(i,k,j)+qc(i,k,j)+qs(i,k,j)
            if (qcld .lt. qcldmin) then
               weight = 0.
            else
               weight = (qi(i,k,j)+qs(i,k,j)) / qcld
            endif
         endif

! mji - for mp options 1 and 3, (qc only)
!  for mp=1, qc = liquid, for mp=3, qc = liquid or ice depending on temperature
         if ( f_qc .and. .not. f_qi .and. .not. f_qs ) then
            qcld = qc(i,k,j)
            if (qcld .lt. qcldmin) then
               weight = 0.
            else
               if (t_phy(i,k,j) .gt. 273.15) weight = 0.
               if (t_phy(i,k,j) .le. 273.15) weight = 1.
            endif
         endif

! mji - for mp option 5; (qc = liquid, qs = ice)
         if ( f_qc .and. .not. f_qi .and. f_qs .and. present(f_ice_phy) ) then

! mixing ratios of cloud water & total ice (cloud ice + snow).
! mixing ratios of rain are not considered in this scheme.
! f_ice is fraction of ice
! f_rain is fraction of rain

           qimid = qs(i,k,j)
           qwmid = qc(i,k,j)
! old method
!           qimid = qc(i,k,j)*f_ice_phy(i,k,j)
!           qwmid = (qc(i,k,j)-qimid)*(1.-f_rain_phy(i,k,j))
!
!--- total "cloud" mixing ratio, qcld.  rain is not part of cloud,
!    only cloud water + cloud ice + snow
!
           qcld=qwmid+qimid
           if (qcld .lt. qcldmin) then
              weight = 0.
           else
              weight = f_ice_phy(i,k,j)
           endif
         endif

      else
         cldfra(i,k,j)=0.

      endif !  if ( f_qi .and. f_qc .and. f_qs)


      qvs_weight = (1-weight)*qvsw + weight*qvsi
      rhum=qv(i,k,j)/qvs_weight   !--- relative humidity
!
!--- determine cloud fraction (modified from original algorithm)
!
      if (qcld .lt. qcldmin) then
!
!--- assume zero cloud fraction if there is no cloud mixing ratio
!
        cldfra(i,k,j)=0.
      elseif(rhum.ge.rhgrid)then
!
!--- assume cloud fraction of unity if near saturation and the cloud
!    mixing ratio is at or above the minimum threshold
!
        cldfra(i,k,j)=1.
      else
!
!--- adaptation of original algorithm (randall, 1994; zhao, 1995)
!    modified based on assumed grid-scale saturation at rh=rhgrid.
!
        subsat=max(1.e-10,rhgrid*qvs_weight-qv(i,k,j))
        denom=(subsat)**gamma
        arg=max(-6.9, -alpha0*qcld/denom)    ! <-- exp(-6.9)=.001
! prevent negative values  (new)
        rhum=max(1.e-10, rhum)
        cldfra(i,k,j)=(rhum/rhgrid)**pexp*(1.-exp(arg))
!!              arg=-1000*qcld/(rhum-rhgrid)
!!              arg=max(arg, argmin)
!!              cldfra(i,k,j)=(rhum/rhgrid)*(1.-exp(arg))
        if (cldfra(i,k,j) .lt. .01) cldfra(i,k,j)=0.
      endif          !--- end if (qcld .lt. qcldmin) ...
    enddo          !--- end do i
    enddo          !--- end do k
    enddo          !--- end do j

   end subroutine cal_cldfra2

end module grist_wrf_radiation_driver
