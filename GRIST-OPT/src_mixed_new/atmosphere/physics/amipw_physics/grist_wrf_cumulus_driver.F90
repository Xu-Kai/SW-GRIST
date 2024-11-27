!----------------------------------------------------------------------------
! Created on 2019
! Author: Yi Zhang
! Version 1.0
! Description: This driver mainly follows the naming convention of WRF, 
!              but the dim conversion requires a transpose from GRIST based 
!              vars, currently only support KF cumulus
! Revision history:
!--------------------------------------------------------------------------

module grist_wrf_cumulus_driver

   use grist_constants,          only: i4, r8, rdry, rvap, cp, cv, gravity, zero,r_v=>rvap
   use grist_wrf_data_structure, only: pstate_wrf, psurf_wrf, ptend_wrf, restart, warm_rain
   use grist_wrf_data_structure, only: p_qv, p_qc, p_qr, p_qi, p_qs, p_qg, param_first_scalar
   use grist_hpe_constants,      only: eta_full
   use grist_wrfphys_nml_module, only: wrfphys_cu_scheme, step_cu, unuse_cu
   use grist_nml_module,         only: nh_dynamics
   use grist_mpi
!
! module physics 
!
   use module_wrfmodel_constants,only: xlv, xlv0, xlv1, xls0, xls1, ep_1, ep_2, svp1, svp2, svp3, svpt0
   use module_cu_ntiedtke_v381,  only: cu_ntiedtke_v381=>cu_ntiedtke, cu_ntiedtke_v381_init=>ntiedtkeinit
   use module_cu_gf_v381,        only: gfdrv_v381=>gfdrv
   use module_cu_gf_v402,        only: gfdrv_v402=>gfdrv

  implicit none

  private

  public  ::  grist_wrf_cumulus_init, &
              grist_wrf_cumulus_run , &
              grist_wrf_cumulus_final

    real(r8), allocatable   :: nca(:, :)
    logical,  allocatable   :: cu_act_flag(:,:)
    integer(i4)             :: stepcu   !# of fundamental timesteps between convection calls

contains

   subroutine grist_wrf_cumulus_init(ncell, nLevel)
    integer(i4),  intent(in) :: ncell, nLevel

      stepcu  = step_cu

      if(.not.allocated(nca))                allocate(nca(1:ncell, 1:1))
      if(.not.allocated(cu_act_flag))        allocate(cu_act_flag(1:ncell, 1:1))
      if(.not.allocated(ptend_wrf%ruucuten)) allocate(ptend_wrf%ruucuten(1:ncell, 1:nLevel, 1:1),source=zero)
      if(.not.allocated(ptend_wrf%rvvcuten)) allocate(ptend_wrf%rvvcuten(1:ncell, 1:nLevel, 1:1),source=zero)
      if(.not.allocated(ptend_wrf%rthcuten)) allocate(ptend_wrf%rthcuten(1:ncell, 1:nLevel, 1:1),source=zero)
      if(.not.allocated(ptend_wrf%rqvcuten)) allocate(ptend_wrf%rqvcuten(1:ncell, 1:nLevel, 1:1),source=zero)
      if(.not.allocated(ptend_wrf%rqccuten)) allocate(ptend_wrf%rqccuten(1:ncell, 1:nLevel, 1:1),source=zero)
      if(.not.allocated(ptend_wrf%rqrcuten)) allocate(ptend_wrf%rqrcuten(1:ncell, 1:nLevel, 1:1),source=zero)
      if(.not.allocated(ptend_wrf%rqicuten)) allocate(ptend_wrf%rqicuten(1:ncell, 1:nLevel, 1:1),source=zero)
      if(.not.allocated(ptend_wrf%rqscuten)) allocate(ptend_wrf%rqscuten(1:ncell, 1:nLevel, 1:1),source=zero)

      select case(trim(wrfphys_cu_scheme))

      case('NTDKV381')
      call cu_ntiedtke_v381_init(rthcuten= ptend_wrf%rthcuten,&
                        rqvcuten= ptend_wrf%rqvcuten,&
                        rqccuten= ptend_wrf%rqccuten,&
                        rqicuten= ptend_wrf%rqicuten,&
                        rucuten = ptend_wrf%ruucuten,&
                        rvcuten = ptend_wrf%rvvcuten,&
                        ! check whether restart needs this later
                        !rthften = ,
                        !rqvften = ,                 &
                        restart = restart,                               &
                        p_qc    = p_qv,                                  &
                        p_qi    = p_qi,                                  &
                        p_first_scalar  = param_first_scalar,            &
                        allowed_to_read = .false.,                       &
                        ids=1,ide=ncell, jds=1,jde=1, kds=1, kde=nLevel, &
                        ims=1,ime=ncell, jms=1,jme=1, kms=1, kme=nLevel, &
                        its=1,ite=ncell, jts=1,jte=1, kts=1, kte=nLevel-1)
      case('GFV381','GFV402')
!
! g3init, gdinit, gfinit just initialize some vars to zero, 
! we can do our own here
!
          ptend_wrf%ruucuten(1:ncell,1:nLevel,1:1) = zero
          ptend_wrf%rvvcuten(1:ncell,1:nLevel,1:1) = zero
          ptend_wrf%rthcuten(1:ncell,1:nLevel,1:1) = zero
          ptend_wrf%rqvcuten(1:ncell,1:nLevel,1:1) = zero
          ptend_wrf%rqccuten(1:ncell,1:nLevel,1:1) = zero
          ptend_wrf%rqrcuten(1:ncell,1:nLevel,1:1) = zero
          ptend_wrf%rqicuten(1:ncell,1:nLevel,1:1) = zero
          ptend_wrf%rqscuten(1:ncell,1:nLevel,1:1) = zero

      case default
           if(mpi_rank().eq.0) print*,"you must set wrfphys_cu_scheme, available options: TDK, NTDKV381"
           call mpi_abort()
      end select

      return
   end subroutine grist_wrf_cumulus_init

   subroutine grist_wrf_cumulus_final

      if(allocated(nca))                deallocate(nca)
      if(allocated(cu_act_flag))        deallocate(cu_act_flag)
      if(allocated(ptend_wrf%ruucuten)) deallocate(ptend_wrf%ruucuten)
      if(allocated(ptend_wrf%rvvcuten)) deallocate(ptend_wrf%rvvcuten)
      if(allocated(ptend_wrf%rthcuten)) deallocate(ptend_wrf%rthcuten)
      if(allocated(ptend_wrf%rqvcuten)) deallocate(ptend_wrf%rqvcuten)
      if(allocated(ptend_wrf%rqccuten)) deallocate(ptend_wrf%rqccuten)
      if(allocated(ptend_wrf%rqrcuten)) deallocate(ptend_wrf%rqrcuten)
      if(allocated(ptend_wrf%rqicuten)) deallocate(ptend_wrf%rqicuten)
      if(allocated(ptend_wrf%rqscuten)) deallocate(ptend_wrf%rqscuten)

    return
   end subroutine grist_wrf_cumulus_final
  
   subroutine grist_wrf_cumulus_run(ncell,nLevel,nspecies,itimestep,dtime)

!======================================================================
! grid structure in physics part of wrf
!----------------------------------------------------------------------
! the horizontal velocities used in the physics are unstaggered
! relative to temperature/moisture variables. all predicted
! variables are carried at half levels except w, which is at full
! levels. some arrays with names (*8w) are at w (full) levels.
!
!----------------------------------------------------------------------
! in wrf, 1 (smallest number) is the bottom level and nLevel (largest
! number) is the top level.  in your scheme, if 1 is at the top level,
! then you have to reverse the order in the k direction.
!
!         nLevel      -   half level (no data at this level)
!         nLevel    ----- full level
!         nLevel-1    -   half level
!         nLevel-1  ----- full level
!         .
!         .
!         .
!         1+2    -   half level
!         1+2  ----- full level
!         1+1    -   half level
!         1+1  ----- full level
!         1      -   half level
!         1    ----- full level
!
!======================================================================
! definitions
!-----------
! rho_d      dry density (kg/m^3)
! theta_m    moist potential temperature (k)
! qv         water vapor mixing ratio (kg/kg)
! qc         cloud water mixing ratio (kg/kg)
! qr         rain water mixing ratio (kg/kg)
! qi         cloud ice mixing ratio (kg/kg)
! qs         snow mixing ratio (kg/kg)
!-----------------------------------------------------------------
!-- dtime            tncell step (second)
!-- itimestep     number of tncell step (integer)   
!-- dx            horizontal space interval (m)
!-- nspecies       number of moisture species
!-- rr            dry air density (kg/m^3)
!
!-- rthcuten      theta tendency due to 
!                 cumulus scheme precipitation (k/s)
!-- rqvcuten      qv tendency due to 
!                 cumulus scheme precipitation (kg/kg/s)
!-- rqrcuten      qr tendency due to 
!                 cumulus scheme precipitation (kg/kg/s)
!-- rqccuten      qc tendency due to 
!                 cumulus scheme precipitation (kg/kg/s)
!-- rqscuten      qs tendency due to 
!                 cumulus scheme precipitation (kg/kg/s)
!-- rqicuten      qi tendency due to 
!                 cumulus scheme precipitation (kg/kg/s)
!
!-- rainc         accumulated total cumulus scheme precipitation (mm)
!-- raincv        cumulus scheme precipitation (mm)
!-- nca           counter of the cloud relaxation 
!                 time in kf cumulus scheme (integer)
!-- u_phy         u-velocity interpolated to theta points (m/s)
!-- v_phy         v-velocity interpolated to theta points (m/s)
!-- th_phy        potential temperature (k)
!-- t_phy         temperature (k)
!-- w             vertical velocity (m/s)
!-- moist         moisture array (4d - last index is species) (kg/kg)
!-- dz8w          dz between full levels (m)
!-- p8w           pressure at full levels (pa)  
!-- p_phy         pressure (pa)
!-- pi_phy        exner function (dimensionless)
!                 points (dimensionless)
!-- rthraten      radiative temp forcing for grell-devenyi scheme
!-- rthblten      pbl temp forcing for grell-devenyi scheme
!-- rqvblten      pbl moisture forcing for grell-devenyi scheme
!-- rthften
!-- rqvften
!-- mass_flux
!-- xf_ens
!-- pr_ens
!-- warm_rain
!-- cu_act_flag
!-- config_flags  
!-- w0avg         average vertical velocity, (for kf scheme) (m/s)
!-- rho_phy       density (kg/m^3)
!-- cldefi        precipitation efficiency (for bmj scheme) (dimensionless)
!-- stepcu        # of fundamental timesteps between convection calls
!-- xland         land-sea mask (1.0 for land; 2.0 for water)
!-- lowlyr        index of lowest model layer above the ground
!-- xlv0          latent heat of vaporization constant 
!                 used in temperature dependent formula (j/kg)
!-- xlv1          latent heat of vaporization constant 
!                 used in temperature dependent formula (j/kg/k)
!-- xls0          latent heat of sublimation constant 
!                 used in temperature dependent formula (j/kg)
!-- xls1          latent heat of sublimation constant
!                 used in temperature dependent formula (j/kg/k)
!-- rdry           gas constant for dry air ( 287. j/kg/k)
!-- rvap           gas constant for water vapor (461 j/k/kg)
!-- cp            specific heat at constant pressure (1004 j/k/kg)
!-- rvovrd        rvap divided by rdry (dimensionless)
!-- g             acceleration due to gravity (m/s^2)
!-- ep_1          constant for virtual temperature 
!                 (rvap/rdry - 1) (dimensionless)
!-- pi_phy        the exner function, (p/p0)**(r/cp) (none unit)
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
!-- 1           start index for i in memory
!-- ncell           end index for i in memory
!-- 1           start index for j in memory
!-- 1           end index for j in memory
!-- 1           start index for k in memory
!-- nLevel           end index for k in memory
!-- i_start       start indices for i in tile
!-- i_end         end indices for i in tile
!-- j_start       start indices for j in tile
!-- j_end         end indices for j in tile
!-- kts           start index for k in tile
!-- kte           end index for k in tile
!-- num_tiles     number of tiles
!-- hbot          index of lowest model layer with convection
!-- htop          index of highest model layer with convection
!-- lbot          index of lowest model layer with convection
!-- ltop          index of highest model layer with convection
!-- kpbl          layer index of the pbl
!
!======================================================================

! io
   integer(i4), intent(in)    :: ncell
   integer(i4), intent(in)    :: nLevel
   integer(i4), intent(in)    :: nspecies
   integer(i4), intent(in)    :: itimestep
   real(r8),    intent(in)    :: dtime
! local
   real(r8)      :: tmppratec(ncell,1)
   real(r8)      :: znu(ncell,nLevel)
   real(r8)      :: cudt, cu_dtime, curr_secs, cudtacttime  ! cumulus dtime in min
   real(r8)      :: rqvften (ncell,nLevel,1)
   real(r8)      :: rthften (ncell,nLevel,1)
   real(r8)      :: rqvblten(ncell,nLevel,1)
! gf internal
   integer(i4)   :: k22_shallow(ncell,1), kbcon_shallow(ncell,1), ktop_shallow(ncell,1)
   real(r8)      :: xmb_shallow(ncell,1), xmb_total(ncell,1), zzh(nLevel)
   integer(i4)   :: ktop_deep(ncell,1)
   real(r8)      :: edt_out(ncell,1), mass_flux(ncell,1), htop(ncell,1), hbot(ncell,1)
   real(r8)      :: apr_gr(1:ncell,1:1), apr_w (1:ncell,1:1),   apr_mc(1:ncell,1:1),   apr_st(1:ncell,1:1)
   real(r8)      :: apr_as(1:ncell,1:1), apr_capma(1:ncell,1:1),apr_capme(1:ncell,1:1),apr_capmi(1:ncell,1:1)
   integer(i4)   :: icell, ilev
    
    do ilev  = 1, nLevel-1
       znu(1:ncell,ilev) = pstate_wrf%p_phy(1:ncell,ilev,1)/psurf_wrf%psfc(1:ncell,1)
    end do

    cudt            = stepcu*dtime/60._r8
    cu_dtime        = stepcu*dtime
    cudtacttime     = zero
    curr_secs       = zero


if( (.not.unuse_cu) .and. (itimestep.eq.0 .or. itimestep .eq. 1 .or.  mod(itimestep,stepcu).eq.0) ) then

    select case(trim(wrfphys_cu_scheme))

    case('NTDKV381')
! This is Tiedtke-Bechtold scheme, 'new-tiedtke' in WRF's terminology
! check whether dyn forcing is necessary for ntdk: tropical rainfall sensitivity is small in APE
!
      rqvften(1:ncell,1:nLevel,1:1) = ptend_wrf%rqvblten(1:ncell,1:nLevel,1:1)+ptend_wrf%rqvdyten(1:ncell,1:nLevel,1:1)
      rthften(1:ncell,1:nLevel,1:1) = ptend_wrf%rthblten(1:ncell,1:nLevel,1:1)+ptend_wrf%rthdyten(1:ncell,1:nLevel,1:1)+ptend_wrf%rthraten(1:ncell,1:nLevel,1:1)
      !rqvften(1:ncell,1:nLevel,1:1) = ptend_wrf%rqvblten(1:ncell,1:nLevel,1:1)
      !rthften(1:ncell,1:nLevel,1:1) = ptend_wrf%rthblten(1:ncell,1:nLevel,1:1)+ptend_wrf%rthraten(1:ncell,1:nLevel,1:1)

      call cu_ntiedtke_v381(&
                dt       = dtime,     &
                itimestep= itimestep, &
                stepcu   = stepcu,    &
                raincv   = psurf_wrf%raincv( 1:ncell,1:1), & ! mm/model_step (fast physics)
                pratec   = tmppratec(        1:ncell,1:1), & ! mm/s
                qfx      = psurf_wrf%qfx(    1:ncell,1:1), &
                hfx      = psurf_wrf%hfx(    1:ncell,1:1), &
                u3d      = pstate_wrf%u_phy( 1:ncell,1:nLevel,1:1), & ! in
                v3d      = pstate_wrf%v_phy( 1:ncell,1:nLevel,1:1), & ! in
                w        = pstate_wrf%www(   1:ncell,1:nLevel,1:1), & ! in, not used
                t3d      = pstate_wrf%t_phy( 1:ncell,1:nLevel,1:1), & ! in
                qv3d     = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qv), & ! in
                qc3d     = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qc), & ! in
                qi3d     = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qi), & ! in
                pi3d     = pstate_wrf%pi_phy(1:ncell,1:nLevel,1:1),  & ! in
                rho3d    = pstate_wrf%rhom(  1:ncell,1:nLevel,1:1),  & ! in
                qvften   = rqvften          (1:ncell,1:nLevel,1:1),  & ! in, qv tend from other
                thften   = rthften          (1:ncell,1:nLevel,1:1),  & ! in, th tend from other
                dz8w     = pstate_wrf%dz8w(  1:ncell,1:nLevel,1:1) , & ! in
                pcps     = pstate_wrf%p_phy( 1:ncell,1:nLevel,1:1) , & ! in
                p8w      = pstate_wrf%p8w(   1:ncell,1:nLevel,1:1) , & ! in
                xland    = psurf_wrf%xland(  1:ncell,         1:1) , &
                cu_act_flag = cu_act_flag,                           & 
                dx       = pstate_wrf%dxmean(1:ncell),               &
                ids=1,ide=ncell, jds=1,jde=1, kds=1,kde=nLevel  ,    &
                ims=1,ime=ncell, jms=1,jme=1, kms=1,kme=nLevel  ,    &
                its=1,ite=ncell, jts=1,jte=1, kts=1,kte=nLevel-1,    &
                rthcuten = ptend_wrf%rthcuten(1:ncell,1:nLevel,1:1), & ! out
                rqvcuten = ptend_wrf%rqvcuten(1:ncell,1:nLevel,1:1), & ! out
                rqccuten = ptend_wrf%rqccuten(1:ncell,1:nLevel,1:1), & ! out
                rqicuten = ptend_wrf%rqicuten(1:ncell,1:nLevel,1:1), & ! out
                rucuten  = ptend_wrf%ruucuten(1:ncell,1:nLevel,1:1), & ! out
                rvcuten  = ptend_wrf%rvvcuten(1:ncell,1:nLevel,1:1), & ! out
                f_qv     = .true., &
                f_qc     = .true., &
                f_qr     = .false.,&
                f_qi     = .true., &
                f_qs     = .false.,&
                ndc  = nh_dynamics,&
                omega=pstate_wrf%omega(   1:ncell,1:nLevel,1:1), &
                cmfmc=pstate_wrf%cmfmc(   1:ncell,1:nLevel,1:1))

    case('GFV381')
      
       rqvften(1:ncell,1:nLevel,1:1) = ptend_wrf%rqvdyten(1:ncell,1:nLevel,1:1)
       rthften(1:ncell,1:nLevel,1:1) = ptend_wrf%rthdyten(1:ncell,1:nLevel,1:1)

       call GFDRV_v381(dt    = cu_dtime,                                 &     ! in
                  itimestep  = itimestep,                                &     ! in
                  dx         = pstate_wrf%dxmean(1:ncell),               &     ! in
                  rho        = pstate_wrf%rhom(  1:ncell,1:nLevel,1:1),  &     ! in
                  RAINCV     = psurf_wrf%raincv( 1:ncell,1:1),           &     ! mm/model_step (fast physics),
                  PRATEC     = tmppratec(        1:ncell,1:1),           &     ! mm/s
                  U          = pstate_wrf%u_phy( 1:ncell,1:nLevel,1:1),  &     ! in
                  V          = pstate_wrf%v_phy( 1:ncell,1:nLevel,1:1),  &     ! in
                  t          = pstate_wrf%t_phy( 1:ncell,1:nLevel,1:1),  &     ! in
                  W          = pstate_wrf%www(   1:ncell,1:nLevel,1:1),  &     ! in, double check, not used yet
                  q          = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qv), & ! in
                  p          = pstate_wrf%p_phy( 1:ncell,1:nLevel,1:1),  &     ! in
                  pi         = pstate_wrf%pi_phy(1:ncell,1:nLevel,1:1),  &     ! in
                  dz8w       = pstate_wrf%dz8w(  1:ncell,1:nLevel,1:1),  &     ! in
                  p8w        = pstate_wrf%p8w(   1:ncell,1:nLevel,1:1),  &     ! in
                  XLV        = xlv,           & 
                  CP         = cp,            &
                  G          = gravity,       &
                  r_v        = r_v,           &
                  htop       = htop(1:ncell,1:1), &
                  hbot       = hbot(1:ncell,1:1), &
                  ktop_deep  = ktop_deep(1:ncell,1:1),            & ! diag var, dum now
                  CU_ACT_FLAG= cu_act_flag(1:ncell,1:1),          & ! not active actually
                  warm_rain  = warm_rain,                         & ! not used actually
                  APR_GR     = apr_gr(1:ncell,1:1),               & !--
                  APR_W      = apr_w (1:ncell,1:1),               &
                  APR_MC     = apr_mc(1:ncell,1:1),               &
                  APR_ST     = apr_st(1:ncell,1:1),               &
                  APR_AS     = apr_as(1:ncell,1:1),               & ! seems unused
                  APR_CAPMA  = apr_capma(1:ncell,1:1),            &
                  APR_CAPME  = apr_capme(1:ncell,1:1),            &
                  APR_CAPMI  = apr_capmi(1:ncell,1:1),            & !--
                  MASS_FLUX  = mass_flux(1:ncell,1:1),            &
                  HT         = psurf_wrf%ht(1:ncell,1:1),         &
                  hfx        = psurf_wrf%hfx(1:ncell,1:1),        &
                  qfx        = psurf_wrf%qfx(1:ncell,1:1),        &
                  XLAND      = psurf_wrf%xland(1:ncell,1:1),      &
                  gsw        = psurf_wrf%gsw(1:ncell,1:1),        &
                  edt_out    = edt_out(1:ncell,1:1),              & ! diag var, dum now
                  !not use now
                  !GDC,
                  !GDC2,
                  !
                  kpbl           = pstate_wrf%kpbl(1:ncell,1:1),  & ! from pbl
                  k22_shallow    = k22_shallow(1:ncell,1:1),      & ! diag var, dum now
                  kbcon_shallow  = kbcon_shallow(1:ncell,1:1),    & ! diag var, dum now
                  ktop_shallow   = ktop_shallow(1:ncell,1:1),     & ! diag var, dum now
                  xmb_shallow    = xmb_shallow(1:ncell,1:1),      & ! diag var, dum now
                  ! not use
                  !cugd_tten,
                  !cugd_qvten,
                  !cugd_qcten,
                  !cugd_ttens,
                  !cugd_qvtens,
                  cugd_avedx  = 999, & ! not used
                  imomentum   = 999, & ! not used
                  ichoice     = 0  , & ! not used
                  ishallow_g3 = 1,   & !
                  ids = 1, ide= ncell, jds=1,jde=1, kds=1,kde=nLevel,     &
                  ims = 1, ime= ncell, jms=1,jme=1, kms=1,kme=nLevel,     &
                  its = 1, ite= ncell, jts=1,jte=1, kts=1,kte=nLevel-1,   &
                  periodic_x  = .false., &
                  periodic_y  = .false., &
                  RQVCUTEN    = ptend_wrf%rqvcuten(1:ncell,1:nLevel,1:1), & ! out
                  RQCCUTEN    = ptend_wrf%rqccuten(1:ncell,1:nLevel,1:1), & ! out
                  RQICUTEN    = ptend_wrf%rqicuten(1:ncell,1:nLevel,1:1), & ! out
                  RQVFTEN     = rqvften(1:ncell,1:nLevel,1:1), & ! in
                  RTHFTEN     = rthften(1:ncell,1:nLevel,1:1), & ! in
                  RTHCUTEN    = ptend_wrf%rthcuten(1:ncell,1:nLevel,1:1), & ! out
                  RTHRATEN    = ptend_wrf%rthraten(1:ncell,1:nLevel,1:1), & ! in
                  rqvblten    = ptend_wrf%rqvblten(1:ncell,1:nLevel,1:1), & ! in
                  rthblten    = ptend_wrf%rthblten(1:ncell,1:nLevel,1:1), &
                  F_QV=.true.,F_QC=.true.,F_QR=.false.,F_QI=.true.,F_QS=.false., &
                  ndc         = nh_dynamics                             , &
                  omega_hdc   = pstate_wrf%omega(  1:ncell,1:nLevel,1:1))
                  ! not used
                  ! Optional CAP suppress option
                  !,do_capsuppress,cap_suppress_loc)

    case default
         if(mpi_rank().eq.0) print*,"you must set wrfphys_cu_scheme, available options: TDK"
         call mpi_abort()
    end select

end if
!
! raincv is mm / step, so we should accumulate it each step, even cu is called
! several steps, implying no-call leads to the same mm/step for raincv
! raincv needs in restart
!
       psurf_wrf%rainc = psurf_wrf%rainc + psurf_wrf%raincv

      return
   end subroutine grist_wrf_cumulus_run
end module grist_wrf_cumulus_driver
