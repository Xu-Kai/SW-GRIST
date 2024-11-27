module module_ra_cam
  use module_ra_cam_support
#ifdef AMIPW_PHYSICS
  use grist_handle_error, only: endrun
#else
  use module_cam_support, only: endrun
#endif

  implicit none
! 
! a. slingo's data for cloud particle radiative properties (from 'a gcm
! parameterization for the shortwave properties of water clouds' jas
! vol. 46 may 1989 pp 1419-1427)
! 
   real(r8) abarl(4)         ! a coefficient for extinction optical depth
   real(r8) bbarl(4)         ! b coefficient for extinction optical depth
   real(r8) cbarl(4)         ! c coefficient for single scat albedo
   real(r8) dbarl(4)         ! d coefficient for single  scat albedo
   real(r8) ebarl(4)         ! e coefficient for asymmetry parameter
   real(r8) fbarl(4)         ! f coefficient for asymmetry parameter

   save abarl, bbarl, cbarl, dbarl, ebarl, fbarl

   data abarl/ 2.817e-02, 2.682e-02,2.264e-02,1.281e-02/
   data bbarl/ 1.305    , 1.346    ,1.454    ,1.641    /
   data cbarl/-5.62e-08 ,-6.94e-06 ,4.64e-04 ,0.201    /
   data dbarl/ 1.63e-07 , 2.35e-05 ,1.24e-03 ,7.56e-03 /
   data ebarl/ 0.829    , 0.794    ,0.754    ,0.826    /
   data fbarl/ 2.482e-03, 4.226e-03,6.560e-03,4.353e-03/

#if 0
! moved and changed to local variables into radcswmx for thread-safety, jm 20100217
   real(r8) abarli           ! a coefficient for current spectral band
   real(r8) bbarli           ! b coefficient for current spectral band
   real(r8) cbarli           ! c coefficient for current spectral band
   real(r8) dbarli           ! d coefficient for current spectral band
   real(r8) ebarli           ! e coefficient for current spectral band
   real(r8) fbarli           ! f coefficient for current spectral band
#endif
! 
! caution... a. slingo recommends no less than 4.0 micro-meters nor
! greater than 20 micro-meters
! 
! ice water coefficients (ebert and curry,1992, jgr, 97, 3831-3836)
! 
   real(r8) abari(4)         ! a coefficient for extinction optical depth
   real(r8) bbari(4)         ! b coefficient for extinction optical depth
   real(r8) cbari(4)         ! c coefficient for single scat albedo
   real(r8) dbari(4)         ! d coefficient for single scat albedo
   real(r8) ebari(4)         ! e coefficient for asymmetry parameter
   real(r8) fbari(4)         ! f coefficient for asymmetry parameter

   save abari, bbari, cbari, dbari, ebari, fbari

   data abari/ 3.448e-03, 3.448e-03,3.448e-03,3.448e-03/
   data bbari/ 2.431    , 2.431    ,2.431    ,2.431    /
   data cbari/ 1.00e-05 , 1.10e-04 ,1.861e-02,.46658   /
   data dbari/ 0.0      , 1.405e-05,8.328e-04,2.05e-05 /
   data ebari/ 0.7661   , 0.7730   ,0.794    ,0.9595   /
   data fbari/ 5.851e-04, 5.665e-04,7.267e-04,1.076e-04/

#if 0
! moved and changed to local variables into radcswmx for thread-safety, jm 20100217
   real(r8) abarii           ! a coefficient for current spectral band
   real(r8) bbarii           ! b coefficient for current spectral band
   real(r8) cbarii           ! c coefficient for current spectral band
   real(r8) dbarii           ! d coefficient for current spectral band
   real(r8) ebarii           ! e coefficient for current spectral band
   real(r8) fbarii           ! f coefficient for current spectral band
#endif
! 
   real(r8) delta            ! pressure (in atm) for stratos. h2o limit
   real(r8) o2mmr            ! o2 mass mixing ratio:

   save delta, o2mmr

!
! update to h2o near-ir: delta optimized for hitran 2k and ckd 2.4
!
   data delta / 0.0014257179260883 /
!
! end update
!
   data o2mmr / .23143 /

! next series depends on spectral interval
! 
   real(r8) frcsol(nspint)   ! fraction of solar flux in spectral interval
   real(r8) wavmin(nspint)   ! min wavelength (micro-meters) of interval
   real(r8) wavmax(nspint)   ! max wavelength (micro-meters) of interval
   real(r8) raytau(nspint)   ! rayleigh scattering optical depth
   real(r8) abh2o(nspint)    ! absorption coefficiant for h2o (cm2/g)
   real(r8) abo3 (nspint)    ! absorption coefficiant for o3  (cm2/g)
   real(r8) abco2(nspint)    ! absorption coefficiant for co2 (cm2/g)
   real(r8) abo2 (nspint)    ! absorption coefficiant for o2  (cm2/g)
   real(r8) ph2o(nspint)     ! weight of h2o in spectral interval
   real(r8) pco2(nspint)     ! weight of co2 in spectral interval
   real(r8) po2 (nspint)     ! weight of o2  in spectral interval
   real(r8) nirwgt(nspint)   ! spectral weights to simulate nimbus-7 filter
   save frcsol ,wavmin ,wavmax ,raytau ,abh2o ,abo3 , &
        abco2  ,abo2   ,ph2o   ,pco2   ,po2   ,nirwgt

   data frcsol / .001488, .001389, .001290, .001686, .002877, &
                 .003869, .026336, .360739, .065392, .526861, &
                 .526861, .526861, .526861, .526861, .526861, &
                 .526861, .006239, .001834, .001834/
! 
! weight for 0.64 - 0.7 microns  appropriate to clear skies over oceans
! 
   data nirwgt /  0.0,   0.0,   0.0,      0.0,   0.0, &
                  0.0,   0.0,   0.0, 0.320518,   1.0,  1.0, &
                  1.0,   1.0,   1.0,      1.0,   1.0, &
                  1.0,   1.0,   1.0 /

   data wavmin / .200,  .245,  .265,  .275,  .285, &
                 .295,  .305,  .350,  .640,  .700,  .701, &
                 .701,  .701,  .701,  .702,  .702, &
                 2.630, 4.160, 4.160/

   data wavmax / .245,  .265,  .275,  .285,  .295, &
                 .305,  .350,  .640,  .700, 5.000, 5.000, &
                 5.000, 5.000, 5.000, 5.000, 5.000, &
                 2.860, 4.550, 4.550/

!
! update to h2o near-ir: rayleigh scattering optimized for hitran 2k & ckd 2.4
!
   real(r8) v_raytau_35
   real(r8) v_raytau_64
   real(r8) v_abo3_35
   real(r8) v_abo3_64
   parameter( &
        v_raytau_35 = 0.155208, &
        v_raytau_64 = 0.0392, &
        v_abo3_35 = 2.4058030e+01, &  
        v_abo3_64 = 2.210e+01 &
        )

   data raytau / 4.020, 2.180, 1.700, 1.450, 1.250, &
                  1.085, 0.730, v_raytau_35, v_raytau_64, &
                  0.02899756, 0.01356763, 0.00537341, &
                  0.00228515, 0.00105028, 0.00046631, &
                  0.00025734, &
                 .0001, .0001, .0001/
!
! end update
!

! 
! absorption coefficients
! 
!
! update to h2o near-ir: abh2o optimized for hitran 2k and ckd 2.4
!
   data abh2o /    .000,     .000,    .000,    .000,    .000, &
                   .000,     .000,    .000,    .000,    &
                   0.00256608,  0.06310504,   0.42287445, 2.45397941, &
                  11.20070807, 47.66091389, 240.19010243, &
                   .000,    .000,    .000/
!
! end update
!

   data abo3  /5.370e+04, 13.080e+04,  9.292e+04, 4.530e+04, 1.616e+04, &
               4.441e+03,  1.775e+02, v_abo3_35, v_abo3_64,      .000, &
               .000,   .000    ,   .000   ,   .000   ,      .000, &
               .000,   .000    ,   .000   ,   .000    /

   data abco2  /   .000,     .000,    .000,    .000,    .000, &
                   .000,     .000,    .000,    .000,    .000, &
                   .000,     .000,    .000,    .000,    .000, &
                   .000,     .094,    .196,   1.963/

   data abo2  /    .000,     .000,    .000,    .000,    .000, &
                   .000,     .000,    .000,1.11e-05,6.69e-05, &
                   .000,     .000,    .000,    .000,    .000, &  
                   .000,     .000,    .000,    .000/
! 
! spectral interval weights
! 
   data ph2o  /    .000,     .000,    .000,    .000,    .000, &
        .000,     .000,    .000,    .000,    .505,     &
        .210,     .120,    .070,    .048,    .029,     &
        .018,     .000,    .000,    .000/

   data pco2  /    .000,     .000,    .000,    .000,    .000, &
        .000,     .000,    .000,    .000,    .000,     &
        .000,     .000,    .000,    .000,    .000,     &
        .000,    1.000,    .640,    .360/

   data po2   /    .000,     .000,    .000,    .000,    .000, &
        .000,     .000,    .000,   1.000,   1.000,     &
        .000,     .000,    .000,    .000,    .000,     &
        .000,     .000,    .000,    .000/

   real(r8) amo                 ! molecular weight of ozone (g/mol)
   save     amo

   data amo   /  48.0000   /

contains
subroutine camrad(rthratenlw,rthratensw,                           &
                     dolw,dosw,                                    &
                     swupt,swuptc,swdnt,swdntc,                    &
                     lwupt,lwuptc,lwdnt,lwdntc,                    &
                     swupb,swupbc,swdnb,swdnbc,                    &
                     lwupb,lwupbc,lwdnb,lwdnbc,                    &
                     swcf,lwcf,olr,cemiss,taucldc,taucldi,coszr,   &
                     gsw,glw,xlat,xlong,                           &
                     albedo,t_phy,tsk,emiss,                       &
                     qv3d,qc3d,qr3d,qi3d,qs3d,qg3d,                &
                     alswvisdir,alswvisdif,                        & !ssib  
                     alswnirdir,alswnirdif,                        & !ssib 
                     swvisdir,swvisdif,                            & !ssib 
                     swnirdir,swnirdif,                            & !ssib 
                     sf_surface_physics,                           & !ssib 
                     f_qv,f_qc,f_qr,f_qi,f_qs,f_qg,                &
                     f_ice_phy,f_rain_phy,                         &
                     p_phy,p8w,z,pi_phy,rho_phy,dz8w,               &
                     cldfra,xland,xice,snow,                        &
                     ozmixm,pin0,levsiz,num_months,                 &
                     m_psp,m_psn,aerosolcp,aerosolcn,m_hybi0,       &
                     cam_abs_dim1, cam_abs_dim2,                    &
                     paerlev,naer_c,                                &
                     gmt,julday,julian,yr,dt,xtime,declin,solcon,         &
                     radt,degrad,n_cldadv,                                  &
                     abstot_3d, absnxt_3d, emstot_3d,              &
                     doabsems,                                     &
                     ids,ide, jds,jde, kds,kde,                    &
                     ims,ime, jms,jme, kms,kme,                    &
                     its,ite, jts,jte, kts,kte                     )
#ifndef AMIPW_PHYSICS
   use module_wrf_error
   use module_state_description, only : ssibscheme              !ssib
#endif
!------------------------------------------------------------------
   implicit none
!------------------------------------------------------------------

   integer,    intent(in   ) ::        ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       its,ite, jts,jte, kts,kte
   logical,    intent(in   ) ::        f_qv,f_qc,f_qr,f_qi,f_qs,f_qg
   logical,    intent(inout) ::        doabsems
   logical,    intent(in   ) ::        dolw,dosw

   integer,    intent(in  )  ::        n_cldadv
   integer,    intent(in  )  ::        julday
   real,       intent(in  )  ::        julian
   integer,    intent(in  )  ::        yr
   real,       intent(in  )  ::        dt
   integer,      intent(in   )    ::   levsiz, num_months
   integer,      intent(in   )    ::   paerlev, naer_c
   integer,      intent(in   )    ::   cam_abs_dim1, cam_abs_dim2


   real, intent(in    )      ::        radt,degrad,             &
                                       xtime,declin,solcon,gmt
!
!
   real, dimension( ims:ime, kms:kme, jms:jme ),                  &
         intent(in    ) ::                                   p_phy, &
                                                           p8w, &
                                                             z, &
                                                            pi_phy, &
                                                           rho_phy, &
                                                              dz8w, &
                                                             t_phy, &
                                                            qv3d, &
                                                            qc3d, &
                                                            qr3d, &
                                                            qi3d, &
                                                            qs3d, &
                                                            qg3d, &
                                                        cldfra

   real, dimension( ims:ime, kms:kme, jms:jme ),                  &
         intent(inout)  ::                              rthratenlw, &
                                                        rthratensw
!
   real, dimension( ims:ime, jms:jme ),                           &
         intent(in   )  ::                                  xlat, &
                                                           xlong, &
                                                           xland, &
                                                           xice, &
                                                           snow, &
                                                           emiss, &
                                                             tsk, &
                                                             albedo

   real,  dimension( ims:ime, levsiz, jms:jme, num_months ),      &
          intent(in   ) ::                                  ozmixm

   real,  dimension(levsiz), intent(in )  ::                   pin0

   real,  dimension(ims:ime,jms:jme), intent(in )  ::      m_psp,m_psn
   real,  dimension(paerlev), intent(in)             ::      m_hybi0
   real,  dimension( ims:ime, paerlev, jms:jme, naer_c ),      &
          intent(in   ) ::                    aerosolcp, aerosolcn

!
   real, dimension( ims:ime, jms:jme ),                           &
         intent(inout)  ::                                   gsw, glw

!---------ssib variables (fds 06/2010)----------------
   real, dimension( ims:ime, jms:jme ),                           &
         intent(in)     ::                            alswvisdir, &
                                                      alswvisdif, &
                                                      alswnirdir, &
                                                      alswnirdif

   real, dimension( ims:ime, jms:jme ),                           &
         intent(out)    ::                              swvisdir, &
                                                        swvisdif, &
                                                        swnirdir, &
                                                        swnirdif
   integer, intent(in) :: sf_surface_physics
!--------------------------------------
! saving arrays for doabsems reduction of radiation calcs

   real, dimension( ims:ime, kms:kme, cam_abs_dim2 , jms:jme ),           &
         intent(inout)  ::                                  abstot_3d
   real, dimension( ims:ime, kms:kme, cam_abs_dim1 , jms:jme ),           &
         intent(inout)  ::                                  absnxt_3d
   real, dimension( ims:ime, kms:kme, jms:jme ),           &
         intent(inout)  ::                                  emstot_3d


! added outputs of total and clearsky fluxes etc
! note that k=1 refers to the half level below the model lowest level (sfc)
!           k=kme refers to the half level above the model highest level (toa)
!
!   real, dimension( ims:ime, kms:kme, jms:jme ),                &
!         intent(inout)  ::                                  swup, &
!                                                       swupclear, &
!                                                            swdn, &
!                                                       swdnclear, &
!                                                            lwup, &
!                                                       lwupclear, &
!                                                            lwdn, &
!                                                       lwdnclear

   real, dimension( ims:ime, jms:jme ), optional, intent(inout) ::&
                    swupt,swuptc,swdnt,swdntc,                    &
                    lwupt,lwuptc,lwdnt,lwdntc,                    &
                    swupb,swupbc,swdnb,swdnbc,                    &
                    lwupb,lwupbc,lwdnb,lwdnbc

   real, dimension( ims:ime, jms:jme ),                           &
         intent(inout)  ::                                  swcf, &
                                                            lwcf, &
                                                             olr, &
                                                            coszr    
   real, dimension( ims:ime, kms:kme, jms:jme )                 , &
         intent(out   )  ::                               cemiss, &        ! cloud emissivity for isccp
                                                         taucldc, &        ! cloud water optical depth for isccp
                                                         taucldi           ! cloud ice optical depth for isccp
!
!
   real, dimension( ims:ime, kms:kme, jms:jme ),                     &
         intent(in   ) ::                                            &
                                                          f_ice_phy, &
                                                         f_rain_phy


! local variables
 
   integer :: lchnk, ncol, pcols, pver, pverp, pverr, pverrp
   integer :: pcnst, pnats, ppcnst, i, j, k, ii, kk, kk1, m, n
   integer :: begchunk, endchunk
   integer :: nyrm, nyrp
   real(r8) doymodel, doydatam, doydatap, deltat, fact1, fact2

   real :: xt24, tloctm, hrang, xxlat, oldxt24
 
   real(r8), dimension( 1:ite-its+1 ) :: coszrs, landfrac, landm, snowh, icefrac, lwups
   real(r8), dimension( 1:ite-its+1 ) :: asdir, asdif, aldir, aldif, ps
   real(r8), dimension( 1:ite-its+1, 1:kte-kts+1 ) :: cld, pmid, lnpmid, pdel, zm, t
   real(r8), dimension( 1:ite-its+1, 1:kte-kts+2 ) ::  pint, lnpint
   real(r8), dimension( 1:ite-its+1, 1:kte-kts+1, n_cldadv) :: q
!   real(r8), dimension( 1:kte-kts+1 ) :: hypm       ! reference pressures at midpoints
!   real(r8), dimension( 1:kte-kts+2 ) :: hypi       ! reference pressures at interfaces
    real(r8), dimension(  1:ite-its+1, 1:kte-kts+1 ) :: cicewp      ! in-cloud cloud ice water path
    real(r8), dimension(  1:ite-its+1, 1:kte-kts+1 ) :: cliqwp      ! in-cloud cloud liquid water path
    real(r8), dimension(  1:ite-its+1, 0:kte-kts+1 ) :: tauxcl      ! cloud water optical depth
    real(r8), dimension(  1:ite-its+1, 0:kte-kts+1 ) :: tauxci      ! cloud ice optical depth
    real(r8), dimension(  1:ite-its+1, 1:kte-kts+1 ) :: emis        ! cloud emissivity
    real(r8), dimension(  1:ite-its+1, 1:kte-kts+1 ) :: rel         ! effective drop radius (microns)
    real(r8), dimension(  1:ite-its+1, 1:kte-kts+1 ) :: rei         ! ice effective drop size (microns)
    real(r8), dimension(  1:ite-its+1, 1:kte-kts+2 ) :: pmxrgn      ! maximum values of pressure for each
    integer , dimension(  1:ite-its+1 ) :: nmxrgn               ! number of maximally overlapped regions

   real(r8), dimension(  1:ite-its+1 ) :: fsns          ! surface absorbed solar flux
   real(r8), dimension(  1:ite-its+1 ) :: fsnt          ! net column abs solar flux at model top
   real(r8), dimension(  1:ite-its+1 ) :: flns          ! srf longwave cooling (up-down) flux
   real(r8), dimension(  1:ite-its+1 ) :: flnt          ! net outgoing lw flux at model top
! added outputs of total and clearsky fluxes etc
   real(r8), dimension(  1:ite-its+1, 1:kte-kts+2 )  :: fsup        ! upward total sky solar
   real(r8), dimension(  1:ite-its+1, 1:kte-kts+2 )  :: fsupc       ! upward clear sky solar
   real(r8), dimension(  1:ite-its+1, 1:kte-kts+2 )  :: fsdn        ! downward total sky solar
   real(r8), dimension(  1:ite-its+1, 1:kte-kts+2 )  :: fsdnc       ! downward clear sky solar
   real(r8), dimension(  1:ite-its+1, 1:kte-kts+2 )  :: flup        ! upward total sky longwave
   real(r8), dimension(  1:ite-its+1, 1:kte-kts+2 )  :: flupc       ! upward clear sky longwave
   real(r8), dimension(  1:ite-its+1, 1:kte-kts+2 )  :: fldn        ! downward total sky longwave
   real(r8), dimension(  1:ite-its+1, 1:kte-kts+2 )  :: fldnc       ! downward clear sky longwave
   real(r8), dimension(  1:ite-its+1 ) :: swcftoa                 ! top of the atmosphere solar cloud forcing
   real(r8), dimension(  1:ite-its+1 ) :: lwcftoa                 ! top of the atmosphere longwave cloud forcing
   real(r8), dimension(  1:ite-its+1 ) :: olrtoa                  ! top of the atmosphere outgoing longwave 
!
   real(r8), dimension(  1:ite-its+1 ) :: sols          ! downward solar rad onto surface (sw direct)
   real(r8), dimension(  1:ite-its+1 ) :: soll          ! downward solar rad onto surface (lw direct)
   real(r8), dimension(  1:ite-its+1 ) :: solsd         ! downward solar rad onto surface (sw diffuse)
   real(r8), dimension(  1:ite-its+1 ) :: solld         ! downward solar rad onto surface (lw diffuse)
   real(r8), dimension(  1:ite-its+1, 1:kte-kts+1 ) :: qrs      ! solar heating rate
   real(r8), dimension(  1:ite-its+1 ) :: fsds          ! flux shortwave downwelling surface
   real(r8), dimension(  1:ite-its+1, 1:kte-kts+1 ) :: qrl      ! longwave cooling rate
   real(r8), dimension(  1:ite-its+1 ) :: flwds          ! surface down longwave flux
   real(r8), dimension(  1:ite-its+1, levsiz, num_months ) :: ozmixmj        ! monthly ozone mixing ratio
   real(r8), dimension(  1:ite-its+1, levsiz ) :: ozmix          ! ozone mixing ratio (time interpolated)
   real(r8), dimension(levsiz)         :: pin            ! ozone pressure level
   real(r8), dimension(1:ite-its+1)    :: m_psjp,m_psjn          ! match surface pressure
   real(r8), dimension(  1:ite-its+1, paerlev, naer_c ) :: aerosoljp        ! monthly aerosol concentrations
   real(r8), dimension(  1:ite-its+1, paerlev, naer_c ) :: aerosoljn        ! monthly aerosol concentrations
   real(r8), dimension(paerlev)                           :: m_hybi
   real(r8), dimension(1:ite-its+1 )          :: clat           ! latitude in radians for columns
   real(r8), dimension(its:ite,kts:kte+1,kts:kte+1) :: abstot ! total absorptivity
   real(r8), dimension(its:ite,kts:kte,4)           :: absnxt ! total nearest layer absorptivity
   real(r8), dimension(its:ite,kts:kte+1)           :: emstot ! total emissivity
   character(len=256) :: msgstr

#if !defined(mac_kludge)
   lchnk = 1
   begchunk = ims
   endchunk = ime
   ncol = ite - its + 1
   pcols= ite - its + 1
   pver = kte - kts + 1
   pverp= pver + 1
   pverr = kte - kts + 1
   pverrp= pverr + 1
! number of advected constituents and non-advected constituents (including water vapor)
   ppcnst = n_cldadv
! number of non-advected constituents
   pnats = 0
   pcnst = ppcnst-pnats

! check the # species defined for the input climatology and naer

!  if(naer_c.ne.naer) then
!            write( wrf_err_message , * ) 'naer_c ne naer ', naer_c, naer
   if(naer_c.ne.naer_all) then
#ifdef AMIPW_PHYSICS
             call endrun('naer_c-1 ne naer_all ')
#else
             write( wrf_err_message , * ) 'naer_c-1 ne naer_all ', naer_c, naer_all
             call wrf_error_fatal ( wrf_err_message )
#endif
   endif 

! update co2 volume mixing ratio (co2vmr)
  
! determine time interpolation factors, check sanity
! of interpolation factors to within 32-bit roundoff
! assume that day of year is 1 for all input data
!
   nyrm     = yr - yrdata(1) + 1
   nyrp     = nyrm + 1
   doymodel = yr*365.    + julian
   doydatam = yrdata(nyrm)*365. + 1.
   doydatap = yrdata(nyrp)*365. + 1.
   deltat   = doydatap - doydatam
   fact1    = (doydatap - doymodel)/deltat
   fact2    = (doymodel - doydatam)/deltat
   co2vmr = (co2(nyrm)*fact1 + co2(nyrp)*fact2)*1.e-06

   co2mmr=co2vmr*mwco2/mwdry
!
!===================================================
! radiation computations
!===================================================

      do k=1,levsiz
      pin(k)=pin0(k)
      enddo

      do k=1,paerlev
      m_hybi(k)=m_hybi0(k)
      enddo

! check for uninitialized arrays
      if(abstot_3d(its,kts,kts,jts) .eq. 0.0 .and. .not.doabsems .and. dolw)then
#ifdef AMIPW_PHYSICS
        call endrun( 'camrad lw: caution: re-calculating abstot, absnxt, emstot on restart')
#else
        call wrf_debug(0, 'camrad lw: caution: re-calculating abstot, absnxt, emstot on restart')
#endif
        doabsems = .true.
      endif

   do j =jts,jte

!
! cosine solar zenith angle for current time step
!

!  call zenith (calday, clat, clon, coszrs, ncol)

      do i = its,ite
      ii = i - its + 1
      ! xt24 is the fractional part of simulation days plus half of radt expressed in 
      ! units of minutes
      ! julian is in days
      ! radt is in minutes
      xt24=mod(xtime+radt*0.5,1440.)
      tloctm=gmt+xt24/60.+xlong(i,j)/15.
      hrang=15.*(tloctm-12.)*degrad
      xxlat=xlat(i,j)*degrad
      clat(ii)=xxlat
      coszrs(ii)=sin(xxlat)*sin(declin)+cos(xxlat)*cos(declin)*cos(hrang)
      enddo

! moist variables

      do k = kts,kte
      kk = kte - k + kts 
      do i = its,ite
      ii = i - its + 1
!    convert to specific humidity
      q(ii,kk,1) = max(1.e-10,qv3d(i,k,j)/(1.+qv3d(i,k,j)))
     if ( f_qi .and. f_qc .and. f_qs ) then
      q(ii,kk,ixcldliq) = max(0.,qc3d(i,k,j)/(1.+qv3d(i,k,j)))
      q(ii,kk,ixcldice) = max(0.,(qi3d(i,k,j)+qs3d(i,k,j))/(1.+qv3d(i,k,j)))
     else if ( f_qc .and. f_qr ) then
! warm rain or simple ice
      q(ii,kk,ixcldliq) = 0.
      q(ii,kk,ixcldice) = 0.
      if(t_phy(i,k,j).gt.273.15)q(ii,kk,ixcldliq) = max(0.,qc3d(i,k,j)/(1.+qv3d(i,k,j)))
      if(t_phy(i,k,j).le.273.15)q(ii,kk,ixcldice) = max(0.,qc3d(i,k,j)/(1.+qv3d(i,k,j)))
     else if ( f_qc .and. f_qs ) then
! for ferrier (note that currently ferrier has qi, so this section will not be used)
      q(ii,kk,ixcldice) = max(0.,qc3d(i,k,j)/(1.+qv3d(i,k,j))*f_ice_phy(i,k,j))
      q(ii,kk,ixcldliq) = max(0.,qc3d(i,k,j)/(1.+qv3d(i,k,j))*(1.-f_ice_phy(i,k,j))*(1.-f_rain_phy(i,k,j)))
     else
      q(ii,kk,ixcldliq) = 0.
      q(ii,kk,ixcldice) = 0.
     endif
      cld(ii,kk) = cldfra(i,k,j)
      enddo
      enddo

      do i = its,ite
      ii = i - its + 1
      landfrac(ii) = 2.-xland(i,j)
      landm(ii) = landfrac(ii)
      snowh(ii) = 0.001*snow(i,j)
      icefrac(ii) = xice(i,j)
      enddo

      do m=1,num_months-1
      do k=1,levsiz
      do i = its,ite
      ii = i - its + 1
      ozmixmj(ii,k,m) = ozmixm(i,k,j,m+1)
      enddo
      enddo
      enddo

      do i = its,ite
      ii = i - its + 1
      m_psjp(ii) = m_psp(i,j)
      m_psjn(ii) = m_psn(i,j)
      enddo

      do n=1,naer_c
      do k=1,paerlev
      do i = its,ite
      ii = i - its + 1
      aerosoljp(ii,k,n) = aerosolcp(i,k,j,n)
      aerosoljn(ii,k,n) = aerosolcn(i,k,j,n)
      enddo
      enddo
      enddo

!
! complete radiation calculations
!
      do i = its,ite
      ii = i - its + 1
      lwups(ii) = stebol*emiss(i,j)*tsk(i,j)**4
      enddo

      do k = kts,kte+1
      kk = kte - k + kts + 1
      do i = its,ite
      ii = i - its + 1
      pint(ii,kk) = p8w(i,k,j)
      if(k.eq.kts)ps(ii)=pint(ii,kk)
      lnpint(ii,kk) = log(pint(ii,kk))
      enddo
      enddo

      if(.not.doabsems .and. dolw)then
!      do kk = kts,kte+1
      do kk = 1,cam_abs_dim2
        do kk1 = kts,kte+1
          do i = its,ite
            abstot(i,kk1,kk) = abstot_3d(i,kk1,kk,j)
          enddo
        enddo
      enddo
!      do kk = 1,4
      do kk = 1,cam_abs_dim1
        do kk1 = kts,kte
          do i = its,ite
            absnxt(i,kk1,kk) = absnxt_3d(i,kk1,kk,j)
          enddo
        enddo
      enddo
      do kk = kts,kte+1
          do i = its,ite
            emstot(i,kk) = emstot_3d(i,kk,j)
          enddo
      enddo
      endif

      do k = kts,kte
      kk = kte - k + kts 
      do i = its,ite
      ii = i - its + 1
      pmid(ii,kk) = p_phy(i,k,j)
      lnpmid(ii,kk) = log(pmid(ii,kk))
      lnpint(ii,kk) = log(pint(ii,kk))
      pdel(ii,kk) = pint(ii,kk+1) - pint(ii,kk)
      t(ii,kk) = t_phy(i,k,j)
      zm(ii,kk) = z(i,k,j)
      enddo
      enddo


! compute cloud water/ice paths and optical properties for input to radiation

      call param_cldoptics_calc(ncol, pcols, pver, pverp, pverr, pverrp, ppcnst, q, cld, landfrac, landm,icefrac, &
                                pdel, t, ps, pmid, pint, cicewp, cliqwp, emis, rel, rei, pmxrgn, nmxrgn, snowh)

!-----fds (06/2010)----------------------------
   select case(sf_surface_physics)
#ifdef AMIPW_PHYSICS
   case (8) ! use 8 for using this part
#else
   case (ssibscheme)
#endif
   if (xtime .gt. 1.0) then
!      call wrf_message("using ssib albedoes for land points")
      do i = its,ite
         ii = i - its + 1
         if (xland(i,j).lt.1.5) then   !land points only
           asdir(ii) = alswvisdir(i,j) ! ssib visdir albedo
           asdif(ii) = alswvisdif(i,j) ! ssib visdif albedo
           aldir(ii) = alswnirdir(i,j) ! ssib nirdir albedo
           aldif(ii) = alswnirdif(i,j) ! ssib nirdif albedo
         else
           asdir(ii) = albedo(i,j)
           asdif(ii) = albedo(i,j)
           aldir(ii) = albedo(i,j)
           aldif(ii) = albedo(i,j)
         endif
      enddo
   else
      do i = its,ite
         ii = i - its + 1
         asdir(ii) = albedo(i,j)
         asdif(ii) = albedo(i,j)
         aldir(ii) = albedo(i,j)
         aldif(ii) = albedo(i,j)
      enddo
   endif
   case default
      do i = its,ite
      ii = i - its + 1
! use same albedo for direct and diffuse
! change this when separate values are provided
      asdir(ii) = albedo(i,j)
      asdif(ii) = albedo(i,j)
      aldir(ii) = albedo(i,j)
      aldif(ii) = albedo(i,j)
      enddo
   end select
!-----------------------------------------------

! wrf allocate space here (not needed if oznini is called)
!  allocate (ozmix(pcols,levsiz,begchunk:endchunk)) ! this line from oznini.f90

      call radctl (j,lchnk, ncol, pcols, pver, pverp, pverr, pverrp, ppcnst, pcnst, lwups, emis, pmid,             &
                   pint, lnpmid, lnpint, pdel, t, q,   &
                   cld, cicewp, cliqwp, tauxcl, tauxci, coszrs, clat, asdir, asdif,               &
                   aldir, aldif, solcon, gmt,julday,julian,dt,xtime,   &
                   pin, ozmixmj, ozmix, levsiz, num_months,  & 
                   m_psjp,m_psjn, aerosoljp, aerosoljn,  m_hybi, paerlev, naer_c, pmxrgn, nmxrgn, &
                   dolw, dosw, doabsems, abstot, absnxt, emstot, &
                   fsup, fsupc, fsdn, fsdnc, flup, flupc, fldn, fldnc, swcftoa, lwcftoa, olrtoa,  &
                   fsns, fsnt    ,flns    ,flnt    , &
                   qrs, qrl, flwds, rel, rei,                       &
                   sols, soll, solsd, solld,                  &
                   landfrac, zm, fsds)

      do k = kts,kte
      kk = kte - k + kts 
      do i = its,ite
      ii = i - its + 1
      if(dolw)rthratenlw(i,k,j) = 1.e4*qrl(ii,kk)/(cpair*pi_phy(i,k,j))
      if(dosw)rthratensw(i,k,j) = 1.e4*qrs(ii,kk)/(cpair*pi_phy(i,k,j))
      cemiss(i,k,j)     = emis(ii,kk)
      taucldc(i,k,j)    = tauxcl(ii,kk)
      taucldi(i,k,j)    = tauxci(ii,kk)
      enddo
      enddo

      if(doabsems .and. dolw)then
!      do kk = kts,kte+1
      do kk = 1,cam_abs_dim2
        do kk1 = kts,kte+1
          do i = its,ite
            abstot_3d(i,kk1,kk,j) = abstot(i,kk1,kk)
          enddo
        enddo
      enddo
!      do kk = 1,4
      do kk = 1,cam_abs_dim1
        do kk1 = kts,kte
          do i = its,ite
            absnxt_3d(i,kk1,kk,j) = absnxt(i,kk1,kk)
          enddo
        enddo
      enddo
      do kk = kts,kte+1
          do i = its,ite
            emstot_3d(i,kk,j) = emstot(i,kk)
          enddo
      enddo
      endif

      if(present(swupt))then
      if(dosw)then
! added shortwave and longwave upward/downward total and clear sky fluxes
      do k = kts,kte+1
      kk = kte +1 - k + kts
      do i = its,ite
      ii = i - its + 1
!      swup(i,k,j)      = fsup(ii,kk)
!      swupclear(i,k,j) = fsupc(ii,kk)
!      swdn(i,k,j)      = fsdn(ii,kk)
!      swdnclear(i,k,j) = fsdnc(ii,kk)
       if(k.eq.kte+1)then
         swupt(i,j)     = fsup(ii,kk)
         swuptc(i,j)    = fsupc(ii,kk)
         swdnt(i,j)     = fsdn(ii,kk)
         swdntc(i,j)    = fsdnc(ii,kk)
       endif
       if(k.eq.kts)then
         swupb(i,j)     = fsup(ii,kk)
         swupbc(i,j)    = fsupc(ii,kk)
         swdnb(i,j)     = fsdn(ii,kk)
         swdnbc(i,j)    = fsdnc(ii,kk)
       endif
!            if(i.eq.30.and.j.eq.30) then
!            print 1234, 'short ', i,ii,k,kk,fsup(ii,kk),fsupc(ii,kk),fsdn(ii,kk),fsdnc(ii,kk)
!            1234 format (a6,4i4,4f10.3)
!            endif
     enddo
      enddo
      endif
      if(dolw)then
! added shortwave and longwave upward/downward total and clear sky fluxes
      do k = kts,kte+1
      kk = kte +1 - k + kts
      do i = its,ite
      ii = i - its + 1
!      lwup(i,k,j)      = flup(ii,kk)
!      lwupclear(i,k,j) = flupc(ii,kk)
!      lwdn(i,k,j)      = fldn(ii,kk)
!      lwdnclear(i,k,j) = fldnc(ii,kk)
       if(k.eq.kte+1)then
         lwupt(i,j)     = flup(ii,kk)
         lwuptc(i,j)    = flupc(ii,kk)
         lwdnt(i,j)     = fldn(ii,kk)
         lwdntc(i,j)    = fldnc(ii,kk)
       endif
       if(k.eq.kts)then
         lwupb(i,j)     = flup(ii,kk)
         lwupbc(i,j)    = flupc(ii,kk)
         lwdnb(i,j)     = fldn(ii,kk)
         lwdnbc(i,j)    = fldnc(ii,kk)
       endif
!            if(i.eq.30.and.j.eq.30) then
!            print 1234, 'long  ', i,ii,k,kk,flup(ii,kk),flupc(ii,kk),fldn(ii,kk),fldnc(ii,kk)
!            1234 format (a6,4i4,4f10.3)
!            endif
      enddo
      enddo
      endif
      endif

      do i = its,ite
      ii = i - its + 1
! added shortwave and longwave cloud forcing at toa and surface
      if(dolw)then
        glw(i,j) = flwds(ii)
        lwcf(i,j) = lwcftoa(ii)
        olr(i,j)  = olrtoa(ii)
      endif
      if(dosw)then
        gsw(i,j) = fsns(ii)
        swcf(i,j) = swcftoa(ii)
        coszr(i,j) = coszrs(ii)
      endif
      enddo
!-------fds (06/2010)---------
   select case(sf_surface_physics)
#ifdef AMIPW_PHYSICS
   case (8)
#else
   case (ssibscheme)
#endif
!   call wrf_message("cam using ssib albedo2")
      if(dosw)then
      do i = its,ite
        ii = i - its + 1
        swvisdir(i,j) = sols(ii)  !ssib 
        swvisdif(i,j) = solsd(ii) !ssib 
        swnirdir(i,j) = soll(ii)  !ssib 
        swnirdif(i,j) = solld(ii) !ssib 
      enddo
      endif
   end select
!-----------------------------

    enddo    ! j-loop

#endif

end subroutine camrad
!====================================================================
   subroutine camradinit(                                           &
                         r_d,r_v,cp,g,stbolt,ep_2,shalf,pptop,               &
#ifdef AMIPW_PHYSICS
                         hypm, &
#endif
                         ozmixm,pin,levsiz,xlat,num_months,         &
                         m_psp,m_psn,m_hybi,aerosolcp,aerosolcn,    &
                         paerlev,naer_c,                            &
                     ids, ide, jds, jde, kds, kde,                  &
                     ims, ime, jms, jme, kms, kme,                  &
                     its, ite, jts, jte, kts, kte                   )
#ifndef AMIPW_PHYSICS
   use module_wrf_error
   use module_state_description
   !use module_configure
#endif

!--------------------------------------------------------------------
   implicit none
!--------------------------------------------------------------------
   integer , intent(in)           :: ids, ide, jds, jde, kds, kde,  &
                                     ims, ime, jms, jme, kms, kme,  &
                                     its, ite, jts, jte, kts, kte
   real, intent(in)               :: pptop
#ifdef AMIPW_PHYSICS
   real, intent(in),  dimension( kms:kme )  :: hypm
#endif
   real, intent(in)               :: r_d,r_v,cp,g,stbolt,ep_2

   real,     dimension( kms:kme )  :: shalf

   integer,      intent(in   )    ::   levsiz, num_months
   integer,      intent(in   )    ::   paerlev, naer_c

   real, dimension( ims:ime, jms:jme ), intent(in   )  :: xlat

   real,  dimension( ims:ime, levsiz, jms:jme, num_months ),      &
          intent(inout   ) ::                                  ozmixm

   real,  dimension(levsiz), intent(inout )  ::                   pin
   real,  dimension(ims:ime, jms:jme), intent(inout )  ::                  m_psp,m_psn
   real,  dimension(paerlev), intent(inout )  ::               m_hybi
   real,  dimension( ims:ime, paerlev, jms:jme, naer_c ),      &
          intent(inout) ::                             aerosolcp,aerosolcn

   real(r8)    :: pstd
   real(r8)    :: rh2o, cpair

! these were made allocatable 20090612 to save static memory allocation. jm
   if ( .not. allocated( ksul   ) ) allocate( ksul( nrh, nspint ) )
   if ( .not. allocated( wsul   ) ) allocate( wsul( nrh, nspint ) )
   if ( .not. allocated( gsul   ) ) allocate( gsul( nrh, nspint ) )
   if ( .not. allocated( ksslt  ) ) allocate( ksslt( nrh, nspint ) )
   if ( .not. allocated( wsslt  ) ) allocate( wsslt( nrh, nspint ) )
   if ( .not. allocated( gsslt  ) ) allocate( gsslt( nrh, nspint ) )
   if ( .not. allocated( kcphil ) ) allocate( kcphil( nrh, nspint ) )
   if ( .not. allocated( wcphil ) ) allocate( wcphil( nrh, nspint ) )
   if ( .not. allocated( gcphil ) ) allocate( gcphil( nrh, nspint ) )

   if ( .not. allocated(ah2onw  ) ) allocate( ah2onw(n_p, n_tp, n_u, n_te, n_rh) )
   if ( .not. allocated(eh2onw  ) ) allocate( eh2onw(n_p, n_tp, n_u, n_te, n_rh) )
   if ( .not. allocated(ah2ow   ) ) allocate( ah2ow(n_p, n_tp, n_u, n_te, n_rh) )
   if ( .not. allocated(cn_ah2ow) ) allocate( cn_ah2ow(n_p, n_tp, n_u, n_te, n_rh) )
   if ( .not. allocated(cn_eh2ow) ) allocate( cn_eh2ow(n_p, n_tp, n_u, n_te, n_rh) )
   if ( .not. allocated(ln_ah2ow) ) allocate( ln_ah2ow(n_p, n_tp, n_u, n_te, n_rh) )
   if ( .not. allocated(ln_eh2ow) ) allocate( ln_eh2ow(n_p, n_tp, n_u, n_te, n_rh) )

#if !defined(mac_kludge)
   ozncyc = .true.
   indirect = .true.
   ixcldliq = 2
   ixcldice = 3
!#if (nmm_core != 1)
! aerosol array is not in the nmm registry 
!   since cam radiation not available to nmm (yet)
!   so this is blocked out to enable cam compilation with nmm
#ifdef AMIPW_PHYSICS
   idxsul = 1 !p_sul
   idxsslt = 2 !p_sslt
   idxdustfirst = 3 ! p_dust1
   idxocpho = 7 !p_ocpho
   idxcarbonfirst = 7 !p_ocpho
   idxbcpho = 8 !p_bcpho
   idxocphi = 9 !p_ocphi
   idxbcphi = 10 !p_bcphi
   idxbg = 11 !p_bg
   idxvolc = 12 !p_volc
#else
   idxsul = p_sul
   idxsslt = p_sslt
   idxdustfirst = p_dust1
   idxocpho = p_ocpho
   idxcarbonfirst = p_ocpho
   idxbcpho = p_bcpho
   idxocphi = p_ocphi
   idxbcphi = p_bcphi
   idxbg = p_bg
   idxvolc = p_volc
#endif

   pstd = 101325.0
! from physconst module
   mwdry = 28.966            ! molecular weight dry air ~ kg/kmole (shr_const_mwdair)
   mwco2 =  44.              ! molecular weight co2
   mwh2o = 18.016            ! molecular weight water vapor (shr_const_mwwv)
   mwch4 =  16.              ! molecular weight ch4
   mwn2o =  44.              ! molecular weight n2o
   mwf11 = 136.              ! molecular weight cfc11
   mwf12 = 120.              ! molecular weight cfc12
   cappa = r_d/cp
   rair = r_d
   tmelt = 273.16            ! freezing t of fresh water ~ k 
   r_universal = 6.02214e26 * stbolt   ! universal gas constant ~ j/k/kmole
   latvap = 2.501e6          ! latent heat of evaporation ~ j/kg
   latice = 3.336e5          ! latent heat of fusion ~ j/kg
   zvir = r_v/r_d - 1.
   rh2o = r_v
   cpair = cp
!
   epsqs = ep_2

   call radini(g, cp, ep_2, stbolt, pstd*10.0 )
   call esinti(epsqs  ,latvap  ,latice  ,rh2o    ,cpair   ,tmelt   )
   call oznini(ozmixm,pin,levsiz,num_months,xlat,                   &
                     ids, ide, jds, jde, kds, kde,                  &
                     ims, ime, jms, jme, kms, kme,                  &
                     its, ite, jts, jte, kts, kte)                   
   call aerosol_init(m_psp,m_psn,m_hybi,aerosolcp,aerosolcn,paerlev,naer_c,shalf,pptop,    &
#ifdef AMIPW_PHYSICS
                     hypm, &
#endif
                     ids, ide, jds, jde, kds, kde,                  &
                     ims, ime, jms, jme, kms, kme,                  &
                     its, ite, jts, jte, kts, kte)

#endif

   end subroutine camradinit
#if !defined(mac_kludge)


subroutine oznint(julday,julian,dt,gmt,xtime,ozmixmj,ozmix,levsiz,num_months,pcols)

      implicit none

   integer,      intent(in   )    ::   levsiz, num_months,pcols

   real(r8),  dimension( pcols, levsiz, num_months ),      &
          intent(in   ) ::                                  ozmixmj 

   real, intent(in    )      ::        xtime,gmt
   integer, intent(in )      ::        julday
   real,    intent(in )      ::        julian
   real,    intent(in )      ::        dt

   real(r8),  dimension( pcols, levsiz ),      &
          intent(out  ) ::                                  ozmix
   !local
   real(r8)  :: intjulian
   integer   :: np1,np,nm,m,k,i
   integer   :: ijul
   integer, dimension(12) ::  date_oz
   data date_oz/16, 45, 75, 105, 136, 166, 197, 228, 258, 289, 319, 350/
   real(r8) :: cdayozp, cdayozm
   real(r8) :: fact1, fact2
   logical  :: finddate
   character(len=256) :: msgstr

   ! julian starts from 0.0 at 0z on 1 jan.
   intjulian = julian + 1.0_r8    ! offset by one day
! jan 1st 00z is julian=1.0 here
   ijul=int(intjulian)
!  note that following will drift. 
!    need to use actual month/day info to compute julian.
   intjulian=intjulian-float(ijul)
   ijul=mod(ijul,365)
   if(ijul.eq.0)ijul=365
   intjulian=intjulian+ijul
   np1=1
   finddate=.false.
!   do m=1,num_months
   do m=1,12
   if(date_oz(m).gt.intjulian.and..not.finddate) then
     np1=m
     finddate=.true.
   endif
   enddo
   cdayozp=date_oz(np1)
   if(np1.gt.1) then
   cdayozm=date_oz(np1-1)
   np=np1
   nm=np-1
   else
   cdayozm=date_oz(12)
   np=np1
   nm=12
   endif
   call getfactors(ozncyc,np1, cdayozm, cdayozp,intjulian, &
                    fact1, fact2) 

!
! time interpolation.
!
      do k=1,levsiz
         do i=1,pcols
            ozmix(i,k) = ozmixmj(i,k,nm)*fact1 + ozmixmj(i,k,np)*fact2
         end do
      end do

end subroutine oznint


subroutine get_aerosol(c, julday, julian, dt, gmt, xtime, m_psp, m_psn, aerosoljp, &
  aerosoljn, m_hybi, paerlev, naer_c, pint, pcols, pver, pverp, pverr, pverrp, aerosolt, scale)
!------------------------------------------------------------------
!
!  input:
!     time at which aerosol mmrs are needed (get_curr_calday())
!     chunk index
!     cam's vertical grid (pint)
!
!  output:
!     values for aerosol mass mixing ratios at specified time
!     on vertical grid specified by cam (aerosolt)
!
!  method:
!     first determine which indexs of aerosols are the bounding data sets
!     interpolate both onto vertical grid aerm(),aerp().
!     from those two, interpolate in time.
!
!------------------------------------------------------------------

!  use volcanicmass, only: get_volcanic_mass
!  use timeinterp, only: getfactors
!
! aerosol fields interpolated to current time step
!   on pressure levels of this time step.
! these should be made read-only for other modules
! is allocation done correctly here?
!
   integer, intent(in) :: c                   ! chunk id.
   integer, intent(in) :: paerlev, naer_c, pcols, pver, pverp, pverr, pverrp
   real(r8), intent(in) :: pint(pcols,pverp)  ! midpoint pres.
   real(r8), intent(in) :: scale(naer_all)    ! scale each aerosol by this amount
   real, intent(in    )      ::        xtime,gmt
   integer, intent(in )      ::        julday
   real, intent(in    )      ::        julian
   real, intent(in    )      ::        dt
   real(r8), intent(in   )      ::        m_psp(pcols),m_psn(pcols)  ! match surface pressure
   real(r8), intent(in   )   ::        aerosoljp(pcols,paerlev,naer_c) 
   real(r8), intent(in   )   ::        aerosoljn(pcols,paerlev,naer_c) 
   real(r8), intent(in   )   ::        m_hybi(paerlev)

   real(r8), intent(out) :: aerosolt(pcols, pver, naer_all) ! aerosols
!
! local workspace
!
   real(r8) caldayloc                     ! calendar day of current timestep
   real(r8) fact1, fact2                  ! time interpolation factors

  integer :: nm = 1                ! index to prv month in array. init to 1 and toggle between 1 and 2
  integer :: np = 2                ! index to nxt month in array. init to 2 and toggle between 1 and 2
  integer :: mo_nxt = bigint       ! index to nxt month in file
  integer :: mo_prv                       ! index to previous month

  real(r8) :: cdaym = inf          ! calendar day of prv month
  real(r8) :: cdayp = inf          ! calendar day of next month
  real(r8) :: mid(12)              ! days into year for mid month date
  data mid/16.5, 46.0, 75.5, 106.0, 136.5, 167.0, 197.5, 228.5, 259.0, 289.5, 320.0, 350.5 /

   integer i, k, j                        ! spatial indices
   integer m                              ! constituent index
   integer lats(pcols),lons(pcols)        ! latitude and longitudes of column
   integer ncol                           ! number of columns
   integer ijul
   real(r8) intjulian

   real(r8) speciesmin(naer)              ! minimal value for each species
!
! values before current time step "the minus month"
! aerosolm(pcols,pver) is value of preceeding month's aerosol mmr
! aerosolp(pcols,pver) is value of next month's aerosol mmr
!  (think minus and plus or values to left and right of point to be interpolated)
!
   real(r8) aerosolm(pcols,pver,naer) ! aerosol mmr from match in column at previous (minus) month
!
! values beyond (or at) current time step "the plus month"
!
   real(r8) aerosolp(pcols,pver,naer) ! aerosol mmr from match in column at next (plus) month
   character(len=256) :: msgstr

   ! julian starts from 0.0 at 0z on 1 jan.
   intjulian = julian + 1.0_r8    ! offset by one day
! jan 1st 00z is julian=1.0 here
   ijul=int(intjulian)
!  note that following will drift. 
!    need to use actual month/day info to compute julian.
   intjulian=intjulian-float(ijul)
   ijul=mod(ijul,365)
   if(ijul.eq.0)ijul=365
   caldayloc=intjulian+ijul

   if (caldayloc < mid(1)) then
      mo_prv = 12
      mo_nxt =  1
   else if (caldayloc >= mid(12)) then
      mo_prv = 12
      mo_nxt =  1
   else
      do i = 2 , 12
         if (caldayloc < mid(i)) then
            mo_prv = i-1
            mo_nxt = i
            exit
         end if
      end do
   end if
!
! set initial calendar day values
!
   cdaym = mid(mo_prv)
   cdayp = mid(mo_nxt)

!
! determine time interpolation factors.  1st arg says we are cycling 1 year of data
!
   call getfactors (.true., mo_nxt, cdaym, cdayp, caldayloc, &
                    fact1, fact2)
!
! interpolate (prv and nxt month) bounding datasets onto cam vertical grid.
! compute mass mixing ratios on cams's pressure coordinate
!  for both the "minus" and "plus" months
!
!  ncol = get_ncols_p(c)
   ncol = pcols

!  call vert_interpolate (m_ps_cam_col(1,c,nm), pint, nm, aerosolm, ncol, c)
!  call vert_interpolate (m_ps_cam_col(1,c,np), pint, np, aerosolp, ncol, c)

   call vert_interpolate (m_psp, aerosoljp, m_hybi, paerlev, naer_c, pint, nm, aerosolm, pcols, pver, pverp, ncol, c)
   call vert_interpolate (m_psn, aerosoljn, m_hybi, paerlev, naer_c, pint, np, aerosolp, pcols, pver, pverp, ncol, c)

!
! time interpolate.
!
   do m=1,naer
      do k=1,pver
         do i=1,ncol
            aerosolt(i,k,m) = aerosolm(i,k,m)*fact1 + aerosolp(i,k,m)*fact2
         end do
      end do
   end do

!  do i=1,ncol
!     match_ps_chunk(i,c) = m_ps(i,nm)*fact1 + m_ps(i,np)*fact2
!  end do
!
! get background aerosol (tuning) field
!
   call background (c, ncol, pint, pcols, pverr, pverrp, aerosolt(:, :, idxbg))

!
! find volcanic aerosol masses
!
! if (strat_volcanic) then
!   call get_volcanic_mass(c, aerosolt(:,:,idxvolc))
! else
    aerosolt(:,:,idxvolc) = 0._r8
! endif

!
! exit if mmr is negative (we have previously set
!  cumulative mass to be a decreasing function.)
!
   speciesmin(:) = 0. ! speciesmin(m) = 0 is minimum mmr for each species

   do m=1,naer
      do k=1,pver
         do i=1,ncol
            if (aerosolt(i, k, m) < speciesmin(m)) then
               write(6,*) 'aerosol_interpolate: negative mass mixing ratio, exiting'
               write(6,*) 'm, column, pver',m, i, k ,aerosolt(i, k, m)
               call endrun ('')
            end if
         end do
      end do
   end do
!
! scale any aerosols as required
!
   call scale_aerosols (aerosolt, pcols, pver, ncol, c, scale)

   return
end subroutine get_aerosol


subroutine aerosol_indirect(ncol,lchnk,pcols,pver,ppcnst,landfrac,pmid,t,qm1,cld,zm,rel)
!--------------------------------------------------------------
! compute effect of sulfate on effective liquid water radius
!  method of martin et. al.
!--------------------------------------------------------------

! use constituents, only: ppcnst, cnst_get_ind
! use history, only: outfld

!#include <comctl.h>

  integer, intent(in) :: ncol                  ! number of atmospheric columns
  integer, intent(in) :: lchnk                 ! chunk identifier
  integer, intent(in) :: pcols,pver,ppcnst

  real(r8), intent(in) :: landfrac(pcols)      ! land fraction
  real(r8), intent(in) :: pmid(pcols,pver)     ! model level pressures
  real(r8), intent(in) :: t(pcols,pver)        ! model level temperatures
  real(r8), intent(in) :: qm1(pcols,pver,ppcnst) ! specific humidity and tracers
  real(r8), intent(in) :: cld(pcols,pver)      ! fractional cloud cover
  real(r8), intent(in) :: zm(pcols,pver)       ! height of midpoints (above surface)
  real(r8), intent(in) :: rel(pcols,pver)      ! liquid effective drop size (microns)
!
! local variables
!
  real(r8) locrhoair(pcols,pver)  ! dry air density            [kg/m^3 ]
  real(r8) lwcwat(pcols,pver)     ! in-cloud liquid water path [kg/m^3 ]
  real(r8) sulfmix(pcols,pver)    ! sulfate mass mixing ratio  [kg/kg  ]
  real(r8) so4mass(pcols,pver)    ! sulfate mass concentration [g/cm^3 ]
  real(r8) aso4(pcols,pver)       ! sulfate # concentration    [#/cm^3 ]
  real(r8) ntot(pcols,pver)       ! ccn # concentration        [#/cm^3 ]
  real(r8) relmod(pcols,pver)     ! effective radius           [microns]

  real(r8) wrel(pcols,pver)       ! weighted effective radius    [microns]
  real(r8) wlwc(pcols,pver)       ! weighted liq. water content  [kg/m^3 ]
  real(r8) cldfrq(pcols,pver)     ! frequency of occurance of...
!                                  ! clouds (cld => 0.01)         [fraction]
  real(r8) locpi                  ! my piece of the pi
  real(r8) rdryair                ! gas constant of dry air   [j/deg/kg]
  real(r8) rhowat                 ! density of water          [kg/m^3  ]
  real(r8) acoef                  ! m->a conversion factor; assumes
!                                  ! dbar=0.10, sigma=2.0      [g^-1    ]
  real(r8) rekappa                ! kappa in evaluation of re(lmod)
  real(r8) recoef                 ! temp. coeficient for calc of re(lmod)
  real(r8) reexp                  ! 1.0/3.0
  real(r8) ntotb                  ! temp var to hold below cloud ccn
! -- parameters for background cdnc (from `ambient' non-sulfate aerosols)...
  real(r8) cmarn                  ! coef for cdnc_marine         [cm^-3]
  real(r8) cland                  ! coef for cdnc_land           [cm^-3]
  real(r8) hmarn                  ! scale height for cdnc_marine [m]
  real(r8) hland                  ! scale height for cdnc_land   [m]
  parameter ( cmarn = 50.0, cland = 100.0 )
  parameter ( hmarn = 1000.0, hland = 2000.0 )
  real(r8) bgaer                  ! temp var to hold background cdnc

  integer i,k     ! loop indices
!
! statement functions
!
  logical land    ! is this a column over land?
  land(i) = nint(landfrac(i)).gt.0.5_r8

  if (indirect) then

!   call endrun ('aerosol_indirect:  indirect effect is obsolete')

!   ramping is not yet resolved so sulfmix is 0.
    sulfmix(1:ncol,1:pver) = 0._r8

    locpi = 3.141592654
    rdryair = 287.04
    rhowat = 1000.0
    acoef = 1.2930e14
    recoef = 3.0/(4.0*locpi*rhowat)
    reexp = 1.0/3.0

!   call cnst_get_ind('cldliq', ixcldliq)
    do k=pver,1,-1
      do i = 1,ncol
        locrhoair(i,k) = pmid(i,k)/( rdryair*t(i,k) )
        lwcwat(i,k) = ( qm1(i,k,ixcldliq)/max(0.01_r8,cld(i,k)) )* &
                      locrhoair(i,k)
!          note: 0.001 converts kg/m3 -> g/cm3
        so4mass(i,k) = sulfmix(i,k)*locrhoair(i,k)*0.001
        aso4(i,k) = so4mass(i,k)*acoef

        if (aso4(i,k) <= 280.0) then
           aso4(i,k) = max(36.0_r8,aso4(i,k))
           ntot(i,k) = -1.15e-3*aso4(i,k)**2 + 0.963*aso4(i,k)+5.30
           rekappa = 0.80
        else
           aso4(i,k) = min(1500.0_r8,aso4(i,k))
           ntot(i,k) = -2.10e-4*aso4(i,k)**2 + 0.568*aso4(i,k)-27.9
           rekappa = 0.67
        end if
        if (land(i)) then ! account for local background aerosol;
           bgaer = cland*exp(-(zm(i,k)/hland))
           ntot(i,k) = max(bgaer,ntot(i,k))
        else
           bgaer = cmarn*exp(-(zm(i,k)/hmarn))
           ntot(i,k) = max(bgaer,ntot(i,k))
        end if

        if (k == pver) then
           ntotb = ntot(i,k)
        else
           ntotb = ntot(i,k+1)
        end if

        relmod(i,k) = (( (recoef*lwcwat(i,k))/(rekappa*ntotb))**reexp)*10000.0
        relmod(i,k) = max(4.0_r8,relmod(i,k))
        relmod(i,k) = min(20.0_r8,relmod(i,k))
        if (cld(i,k) >= 0.01) then
           cldfrq(i,k) = 1.0
        else
           cldfrq(i,k) = 0.0
        end if
        wrel(i,k) = relmod(i,k)*cldfrq(i,k)
        wlwc(i,k) = lwcwat(i,k)*cldfrq(i,k)
      end do
    end do
!   call outfld('mso4    ',so4mass,pcols,lchnk)
!   call outfld('lwc     ',lwcwat ,pcols,lchnk)
!   call outfld('cldfrq  ',cldfrq ,pcols,lchnk)
!   call outfld('wrel    ',wrel   ,pcols,lchnk)
!   call outfld('wlwc    ',wlwc   ,pcols,lchnk)
!   write(6,*)'warning: indirect calculation has no effects'
  else
    do k = 1, pver
      do i = 1, ncol
        relmod(i,k) = rel(i,k)
      end do
    end do
  endif

! call outfld('rel     ',relmod ,pcols,lchnk)

  return
end subroutine aerosol_indirect


      subroutine aer_trn(aer_mpp, aer_trn_ttl, pcols, plev, plevp )
!
!     purpose: compute strat. aerosol transmissions needed in absorptivity/
!              emissivity calculations
!              aer_trn() is called by radclw() when doabsems is .true.
!
!     use shr_kind_mod, only: r8 => shr_kind_r8
!     use pmgrid
!     use ppgrid
!     use prescribed_aerosols, only: strat_volcanic
      implicit none

!     input arguments
!
!       [kg m-2] volcanics path above kth interface level
!
      integer, intent(in)         :: pcols, plev, plevp
      real(r8), intent(in) :: aer_mpp(pcols,plevp)

!     output arguments
!
!       [fraction] total volcanic transmission between interfaces k1 and k2
!
      real(r8), intent(out) ::  aer_trn_ttl(pcols,plevp,plevp,bnd_nbr_lw)

!-------------------------------------------------------------------------
!     local variables

      integer bnd_idx           ! lw band index
      integer i                 ! lon index
      integer k1                ! lev index
      integer k2                ! lev index
      real(r8) aer_pth_dlt      ! [kg m-2] volcanics path between interface
                                !          levels k1 and k2
      real(r8) odap_aer_ttl     ! [fraction] total path absorption optical
                                !            depth

!-------------------------------------------------------------------------

      if (strat_volcanic) then
        do bnd_idx=1,bnd_nbr_lw
           do i=1,pcols
              aer_trn_ttl(i,1,1,bnd_idx)=1.0
           end do
           do k1=2,plevp
              do i=1,pcols
                 aer_trn_ttl(i,k1,k1,bnd_idx)=1.0

                 aer_pth_dlt  = abs(aer_mpp(i,k1) - aer_mpp(i,1))
                 odap_aer_ttl = abs_cff_mss_aer(bnd_idx) * aer_pth_dlt

                 aer_trn_ttl(i,1,k1,bnd_idx) = exp(-1.66 * odap_aer_ttl)
              end do
           end do

           do k1=2,plev
              do k2=k1+1,plevp
                 do i=1,pcols
                    aer_trn_ttl(i,k1,k2,bnd_idx) = &
                         aer_trn_ttl(i,1,k2,bnd_idx) / &
                         aer_trn_ttl(i,1,k1,bnd_idx)
                 end do
              end do
           end do

           do k1=2,plevp
              do k2=1,k1-1
                 do i=1,pcols
                    aer_trn_ttl(i,k1,k2,bnd_idx)=aer_trn_ttl(i,k2,k1,bnd_idx)
                 end do
              end do
           end do
        end do
      else
        aer_trn_ttl = 1.0
      endif

      return
      end subroutine aer_trn

      subroutine aer_pth(aer_mass, aer_mpp, ncol, pcols, plev, plevp)
!------------------------------------------------------
!     purpose: convert mass per layer to cumulative mass from top
!------------------------------------------------------
!     use shr_kind_mod, only: r8 => shr_kind_r8
!     use ppgrid
!     use pmgrid
      implicit none
!#include <crdcon.h>

!     parameters
!     input
      integer, intent(in)        :: pcols, plev, plevp
      real(r8), intent(in):: aer_mass(pcols,plev)  ! rad level aerosol mass mixing ratio
      integer, intent(in):: ncol
!
!     output
      real(r8), intent(out):: aer_mpp(pcols,plevp) ! [kg m-2] volcanics path above kth interface
!
!     local
      integer i      ! column index
      integer k      ! level index
!------------------------------------------------------
!------------------------------------------------------

      aer_mpp(1:ncol,1) =  0._r8
      do k=2,plevp
          aer_mpp(1:ncol,k) = aer_mpp(1:ncol,k-1) + aer_mass(1:ncol,k-1)
      enddo
!
      return
      end subroutine aer_pth

subroutine radctl(j, lchnk   ,ncol    , pcols, pver, pverp, pverr, pverrp, ppcnst, pcnst,  &
                  lwups   ,emis    ,          &
                  pmid    ,pint    ,pmln    ,piln    ,pdel    ,t       , &
!                 qm1     ,cld     ,cicewp  ,cliqwp  ,coszrs,  clat, &
                  qm1     ,cld     ,cicewp  ,cliqwp  ,tauxcl, tauxci, coszrs,  clat, &
                  asdir   ,asdif   ,aldir   ,aldif   ,solcon, gmt,julday,julian,dt,xtime,  &
                  pin, ozmixmj, ozmix, levsiz, num_months,      &
                  m_psp, m_psn,  aerosoljp, aerosoljn, m_hybi, paerlev, naer_c, pmxrgn  , &
                  nmxrgn  ,                   &
                  dolw, dosw, doabsems, abstot, absnxt, emstot, &
                  fsup    ,fsupc   ,fsdn    ,fsdnc   , &
                  flup    ,flupc   ,fldn    ,fldnc   , &
                  swcf    ,lwcf    ,flut    ,          &
                  fsns    ,fsnt    ,flns    ,flnt    , &
                  qrs     ,qrl     ,flwds   ,rel     ,rei     , &
                  sols    ,soll    ,solsd   ,solld   , &
                  landfrac,zm      ,fsds     )
!----------------------------------------------------------------------- 
! 
! purpose: 
! driver for radiation computation.
! 
! method: 
! radiation uses cgs units, so conversions must be done from
! model fields to radiation fields.
!
! author: ccm1,  cms contact: j. truesdale
! 
!-----------------------------------------------------------------------
!  use shr_kind_mod, only: r8 => shr_kind_r8
!  use ppgrid
!  use pspect
!  use commap
!  use history, only: outfld
!  use constituents, only: ppcnst, cnst_get_ind
!  use prescribed_aerosols, only: get_aerosol, naer_all, aerosol_diagnostics, &
!     aerosol_indirect, get_rf_scales, get_int_scales, radforce, idxvolc
!  use physics_types, only: physics_state
!  use wv_saturation, only: aqsat
!  use chemistry,    only: trace_gas
!  use physconst, only: cpair, epsilo
!  use aer_optics, only: idxvis
!  use aerosol_intr, only: set_aerosol_from_prognostics


   implicit none

!
! input arguments
!
   integer, intent(in) :: lchnk,j                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns
   integer, intent(in) :: levsiz                ! number of ozone data levels
   integer, intent(in) :: num_months            ! 12 months
   integer, intent(in) :: paerlev,naer_c          ! aerosol vertical level and # species
   integer, intent(in) :: pcols, pver, pverp, pverr, pverrp, ppcnst, pcnst
   logical, intent(in) :: dolw,dosw,doabsems


   integer nspint            ! num of spctrl intervals across solar spectrum
   integer naer_groups       ! num of aerosol groups for optical diagnostics
   parameter ( nspint = 19 )
   parameter ( naer_groups = 7 )    ! current groupings are sul, sslt, all carbons, all dust, background, and all aerosols


   real(r8), intent(in) :: lwups(pcols)         ! longwave up flux at surface
   real(r8), intent(in) :: emis(pcols,pver)     ! cloud emissivity
   real(r8), intent(in) :: pmid(pcols,pver)     ! model level pressures
   real(r8), intent(in) :: pint(pcols,pverp)    ! model interface pressures
   real(r8), intent(in) :: pmln(pcols,pver)     ! natural log of pmid
   real(r8), intent(in) :: rel(pcols,pver)      ! liquid effective drop size (microns)
   real(r8), intent(in) :: rei(pcols,pver)      ! ice effective drop size (microns)
   real(r8), intent(in) :: piln(pcols,pverp)    ! natural log of pint
   real(r8), intent(in) :: pdel(pcols,pverp)    ! pressure difference across layer 
   real(r8), intent(in) :: t(pcols,pver)        ! model level temperatures
   real(r8), intent(in) :: qm1(pcols,pver,ppcnst) ! specific humidity and tracers
   real(r8), intent(in) :: cld(pcols,pver)      ! fractional cloud cover
   real(r8), intent(in) :: cicewp(pcols,pver)   ! in-cloud cloud ice water path
   real(r8), intent(in) :: cliqwp(pcols,pver)   ! in-cloud cloud liquid water path
   real(r8), intent(inout) :: tauxcl(pcols,0:pver) ! cloud water optical depth
   real(r8), intent(inout) :: tauxci(pcols,0:pver) ! cloud ice optical depth
   real(r8), intent(in) :: coszrs(pcols)        ! cosine solar zenith angle
   real(r8), intent(in) :: clat(pcols)          ! latitude in radians for columns 
   real(r8), intent(in) :: asdir(pcols)         ! albedo shortwave direct
   real(r8), intent(in) :: asdif(pcols)         ! albedo shortwave diffuse
   real(r8), intent(in) :: aldir(pcols)         ! albedo longwave direct
   real(r8), intent(in) :: aldif(pcols)         ! albedo longwave diffuse
   real(r8), intent(in) :: landfrac(pcols)      ! land fraction
   real(r8), intent(in) :: zm(pcols,pver)       ! height of midpoints (above surface)
   real(r8), intent(in) :: pin(levsiz)          ! pressure levels of ozone data
   real(r8), intent(in) :: ozmixmj(pcols,levsiz,num_months)  ! monthly ozone mixing ratio
   real(r8), intent(inout) :: ozmix(pcols,levsiz)  ! ozone data
   real, intent(in) :: solcon               ! solar constant with eccentricity factor
   real, intent(in    )      ::        xtime,gmt              
   integer, intent(in )      ::        julday
   real,    intent(in )      ::        julian
   real,    intent(in )      ::        dt
   real(r8), intent(in)     :: m_psp(pcols),m_psn(pcols)       ! match surface pressure
   real(r8), intent(in)     :: aerosoljp(pcols,paerlev,naer_c)   ! aerosol concentrations
   real(r8), intent(in)     :: aerosoljn(pcols,paerlev,naer_c)   ! aerosol concentrations
   real(r8), intent(in)     :: m_hybi(paerlev)
!  type(physics_state), intent(in) :: state     
   real(r8), intent(inout) :: pmxrgn(pcols,pverp) ! maximum values of pmid for each
!    maximally overlapped region.
!    0->pmxrgn(i,1) is range of pmid for
!    1st region, pmxrgn(i,1)->pmxrgn(i,2) for
!    2nd region, etc
   integer, intent(inout) :: nmxrgn(pcols)     ! number of maximally overlapped regions

    real(r8) :: pmxrgnrf(pcols,pverp)             ! temporary copy of pmxrgn
    integer  :: nmxrgnrf(pcols)     ! temporary copy of nmxrgn

!
! output solar arguments
!
   real(r8), intent(out) :: fsns(pcols)          ! surface absorbed solar flux
   real(r8), intent(out) :: fsnt(pcols)          ! net column abs solar flux at model top
   real(r8), intent(out) :: flns(pcols)          ! srf longwave cooling (up-down) flux
   real(r8), intent(out) :: flnt(pcols)          ! net outgoing lw flux at model top
   real(r8), intent(out) :: sols(pcols)          ! downward solar rad onto surface (sw direct)
   real(r8), intent(out) :: soll(pcols)          ! downward solar rad onto surface (lw direct)
   real(r8), intent(out) :: solsd(pcols)         ! downward solar rad onto surface (sw diffuse)
   real(r8), intent(out) :: solld(pcols)         ! downward solar rad onto surface (lw diffuse)
   real(r8), intent(out) :: qrs(pcols,pver)      ! solar heating rate
   real(r8), intent(out) :: fsds(pcols)          ! flux shortwave downwelling surface
! added outputs of total and clearsky fluxes etc
   real(r8), intent(out) :: fsup(pcols,pverp)    ! upward total sky solar
   real(r8), intent(out) :: fsupc(pcols,pverp)   ! upward clear sky solar
   real(r8), intent(out) :: fsdn(pcols,pverp)    ! downward total sky solar
   real(r8), intent(out) :: fsdnc(pcols,pverp)   ! downward clear sky solar
   real(r8), intent(out) :: flup(pcols,pverp)    ! upward total sky longwave
   real(r8), intent(out) :: flupc(pcols,pverp)   ! upward clear sky longwave
   real(r8), intent(out) :: fldn(pcols,pverp)    ! downward total sky longwave
   real(r8), intent(out) :: fldnc(pcols,pverp)   ! downward clear sky longwave
   real(r8), intent(out) :: swcf(pcols)          ! top of the atmosphere solar cloud forcing
   real(r8), intent(out) :: lwcf(pcols)          ! top of the atmosphere longwave cloud forcing
   real(r8), intent(out) :: flut(pcols)          ! top of the atmosphere outgoing longwave
!
! output longwave arguments
!
   real(r8), intent(out) :: qrl(pcols,pver)      ! longwave cooling rate
   real(r8), intent(out) :: flwds(pcols)         ! surface down longwave flux

   real(r8), intent(inout) :: abstot(pcols,pverp,pverp) ! total absorptivity
   real(r8), intent(inout) :: absnxt(pcols,pver,4)      ! total nearest layer absorptivity
   real(r8), intent(inout) :: emstot(pcols,pverp)     ! total emissivity


!
!---------------------------local variables-----------------------------
!
   integer i, k              ! index

   integer :: in2o, ich4, if11, if12 ! indexes of gases in constituent array

   real(r8) solin(pcols)         ! solar incident flux
!  real(r8) fsds(pcols)          ! flux shortwave downwelling surface
   real(r8) fsntoa(pcols)        ! net solar flux at toa
   real(r8) fsntoac(pcols)       ! clear sky net solar flux at toa
   real(r8) fsnirt(pcols)        ! near-ir flux absorbed at toa
   real(r8) fsnrtc(pcols)        ! clear sky near-ir flux absorbed at toa
   real(r8) fsnirtsq(pcols)      ! near-ir flux absorbed at toa >= 0.7 microns
   real(r8) fsntc(pcols)         ! clear sky total column abs solar flux
   real(r8) fsnsc(pcols)         ! clear sky surface abs solar flux
   real(r8) fsdsc(pcols)         ! clear sky surface downwelling solar flux
!  real(r8) flut(pcols)          ! upward flux at top of model
!  real(r8) lwcf(pcols)          ! longwave cloud forcing
!  real(r8) swcf(pcols)          ! shortwave cloud forcing
   real(r8) flutc(pcols)         ! upward clear sky flux at top of model
   real(r8) flntc(pcols)         ! clear sky lw flux at model top
   real(r8) flnsc(pcols)         ! clear sky lw flux at srf (up-down)
   real(r8) ftem(pcols,pver)     ! temporary array for outfld

   real(r8) pbr(pcols,pverr)     ! model mid-level pressures (dynes/cm2)
   real(r8) pnm(pcols,pverrp)    ! model interface pressures (dynes/cm2)
   real(r8) o3vmr(pcols,pverr)   ! ozone volume mixing ratio
   real(r8) o3mmr(pcols,pverr)   ! ozone mass mixing ratio
   real(r8) eccf                 ! earth/sun distance factor
   real(r8) n2o(pcols,pver)      ! nitrous oxide mass mixing ratio
   real(r8) ch4(pcols,pver)      ! methane mass mixing ratio
   real(r8) cfc11(pcols,pver)    ! cfc11 mass mixing ratio
   real(r8) cfc12(pcols,pver)    ! cfc12 mass mixing ratio
   real(r8) rh(pcols,pverr)      ! level relative humidity (fraction)
   real(r8) lwupcgs(pcols)       ! upward longwave flux in cgs units

   real(r8) esat(pcols,pverr)    ! saturation vapor pressure
   real(r8) qsat(pcols,pverr)    ! saturation specific humidity

   real(r8) :: frc_day(pcols) ! = 1 for daylight, =0 for night colums
   real(r8) :: aertau(pcols,nspint,naer_groups) ! aerosol column optical depth
   real(r8) :: aerssa(pcols,nspint,naer_groups) ! aerosol column averaged single scattering albedo
   real(r8) :: aerasm(pcols,nspint,naer_groups) ! aerosol column averaged asymmetry parameter
   real(r8) :: aerfwd(pcols,nspint,naer_groups) ! aerosol column averaged forward scattering

   real(r8) aerosol(pcols, pver, naer_all) ! aerosol mass mixing ratios
   real(r8) scales(naer_all)               ! scaling factors for aerosols


!
! interpolate ozone volume mixing ratio to model levels
!
! wrf: added pin, levsiz, ozmix here
   call oznint(julday,julian,dt,gmt,xtime,ozmixmj,ozmix,levsiz,num_months,pcols)

   call radozn(lchnk   ,ncol    &
              ,pcols, pver &
              ,pmid    ,pin, levsiz, ozmix, o3vmr   )

!  call outfld('o3vmr   ',o3vmr ,pcols, lchnk)

!
! set chunk dependent radiation input
!
   call radinp(lchnk   ,ncol    ,pcols, pver, pverp,      &
               pmid    ,pint    ,o3vmr   , pbr     ,&
               pnm     ,eccf    ,o3mmr   )

!
! solar radiation computation
!
   if (dosw) then

!
! calculate heating with aerosols
!
      call aqsat(t, pmid, esat, qsat, pcols, &
                 ncol, pver, 1, pver)

      ! calculate relative humidity
!     rh(1:ncol,1:pver) = q(1:ncol,1:pver,1) / qsat(1:ncol,1:pver) * &
!        ((1.0 - epsilo) * qsat(1:ncol,1:pver) + epsilo) / &
!        ((1.0 - epsilo) * q(1:ncol,1:pver,1) + epsilo)
      rh(1:ncol,1:pver) = qm1(1:ncol,1:pver,1) / qsat(1:ncol,1:pver) * &
         ((1.0 - epsilo) * qsat(1:ncol,1:pver) + epsilo) / &
         ((1.0 - epsilo) * qm1(1:ncol,1:pver,1) + epsilo)

      if (radforce) then

         pmxrgnrf = pmxrgn
         nmxrgnrf = nmxrgn

         call get_rf_scales(scales)

         call get_aerosol(lchnk, julday, julian, dt, gmt, xtime, m_psp, m_psn, aerosoljp, &
           aerosoljn, m_hybi, paerlev, naer, pint, pcols, pver, pverp, pverr, pverrp, aerosol, scales)

         ! overwrite with prognostics aerosols

!   no feedback from prognostic aerosols 
!        call set_aerosol_from_prognostics (ncol, q, aerosol)

         call aerosol_indirect(ncol,lchnk,pcols,pver,ppcnst,landfrac,pmid,t,qm1,cld,zm,rel)
   
!        call t_startf('radcswmx_rf')
         call radcswmx(j, lchnk   ,ncol ,pcols, pver, pverp,         &
                    pnm     ,pbr     ,qm1     ,rh      ,o3mmr   , &
                    aerosol ,cld     ,cicewp  ,cliqwp  ,rel     , &
!                   rei     ,eccf    ,coszrs  ,scon    ,solin   ,solcon , &
                    rei     ,tauxcl  ,tauxci  ,eccf    ,coszrs  ,scon    ,solin   ,solcon , &
                    asdir   ,asdif   ,aldir   ,aldif   ,nmxrgnrf, &
                    pmxrgnrf,qrs     ,fsnt    ,fsntc   ,fsntoa  , &
                    fsntoac ,fsnirt  ,fsnrtc  ,fsnirtsq,fsns    , &
                    fsnsc   ,fsdsc   ,fsds    ,sols    ,soll    , &
                    solsd   ,solld   ,frc_day ,                   &
                    fsup    ,fsupc   ,fsdn    ,fsdnc   ,          &
                    aertau  ,aerssa  ,aerasm  ,aerfwd             )
!        call t_stopf('radcswmx_rf')

!
! convert units of shortwave fields needed by rest of model from cgs to mks
!

            do i = 1, ncol
            solin(i) = solin(i)*1.e-3
            fsnt(i)  = fsnt(i) *1.e-3
            fsns(i)  = fsns(i) *1.e-3
            fsntc(i) = fsntc(i)*1.e-3
            fsnsc(i) = fsnsc(i)*1.e-3
            end do
         ftem(:ncol,:pver) = qrs(:ncol,:pver)/cpair
!
! dump shortwave radiation information to history tape buffer (diagnostics)
!
!        call outfld('qrs_rf  ',ftem  ,pcols,lchnk)
!        call outfld('fsnt_rf ',fsnt  ,pcols,lchnk)
!        call outfld('fsns_rf ',fsns  ,pcols,lchnk)
!        call outfld('fsntc_rf',fsntc ,pcols,lchnk)
!        call outfld('fsnsc_rf',fsnsc ,pcols,lchnk)
 
      endif ! if (radforce)

      call get_int_scales(scales)

      call get_aerosol(lchnk, julday, julian, dt, gmt, xtime, m_psp, m_psn, aerosoljp, aerosoljn, &
             m_hybi, paerlev, naer, pint, pcols, pver, pverp, pverr, pverrp, aerosol, scales)

      ! overwrite with prognostics aerosols
!     call set_aerosol_from_prognostics (ncol, q, aerosol)

      call aerosol_indirect(ncol,lchnk,pcols,pver,ppcnst,landfrac,pmid,t,qm1,cld,zm,rel)
!     call t_startf('radcswmx')

      call radcswmx(j, lchnk   ,ncol    ,pcols, pver, pverp,         &
                    pnm     ,pbr     ,qm1     ,rh      ,o3mmr   , &
                    aerosol ,cld     ,cicewp  ,cliqwp  ,rel     , &
!                   rei     ,eccf    ,coszrs  ,scon    ,solin   ,solcon , &
                    rei     ,tauxcl  ,tauxci  ,eccf    ,coszrs  ,scon    ,solin   ,solcon , &
                    asdir   ,asdif   ,aldir   ,aldif   ,nmxrgn  , &
                    pmxrgn  ,qrs     ,fsnt    ,fsntc   ,fsntoa  , &
                    fsntoac ,fsnirt  ,fsnrtc  ,fsnirtsq,fsns    , &
                    fsnsc   ,fsdsc   ,fsds    ,sols    ,soll    , &
                    solsd   ,solld   ,frc_day ,                   &
                    fsup    ,fsupc   ,fsdn    ,fsdnc   ,          &
                    aertau  ,aerssa  ,aerasm  ,aerfwd             )
!     call t_stopf('radcswmx')

! -- tls ---------------------------------------------------------------2
!
! convert units of shortwave fields needed by rest of model from cgs to mks
!
      do i=1,ncol
         solin(i) = solin(i)*1.e-3
         fsds(i)  = fsds(i)*1.e-3
         fsnirt(i)= fsnirt(i)*1.e-3
         fsnrtc(i)= fsnrtc(i)*1.e-3
         fsnirtsq(i)= fsnirtsq(i)*1.e-3
         fsnt(i)  = fsnt(i) *1.e-3
         fsns(i)  = fsns(i) *1.e-3
         fsntc(i) = fsntc(i)*1.e-3
         fsnsc(i) = fsnsc(i)*1.e-3
         fsdsc(i) = fsdsc(i)*1.e-3
         fsntoa(i)=fsntoa(i)*1.e-3
         fsntoac(i)=fsntoac(i)*1.e-3
         swcf(i)  = fsntoa(i) - fsntoac(i)
      end do
      ftem(:ncol,:pver) = qrs(:ncol,:pver)/cpair

! added upward/downward total and clear sky fluxes
         do k = 1, pverp
            do i = 1, ncol
            fsup(i,k)  = fsup(i,k)*1.e-3
            fsupc(i,k) = fsupc(i,k)*1.e-3
            fsdn(i,k)  = fsdn(i,k)*1.e-3
            fsdnc(i,k) = fsdnc(i,k)*1.e-3
            end do
         end do

!
! dump shortwave radiation information to history tape buffer (diagnostics)
!

!     call outfld('frc_day ', frc_day, pcols, lchnk)
!     call outfld('sulod_v ', aertau(:,idxvis,1) ,pcols,lchnk)
!     call outfld('ssltod_v', aertau(:,idxvis,2) ,pcols,lchnk)
!     call outfld('carod_v ', aertau(:,idxvis,3) ,pcols,lchnk)
!     call outfld('dustod_v', aertau(:,idxvis,4) ,pcols,lchnk)
!     call outfld('bgod_v  ', aertau(:,idxvis,5) ,pcols,lchnk)
!     call outfld('volcod_v', aertau(:,idxvis,6) ,pcols,lchnk)
!     call outfld('aerod_v ', aertau(:,idxvis,7) ,pcols,lchnk)
!     call outfld('aerssa_v', aerssa(:,idxvis,7) ,pcols,lchnk)
!     call outfld('aerasm_v', aerasm(:,idxvis,7) ,pcols,lchnk)
!     call outfld('aerfwd_v', aerfwd(:,idxvis,7) ,pcols,lchnk)
!     call aerosol_diagnostics (lchnk, ncol, pdel, aerosol)

!     call outfld('qrs     ',ftem  ,pcols,lchnk)
!     call outfld('solin   ',solin ,pcols,lchnk)
!     call outfld('fsds    ',fsds  ,pcols,lchnk)
!     call outfld('fsnirtoa',fsnirt,pcols,lchnk)
!     call outfld('fsnrtoac',fsnrtc,pcols,lchnk)
!     call outfld('fsnrtoas',fsnirtsq,pcols,lchnk)
!     call outfld('fsnt    ',fsnt  ,pcols,lchnk)
!     call outfld('fsns    ',fsns  ,pcols,lchnk)
!     call outfld('fsntc   ',fsntc ,pcols,lchnk)
!     call outfld('fsnsc   ',fsnsc ,pcols,lchnk)
!     call outfld('fsdsc   ',fsdsc ,pcols,lchnk)
!     call outfld('fsntoa  ',fsntoa,pcols,lchnk)
!     call outfld('fsntoac ',fsntoac,pcols,lchnk)
!     call outfld('sols    ',sols  ,pcols,lchnk)
!     call outfld('soll    ',soll  ,pcols,lchnk)
!     call outfld('solsd   ',solsd ,pcols,lchnk)
!     call outfld('solld   ',solld ,pcols,lchnk)

   end if
!
! longwave radiation computation
!
   if (dolw) then

      call get_int_scales(scales)

      call get_aerosol(lchnk, julday, julian, dt, gmt, xtime, m_psp, m_psn, aerosoljp, aerosoljn, &
             m_hybi, paerlev, naer, pint, pcols, pver, pverp, pverr, pverrp, aerosol, scales)

!
! convert upward longwave flux units to cgs
!
      do i=1,ncol
!        lwupcgs(i) = lwup(i)*1000.
         lwupcgs(i) = lwups(i)
      end do
!
! do longwave computation. if not implementing greenhouse gas code then
! first specify trace gas mixing ratios. if greenhouse gas code then:
!  o ixtrcg   => indx of advected n2o tracer
!  o ixtrcg+1 => indx of advected ch4 tracer
!  o ixtrcg+2 => indx of advected cfc11 tracer
!  o ixtrcg+3 => indx of advected cfc12 tracer
!
      if (trace_gas) then
!        call cnst_get_ind('n2o'  , in2o)
!        call cnst_get_ind('ch4'  , ich4)
!        call cnst_get_ind('cfc11', if11)
!        call cnst_get_ind('cfc12', if12)
!        call t_startf("radclwmx")
         call radclwmx(lchnk   ,ncol    ,pcols, pver, pverp ,        & 
                       lwupcgs ,t       ,qm1(1,1,1)       ,o3vmr ,   &
                       pbr     ,pnm     ,pmln    ,piln    ,          &
                       qm1(1,1,in2o)    ,qm1(1,1,ich4)    ,          &
                       qm1(1,1,if11)    ,qm1(1,1,if12)    ,          &
                       cld     ,emis    ,pmxrgn  ,nmxrgn  ,qrl     , &
                       doabsems, abstot, absnxt, emstot,             &
                       flns    ,flnt    ,flnsc   ,flntc   ,flwds   , &
                       flut    ,flutc   ,                            &
                       flup    ,flupc   ,fldn    ,fldnc   ,          &
                       aerosol(:,:,idxvolc))
!        call t_stopf("radclwmx")
      else
         call trcmix(lchnk   ,ncol    ,pcols, pver,  &
                     pmid    ,clat, n2o     ,ch4     ,                     &
                     cfc11   ,cfc12   )

!        call t_startf("radclwmx")
         call radclwmx(lchnk     ,ncol    ,pcols, pver, pverp ,        &
                       lwupcgs   ,t       ,qm1(1,1,1)       ,o3vmr ,   &
                       pbr       ,pnm     ,pmln    ,piln    ,          &
                       n2o       ,ch4     ,cfc11   ,cfc12   ,          &
                       cld       ,emis    ,pmxrgn  ,nmxrgn  ,qrl     , &
                       doabsems, abstot, absnxt, emstot,             &
                       flns      ,flnt    ,flnsc   ,flntc   ,flwds   , &
                       flut      ,flutc   ,                            &
                       flup      ,flupc   ,fldn    ,fldnc   ,          &
                       aerosol(:,:,idxvolc))
!        call t_stopf("radclwmx")
      endif
!
! convert units of longwave fields needed by rest of model from cgs to mks
!
      do i=1,ncol
         flnt(i)  = flnt(i)*1.e-3
         flut(i)  = flut(i)*1.e-3
         flutc(i) = flutc(i)*1.e-3
         flns(i)  = flns(i)*1.e-3
         flntc(i) = flntc(i)*1.e-3
         flnsc(i) = flnsc(i)*1.e-3
         flwds(i) = flwds(i)*1.e-3
         lwcf(i)  = flutc(i) - flut(i)
      end do

! added upward/downward total and clear sky fluxes
         do k = 1, pverp
            do i = 1, ncol
            flup(i,k)  = flup(i,k)*1.e-3
            flupc(i,k) = flupc(i,k)*1.e-3
            fldn(i,k)  = fldn(i,k)*1.e-3
            fldnc(i,k) = fldnc(i,k)*1.e-3
            end do
         end do
!
! dump longwave radiation information to history tape buffer (diagnostics)
!
!     call outfld('qrl     ',qrl(:ncol,:)/cpair,ncol,lchnk)
!     call outfld('flnt    ',flnt  ,pcols,lchnk)
!     call outfld('flut    ',flut  ,pcols,lchnk)
!     call outfld('flutc   ',flutc ,pcols,lchnk)
!     call outfld('flntc   ',flntc ,pcols,lchnk)
!     call outfld('flns    ',flns  ,pcols,lchnk)
!     call outfld('flnsc   ',flnsc ,pcols,lchnk)
!     call outfld('lwcf    ',lwcf  ,pcols,lchnk)
!     call outfld('swcf    ',swcf  ,pcols,lchnk)
!
   end if
!
   return
end subroutine radctl
  subroutine param_cldoptics_calc(ncol, pcols, pver, pverp, pverr, pverrp, ppcnst, &
                                  q, cldn, landfrac, landm,icefrac, &
        pdel,  t, ps, pmid, pint, cicewp, cliqwp, emis, rel, rei, pmxrgn, nmxrgn, snowh )
!
! compute (liquid+ice) water path and cloud water/ice diagnostics
! *** soon this code will compute liquid and ice paths from input liquid and ice mixing ratios
! 
! **** mixes interface and physics code temporarily
!-----------------------------------------------------------------------
!   use physics_types, only: physics_state
!   use history,       only: outfld
!   use pkg_cldoptics, only: cldefr, cldems, cldovrlap, cldclw

    implicit none

! arguments
    integer, intent(in) :: ncol, pcols, pver, pverp, pverr, pverrp, ppcnst
    real(r8), intent(in)  :: q(pcols,pver,ppcnst)     ! moisture arrays
    real(r8), intent(in)  :: cldn(pcols,pver)        ! new cloud fraction
    real(r8), intent(in)  :: pdel(pcols,pver)        ! pressure thickness
    real(r8), intent(in)  :: t(pcols,pver)           ! temperature
    real(r8), intent(in)  :: pmid(pcols,pver)        ! pressure 
    real(r8), intent(in)  :: pint(pcols,pverp)       ! pressure 
    real(r8), intent(in)  :: ps(pcols)               ! surface pressure 
    real(r8), intent(in)  :: landfrac(pcols)         ! land fraction
    real(r8), intent(in)  :: icefrac(pcols)          ! ice fraction
    real(r8), intent(in)  :: landm(pcols)            ! land fraction ramped
    real(r8), intent(in) :: snowh(pcols)         ! snow depth over land, water equivalent (m)

!!$    real(r8), intent(out) :: cwp   (pcols,pver)      ! in-cloud cloud (total) water path
    real(r8), intent(out) :: cicewp(pcols,pver)      ! in-cloud cloud ice water path
    real(r8), intent(out) :: cliqwp(pcols,pver)      ! in-cloud cloud liquid water path
    real(r8), intent(out) :: emis  (pcols,pver)      ! cloud emissivity
    real(r8), intent(out) :: rel   (pcols,pver)      ! effective drop radius (microns)
    real(r8), intent(out) :: rei   (pcols,pver)      ! ice effective drop size (microns)
    real(r8), intent(out) :: pmxrgn(pcols,pver+1)    ! maximum values of pressure for each
    integer , intent(out) :: nmxrgn(pcols)           ! number of maximally overlapped regions

! local variables
    real(r8) :: cwp   (pcols,pver)      ! in-cloud cloud (total) water path
!!$    real(r8) :: cicewp(pcols,pver)      ! in-cloud cloud ice water path
!!$    real(r8) :: cliqwp(pcols,pver)      ! in-cloud cloud liquid water path
    real(r8) :: effcld(pcols,pver)                   ! effective cloud=cld*emis
    real(r8) :: gicewp(pcols,pver)                   ! grid-box cloud ice water path
    real(r8) :: gliqwp(pcols,pver)                   ! grid-box cloud liquid water path
    real(r8) :: gwp   (pcols,pver)                   ! grid-box cloud (total) water path
    real(r8) :: hl     (pcols)                       ! liquid water scale height
    real(r8) :: tgicewp(pcols)                       ! vertically integrated ice water path
    real(r8) :: tgliqwp(pcols)                       ! vertically integrated liquid water path
    real(r8) :: tgwp   (pcols)                       ! vertically integrated (total) cloud water path
    real(r8) :: tpw    (pcols)                       ! total precipitable water
    real(r8) :: clwpold(pcols,pver)                  ! presribed cloud liq. h2o path
    real(r8) :: ficemr (pcols,pver)                  ! ice fraction from ice and liquid mixing ratios

    real(r8) :: rgrav                ! inverse gravitational acceleration

    integer :: i,k                                   ! loop indexes
    integer :: lchnk

!-----------------------------------------------------------------------

! compute liquid and ice water paths
    tgicewp(:ncol) = 0.
    tgliqwp(:ncol) = 0.
    do k=1,pver
       do i = 1,ncol
          gicewp(i,k) = q(i,k,ixcldice)*pdel(i,k)/gravmks*1000.0  ! grid box ice water path.
          gliqwp(i,k) = q(i,k,ixcldliq)*pdel(i,k)/gravmks*1000.0  ! grid box liquid water path.
!!$          gwp   (i,k) = gicewp(i,k) + gliqwp(i,k)
          cicewp(i,k) = gicewp(i,k) / max(0.01_r8,cldn(i,k))                 ! in-cloud ice water path.
          cliqwp(i,k) = gliqwp(i,k) / max(0.01_r8,cldn(i,k))                 ! in-cloud liquid water path.
!!$          cwp   (i,k) = gwp   (i,k) / max(0.01_r8,cldn(i,k))
          ficemr(i,k) = q(i,k,ixcldice) /                 &
               max(1.e-10_r8,(q(i,k,ixcldice)+q(i,k,ixcldliq)))
          
          tgicewp(i)  = tgicewp(i) + gicewp(i,k)
          tgliqwp(i)  = tgliqwp(i) + gliqwp(i,k)
       end do
    end do
    tgwp(:ncol) = tgicewp(:ncol) + tgliqwp(:ncol)
    gwp(:ncol,:pver) = gicewp(:ncol,:pver) + gliqwp(:ncol,:pver) 
    cwp(:ncol,:pver) = cicewp(:ncol,:pver) + cliqwp(:ncol,:pver) 

! compute total preciptable water in column (in mm)
    tpw(:ncol) = 0.0
    rgrav = 1.0/gravmks
    do k=1,pver
       do i=1,ncol
          tpw(i) = tpw(i) + pdel(i,k)*q(i,k,1)*rgrav
       end do
    end do

! diagnostic liquid water path (old specified form)
!   call cldclw(lchnk, ncol, pcols, pver, pverp, state%zi, clwpold, tpw, hl)

! cloud water and ice particle sizes
    call cldefr(lchnk, ncol, pcols, pver, pverp, landfrac, t, rel, rei, ps, pmid, landm, icefrac, snowh)

! cloud emissivity.
    call cldems(lchnk, ncol, pcols, pver, pverp, cwp, ficemr, rei, emis)

! effective cloud cover
    do k=1,pver
       do i=1,ncol
          effcld(i,k) = cldn(i,k)*emis(i,k)
       end do
    end do

! determine parameters for maximum/random overlap
    call cldovrlap(lchnk, ncol, pcols, pver, pverp, pint, cldn, nmxrgn, pmxrgn)

!   call outfld('gcldlwp' ,gwp    , pcols,lchnk)
!   call outfld('tgcldcwp',tgwp   , pcols,lchnk)
!   call outfld('tgcldlwp',tgliqwp, pcols,lchnk)
!   call outfld('tgcldiwp',tgicewp, pcols,lchnk)
!   call outfld('icldlwp' ,cwp    , pcols,lchnk)
!   call outfld('setlwp'  ,clwpold, pcols,lchnk)
!   call outfld('effcld'  ,effcld , pcols,lchnk)
!   call outfld('lwsh'    ,hl     , pcols,lchnk)

  end subroutine param_cldoptics_calc

subroutine radabs(lchnk   ,ncol    ,pcols, pver, pverp,   &
   pbr    ,pnm     ,co2em    ,co2eml  ,tplnka  , &
   s2c    ,tcg     ,w        ,h2otr   ,plco2   , &
   plh2o  ,co2t    ,tint     ,tlayr   ,plol    , &
   plos   ,pmln    ,piln     ,ucfc11  ,ucfc12  , &
   un2o0  ,un2o1   ,uch4     ,uco211  ,uco212  , &
   uco213 ,uco221  ,uco222   ,uco223  ,uptype  , &
   bn2o0  ,bn2o1   ,bch4    ,abplnk1  ,abplnk2 , &
   abstot ,absnxt  ,plh2ob  ,wb       , &
   aer_mpp ,aer_trn_ttl)
!----------------------------------------------------------------------- 
! 
! purpose: 
! compute absorptivities for h2o, co2, o3, ch4, n2o, cfc11 and cfc12
! 
! method: 
! h2o  ....  uses nonisothermal emissivity method for water vapor from
!            ramanathan, v. and  p.downey, 1986: a nonisothermal
!            emissivity and absorptivity formulation for water vapor
!            journal of geophysical research, vol. 91., d8, pp 8649-8666
!
!            implementation updated by collins, hackney, and edwards (2001)
!               using line-by-line calculations based upon hitran 1996 and
!               ckd 2.1 for absorptivity and emissivity
!
!            implementation updated by collins, lee-taylor, and edwards (2003)
!               using line-by-line calculations based upon hitran 2000 and
!               ckd 2.4 for absorptivity and emissivity
!
! co2  ....  uses absorptance parameterization of the 15 micro-meter
!            (500 - 800 cm-1) band system of carbon dioxide, from
!            kiehl, j.t. and b.p.briegleb, 1991: a new parameterization
!            of the absorptance due to the 15 micro-meter band system
!            of carbon dioxide jouranl of geophysical research,
!            vol. 96., d5, pp 9013-9019.
!            parameterizations for the 9.4 and 10.4 mircon bands of co2
!            are also included.
!
! o3   ....  uses absorptance parameterization of the 9.6 micro-meter
!            band system of ozone, from ramanathan, v. and r.dickinson,
!            1979: the role of stratospheric ozone in the zonal and
!            seasonal radiative energy balance of the earth-troposphere
!            system. journal of the atmospheric sciences, vol. 36,
!            pp 1084-1104
!
! ch4  ....  uses a broad band model for the 7.7 micron band of methane.
!
! n20  ....  uses a broad band model for the 7.8, 8.6 and 17.0 micron
!            bands of nitrous oxide
!
! cfc11 ...  uses a quasi-linear model for the 9.2, 10.7, 11.8 and 12.5
!            micron bands of cfc11
!
! cfc12 ...  uses a quasi-linear model for the 8.6, 9.1, 10.8 and 11.2
!            micron bands of cfc12
!
!
! computes individual absorptivities for non-adjacent layers, accounting
! for band overlap, and sums to obtain the total; then, computes the
! nearest layer contribution.
! 
! author: w. collins (h2o absorptivity) and j. kiehl
! 
!-----------------------------------------------------------------------
!------------------------------arguments--------------------------------
!
! input arguments
!
   integer, intent(in) :: lchnk                       ! chunk identifier
   integer, intent(in) :: ncol                        ! number of atmospheric columns
   integer, intent(in) :: pcols, pver, pverp

   real(r8), intent(in) :: pbr(pcols,pver)            ! prssr at mid-levels (dynes/cm2)
   real(r8), intent(in) :: pnm(pcols,pverp)           ! prssr at interfaces (dynes/cm2)
   real(r8), intent(in) :: co2em(pcols,pverp)         ! co2 emissivity function
   real(r8), intent(in) :: co2eml(pcols,pver)         ! co2 emissivity function
   real(r8), intent(in) :: tplnka(pcols,pverp)        ! planck fnctn level temperature
   real(r8), intent(in) :: s2c(pcols,pverp)           ! h2o continuum path length
   real(r8), intent(in) :: tcg(pcols,pverp)           ! h2o-mass-wgted temp. (curtis-godson approx.)
   real(r8), intent(in) :: w(pcols,pverp)             ! h2o prs wghted path
   real(r8), intent(in) :: h2otr(pcols,pverp)         ! h2o trnsmssn fnct for o3 overlap
   real(r8), intent(in) :: plco2(pcols,pverp)         ! co2 prs wghted path length
   real(r8), intent(in) :: plh2o(pcols,pverp)         ! h2o prs wfhted path length
   real(r8), intent(in) :: co2t(pcols,pverp)          ! tmp and prs wghted path length
   real(r8), intent(in) :: tint(pcols,pverp)          ! interface temperatures
   real(r8), intent(in) :: tlayr(pcols,pverp)         ! k-1 level temperatures
   real(r8), intent(in) :: plol(pcols,pverp)          ! ozone prs wghted path length
   real(r8), intent(in) :: plos(pcols,pverp)          ! ozone path length
   real(r8), intent(in) :: pmln(pcols,pver)           ! ln(pmidm1)
   real(r8), intent(in) :: piln(pcols,pverp)          ! ln(pintm1)
   real(r8), intent(in) :: plh2ob(nbands,pcols,pverp) ! pressure weighted h2o path with 
                                                      !    hulst-curtis-godson temp. factor 
                                                      !    for h2o bands 
   real(r8), intent(in) :: wb(nbands,pcols,pverp)     ! h2o path length with 
                                                      !    hulst-curtis-godson temp. factor 
                                                      !    for h2o bands 

   real(r8), intent(in) :: aer_mpp(pcols,pverp) ! straer path above kth interface level
   real(r8), intent(in) :: aer_trn_ttl(pcols,pverp,pverp,bnd_nbr_lw) ! aer trn.


!
! trace gas variables
!
   real(r8), intent(in) :: ucfc11(pcols,pverp)        ! cfc11 path length
   real(r8), intent(in) :: ucfc12(pcols,pverp)        ! cfc12 path length
   real(r8), intent(in) :: un2o0(pcols,pverp)         ! n2o path length
   real(r8), intent(in) :: un2o1(pcols,pverp)         ! n2o path length (hot band)
   real(r8), intent(in) :: uch4(pcols,pverp)          ! ch4 path length
   real(r8), intent(in) :: uco211(pcols,pverp)        ! co2 9.4 micron band path length
   real(r8), intent(in) :: uco212(pcols,pverp)        ! co2 9.4 micron band path length
   real(r8), intent(in) :: uco213(pcols,pverp)        ! co2 9.4 micron band path length
   real(r8), intent(in) :: uco221(pcols,pverp)        ! co2 10.4 micron band path length
   real(r8), intent(in) :: uco222(pcols,pverp)        ! co2 10.4 micron band path length
   real(r8), intent(in) :: uco223(pcols,pverp)        ! co2 10.4 micron band path length
   real(r8), intent(in) :: uptype(pcols,pverp)        ! continuum path length
   real(r8), intent(in) :: bn2o0(pcols,pverp)         ! pressure factor for n2o
   real(r8), intent(in) :: bn2o1(pcols,pverp)         ! pressure factor for n2o
   real(r8), intent(in) :: bch4(pcols,pverp)          ! pressure factor for ch4
   real(r8), intent(in) :: abplnk1(14,pcols,pverp)    ! non-nearest layer planck factor
   real(r8), intent(in) :: abplnk2(14,pcols,pverp)    ! nearest layer factor
!
! output arguments
!
   real(r8), intent(out) :: abstot(pcols,pverp,pverp) ! total absorptivity
   real(r8), intent(out) :: absnxt(pcols,pver,4)      ! total nearest layer absorptivity
!
!---------------------------local variables-----------------------------
!
   integer i                   ! longitude index
   integer k                   ! level index
   integer k1                  ! level index
   integer k2                  ! level index
   integer kn                  ! nearest level index
   integer wvl                 ! wavelength index

   real(r8) abstrc(pcols)              ! total trace gas absorptivity
   real(r8) bplnk(14,pcols,4)          ! planck functions for sub-divided layers
   real(r8) pnew(pcols)        ! effective pressure for h2o vapor linewidth
   real(r8) pnewb(nbands)      ! effective pressure for h2o linewidth w/
                               !    hulst-curtis-godson correction for
                               !    each band
   real(r8) u(pcols)           ! pressure weighted h2o path length
   real(r8) ub(nbands)         ! pressure weighted h2o path length with
                               !    hulst-curtis-godson correction for
                               !    each band
   real(r8) tbar(pcols,4)      ! mean layer temperature
   real(r8) emm(pcols,4)       ! mean co2 emissivity
   real(r8) o3emm(pcols,4)     ! mean o3 emissivity
   real(r8) o3bndi             ! ozone band parameter
   real(r8) temh2o(pcols,4)    ! mean layer temperature equivalent to tbar
   real(r8) k21                ! exponential coefficient used to calculate
!                              !  rotation band transmissvty in the 650-800
!                              !  cm-1 region (tr1)
   real(r8) k22                ! exponential coefficient used to calculate
!                              !  rotation band transmissvty in the 500-650
!                              !  cm-1 region (tr2)
   real(r8) uc1(pcols)         ! h2o continuum pathlength in 500-800 cm-1
   real(r8) to3h2o(pcols)      ! h2o trnsmsn for overlap with o3
   real(r8) pi                 ! for co2 absorptivity computation
   real(r8) sqti(pcols)        ! used to store sqrt of mean temperature
   real(r8) et                 ! co2 hot band factor
   real(r8) et2                ! co2 hot band factor squared
   real(r8) et4                ! co2 hot band factor to fourth power
   real(r8) omet               ! co2 stimulated emission term
   real(r8) f1co2              ! co2 central band factor
   real(r8) f2co2(pcols)       ! co2 weak band factor
   real(r8) f3co2(pcols)       ! co2 weak band factor
   real(r8) t1co2(pcols)       ! overlap factr weak bands on strong band
   real(r8) sqwp               ! sqrt of co2 pathlength
   real(r8) f1sqwp(pcols)      ! main co2 band factor
   real(r8) oneme              ! co2 stimulated emission term
   real(r8) alphat             ! part of the co2 stimulated emission term
   real(r8) wco2               ! constants used to define co2 pathlength
   real(r8) posqt              ! effective pressure for co2 line width
   real(r8) u7(pcols)          ! co2 hot band path length
   real(r8) u8                 ! co2 hot band path length
   real(r8) u9                 ! co2 hot band path length
   real(r8) u13                ! co2 hot band path length
   real(r8) rbeta7(pcols)      ! inverse of co2 hot band line width par
   real(r8) rbeta8             ! inverse of co2 hot band line width par
   real(r8) rbeta9             ! inverse of co2 hot band line width par
   real(r8) rbeta13            ! inverse of co2 hot band line width par
   real(r8) tpatha             ! for absorptivity computation
   real(r8) abso(pcols,4)      ! absorptivity for various gases/bands
   real(r8) dtx(pcols)         ! planck temperature minus 250 k
   real(r8) dty(pcols)         ! path temperature minus 250 k
   real(r8) term7(pcols,2)     ! kl_inf(i) in eq(r8) of table a3a of r&d
   real(r8) term8(pcols,2)     ! delta kl_inf(i) in eq(r8)
   real(r8) tr1                ! eqn(6) in table a2 of r&d for 650-800
   real(r8) tr10(pcols)        ! eqn (6) times eq(4) in table a2
!                              !  of r&d for 500-650 cm-1 region
   real(r8) tr2                ! eqn(6) in table a2 of r&d for 500-650
   real(r8) tr5                ! eqn(4) in table a2 of r&d for 650-800
   real(r8) tr6                ! eqn(4) in table a2 of r&d for 500-650
   real(r8) tr9(pcols)         ! equation (6) times eq(4) in table a2
!                              !  of r&d for 650-800 cm-1 region
   real(r8) sqrtu(pcols)       ! sqrt of pressure weighted h20 pathlength
   real(r8) fwk(pcols)         ! equation(33) in r&d far wing correction
   real(r8) fwku(pcols)        ! gu term in eqs(1) and (6) in table a2
   real(r8) to3co2(pcols)      ! p weighted temp in ozone band model
   real(r8) dpnm(pcols)        ! pressure difference between two levels
   real(r8) pnmsq(pcols,pverp) ! pressure squared
   real(r8) dw(pcols)          ! amount of h2o between two levels
   real(r8) uinpl(pcols,4)     ! nearest layer subdivision factor
   real(r8) winpl(pcols,4)     ! nearest layer subdivision factor
   real(r8) zinpl(pcols,4)     ! nearest layer subdivision factor
   real(r8) pinpl(pcols,4)     ! nearest layer subdivision factor
   real(r8) dplh2o(pcols)      ! difference in press weighted h2o amount
   real(r8) r293               ! 1/293
   real(r8) r250               ! 1/250
   real(r8) r3205              ! line width factor for o3 (see r&di)
   real(r8) r300               ! 1/300
   real(r8) rsslp              ! reciprocal of sea level pressure
   real(r8) r2sslp             ! 1/2 of rsslp
   real(r8) ds2c               ! y in eq(7) in table a2 of r&d
   real(r8)  dplos             ! ozone pathlength eq(a2) in r&di
   real(r8) dplol              ! presure weighted ozone pathlength
   real(r8) tlocal             ! local interface temperature
   real(r8) beta               ! ozone mean line parameter eq(a3) in r&di
!                               (includes voigt line correction factor)
   real(r8) rphat              ! effective pressure for ozone beta
   real(r8) tcrfac             ! ozone temperature factor table 1 r&di
   real(r8) tmp1               ! ozone band factor see eq(a1) in r&di
   real(r8) u1                 ! effective ozone pathlength eq(a2) in r&di
   real(r8) realnu             ! 1/beta factor in ozone band model eq(a1)
   real(r8) tmp2               ! ozone band factor see eq(a1) in r&di
   real(r8) u2                 ! effective ozone pathlength eq(a2) in r&di
   real(r8) rsqti              ! reciprocal of sqrt of path temperature
   real(r8) tpath              ! path temperature used in co2 band model
   real(r8) tmp3               ! weak band factor see k&b
   real(r8) rdpnmsq            ! reciprocal of difference in press^2
   real(r8) rdpnm              ! reciprocal of difference in press
   real(r8) p1                 ! mean pressure factor
   real(r8) p2                 ! mean pressure factor
   real(r8) dtym10             ! t - 260 used in eq(9) and (10) table a3a
   real(r8) dplco2             ! co2 path length
   real(r8) te                 ! a_0 t factor in ozone model table 1 of r&di
   real(r8) denom              ! denominator in eq(r8) of table a3a of r&d
   real(r8) th2o(pcols)        ! transmission due to h2o
   real(r8) tco2(pcols)        ! transmission due to co2
   real(r8) to3(pcols)         ! transmission due to o3
!
! transmission terms for various spectral intervals:
!
   real(r8) trab2(pcols)       ! h2o   500 -  800 cm-1
   real(r8) absbnd             ! proportional to co2 band absorptance
   real(r8) dbvtit(pcols,pverp)! intrfc drvtv plnck fnctn for o3
   real(r8) dbvtly(pcols,pver) ! level drvtv plnck fnctn for o3
!
! variables for collins/hackney/edwards (c/h/e) & 
!       collins/lee-taylor/edwards (c/lt/e) h2o parameterization

!
! notation:
! u   = integral (p/p_0 dw)  eq. 15 in ramanathan/downey 1986
! p   = atmospheric pressure
! p_0 = reference atmospheric pressure
! w   = precipitable water path
! t_e = emission temperature
! t_p = path temperature
! rh  = path relative humidity
!
   real(r8) fa               ! asymptotic value of abs. as u->infinity
   real(r8) a_star           ! normalized absorptivity for non-window
   real(r8) l_star           ! interpolated line transmission
   real(r8) c_star           ! interpolated continuum transmission

   real(r8) te1              ! emission temperature
   real(r8) te2              ! te^2
   real(r8) te3              ! te^3
   real(r8) te4              ! te^4
   real(r8) te5              ! te^5

   real(r8) log_u            ! log base 10 of u 
   real(r8) log_uc           ! log base 10 of h2o continuum path
   real(r8) log_p            ! log base 10 of p
   real(r8) t_p              ! t_p
   real(r8) t_e              ! t_e (offset by t_p)

   integer iu                ! index for log10(u)
   integer iu1               ! iu + 1
   integer iuc               ! index for log10(h2o continuum path)
   integer iuc1              ! iuc + 1
   integer ip                ! index for log10(p)
   integer ip1               ! ip + 1
   integer itp               ! index for t_p
   integer itp1              ! itp + 1
   integer ite               ! index for t_e
   integer ite1              ! ite + 1
   integer irh               ! index for rh
   integer irh1              ! irh + 1

   real(r8) dvar             ! normalized variation in t_p/t_e/p/u
   real(r8) uvar             ! u * diffusivity factor
   real(r8) uscl             ! factor for lineary scaling as u->0

   real(r8) wu               ! weight for u
   real(r8) wu1              ! 1 - wu
   real(r8) wuc              ! weight for h2o continuum path
   real(r8) wuc1             ! 1 - wuc
   real(r8) wp               ! weight for p
   real(r8) wp1              ! 1 - wp
   real(r8) wtp              ! weight for t_p
   real(r8) wtp1             ! 1 - wtp
   real(r8) wte              ! weight for t_e
   real(r8) wte1             ! 1 - wte
   real(r8) wrh              ! weight for rh
   real(r8) wrh1             ! 1 - wrh

   real(r8) w_0_0_           ! weight for tp/te combination
   real(r8) w_0_1_           ! weight for tp/te combination
   real(r8) w_1_0_           ! weight for tp/te combination
   real(r8) w_1_1_           ! weight for tp/te combination

   real(r8) w_0_00           ! weight for tp/te/rh combination
   real(r8) w_0_01           ! weight for tp/te/rh combination
   real(r8) w_0_10           ! weight for tp/te/rh combination
   real(r8) w_0_11           ! weight for tp/te/rh combination
   real(r8) w_1_00           ! weight for tp/te/rh combination
   real(r8) w_1_01           ! weight for tp/te/rh combination
   real(r8) w_1_10           ! weight for tp/te/rh combination
   real(r8) w_1_11           ! weight for tp/te/rh combination

   real(r8) w00_00           ! weight for p/tp/te/rh combination
   real(r8) w00_01           ! weight for p/tp/te/rh combination
   real(r8) w00_10           ! weight for p/tp/te/rh combination
   real(r8) w00_11           ! weight for p/tp/te/rh combination
   real(r8) w01_00           ! weight for p/tp/te/rh combination
   real(r8) w01_01           ! weight for p/tp/te/rh combination
   real(r8) w01_10           ! weight for p/tp/te/rh combination
   real(r8) w01_11           ! weight for p/tp/te/rh combination
   real(r8) w10_00           ! weight for p/tp/te/rh combination
   real(r8) w10_01           ! weight for p/tp/te/rh combination
   real(r8) w10_10           ! weight for p/tp/te/rh combination
   real(r8) w10_11           ! weight for p/tp/te/rh combination
   real(r8) w11_00           ! weight for p/tp/te/rh combination
   real(r8) w11_01           ! weight for p/tp/te/rh combination
   real(r8) w11_10           ! weight for p/tp/te/rh combination
   real(r8) w11_11           ! weight for p/tp/te/rh combination

   integer ib                ! spectral interval:
                             !   1 = 0-800 cm^-1 and 1200-2200 cm^-1
                             !   2 = 800-1200 cm^-1


   real(r8) pch2o            ! h2o continuum path
   real(r8) fch2o            ! temp. factor for continuum
   real(r8) uch2o            ! u corresponding to h2o cont. path (window)

   real(r8) fdif             ! secant(zenith angle) for diffusivity approx.

   real(r8) sslp_mks         ! sea-level pressure in mks units
   real(r8) esx              ! saturation vapor pressure returned by vqsatd
   real(r8) qsx              ! saturation mixing ratio returned by vqsatd
   real(r8) pnew_mks         ! pnew in mks units
   real(r8) q_path           ! effective specific humidity along path
   real(r8) rh_path          ! effective relative humidity along path
   real(r8) omeps            ! 1 - epsilo

   integer  iest             ! index in estblh2o

      integer bnd_idx        ! lw band index
      real(r8) aer_pth_dlt   ! [kg m-2] straer path between interface levels k1 and k2
      real(r8) aer_pth_ngh(pcols)
                             ! [kg m-2] straer path between neighboring layers
      real(r8) odap_aer_ttl  ! [fraction] total path absorption optical depth
      real(r8) aer_trn_ngh(pcols,bnd_nbr_lw) 
                             ! [fraction] total transmission between 
                             !            nearest neighbor sub-levels
!
!--------------------------statement function---------------------------
!
   real(r8) dbvt,t             ! planck fnctn tmp derivative for o3
!
   dbvt(t)=(-2.8911366682e-4+(2.3771251896e-6+1.1305188929e-10*t)*t)/ &
      (1.0+(-6.1364820707e-3+1.5550319767e-5*t)*t)
!
!
!-----------------------------------------------------------------------
!
! initialize
!
   do k2=1,ntoplw-1
      do k1=1,ntoplw-1
         abstot(:,k1,k2) = inf    ! set unused portions for lf95 restart write
      end do
   end do
   do k2=1,4
      do k1=1,ntoplw-1
         absnxt(:,k1,k2) = inf    ! set unused portions for lf95 restart write
      end do
   end do

   do k=ntoplw,pverp
      abstot(:,k,k) = inf         ! set unused portions for lf95 restart write
   end do

   do k=ntoplw,pver
      do i=1,ncol
         dbvtly(i,k) = dbvt(tlayr(i,k+1))
         dbvtit(i,k) = dbvt(tint(i,k))
      end do
   end do
   do i=1,ncol
      dbvtit(i,pverp) = dbvt(tint(i,pverp))
   end do
!
   r293    = 1./293.
   r250    = 1./250.
   r3205   = 1./.3205
   r300    = 1./300.
   rsslp   = 1./sslp
   r2sslp  = 1./(2.*sslp)
!
!constants for computing u corresponding to h2o cont. path
!
   fdif       = 1.66
   sslp_mks   = sslp / 10.0
   omeps      = 1.0 - epsilo
!
! non-adjacent layer absorptivity:
!
! abso(i,1)     0 -  800 cm-1   h2o rotation band
! abso(i,1)  1200 - 2200 cm-1   h2o vibration-rotation band
! abso(i,2)   800 - 1200 cm-1   h2o window
!
! separation between rotation and vibration-rotation dropped, so
!                only 2 slots needed for h2o absorptivity
!
! 500-800 cm^-1 h2o continuum/line overlap already included
!                in abso(i,1).  this used to be in abso(i,4)
!
! abso(i,3)   o3  9.6 micrometer band (nu3 and nu1 bands)
! abso(i,4)   co2 15  micrometer band system
!
   do k=ntoplw,pverp
      do i=1,ncol
         pnmsq(i,k) = pnm(i,k)**2
         dtx(i) = tplnka(i,k) - 250.
      end do
   end do
!
! non-nearest layer level loops
!
   do k1=pverp,ntoplw,-1
      do k2=pverp,ntoplw,-1
         if (k1 == k2) cycle
         do i=1,ncol
            dplh2o(i) = plh2o(i,k1) - plh2o(i,k2)
            u(i)      = abs(dplh2o(i))
            sqrtu(i)  = sqrt(u(i))
            ds2c      = abs(s2c(i,k1) - s2c(i,k2))
            dw(i)     = abs(w(i,k1) - w(i,k2))
            uc1(i)    = (ds2c + 1.7e-3*u(i))*(1. +  2.*ds2c)/(1. + 15.*ds2c)
            pch2o     = ds2c
            pnew(i)   = u(i)/dw(i)
            pnew_mks  = pnew(i) * sslp_mks
!
! changed effective path temperature to std. curtis-godson form
!
            tpatha = abs(tcg(i,k1) - tcg(i,k2))/dw(i)
            t_p = min(max(tpatha, min_tp_h2o), max_tp_h2o)
            iest = floor(t_p) - min_tp_h2o
            esx = estblh2o(iest) + (estblh2o(iest+1)-estblh2o(iest)) * &
                 (t_p - min_tp_h2o - iest)
            qsx = epsilo * esx / (pnew_mks - omeps * esx)
!
! compute effective rh along path
!
            q_path = dw(i) / abs(pnm(i,k1) - pnm(i,k2)) / rga
!
! calculate effective u, pnew for each band using
!        hulst-curtis-godson approximation:
! formulae: goody and yung, atmospheric radiation: theoretical basis, 
!           2nd edition, oxford university press, 1989.
! effective h2o path (w)
!      eq. 6.24, p. 228
! effective h2o path pressure (pnew = u/w):
!      eq. 6.29, p. 228
!
            ub(1) = abs(plh2ob(1,i,k1) - plh2ob(1,i,k2)) / psi(t_p,1)
            ub(2) = abs(plh2ob(2,i,k1) - plh2ob(2,i,k2)) / psi(t_p,2)
            
            pnewb(1) = ub(1) / abs(wb(1,i,k1) - wb(1,i,k2)) * phi(t_p,1)
            pnewb(2) = ub(2) / abs(wb(2,i,k1) - wb(2,i,k2)) * phi(t_p,2)

            dtx(i)      = tplnka(i,k2) - 250.
            dty(i)      = tpatha       - 250.

            fwk(i)  = fwcoef + fwc1/(1. + fwc2*u(i))
            fwku(i) = fwk(i)*u(i)
!
! define variables for c/h/e (now c/lt/e) fit
!
! abso(i,1)     0 -  800 cm-1   h2o rotation band
! abso(i,1)  1200 - 2200 cm-1   h2o vibration-rotation band
! abso(i,2)   800 - 1200 cm-1   h2o window
!
! separation between rotation and vibration-rotation dropped, so
!                only 2 slots needed for h2o absorptivity
!
! notation:
! u   = integral (p/p_0 dw)  
! p   = atmospheric pressure
! p_0 = reference atmospheric pressure
! w   = precipitable water path
! t_e = emission temperature
! t_p = path temperature
! rh  = path relative humidity
!
!
! terms for asymptotic value of emissivity
!
            te1  = tplnka(i,k2)
            te2  = te1 * te1
            te3  = te2 * te1
            te4  = te3 * te1
            te5  = te4 * te1

!
!  band-independent indices for lines and continuum tables
!
            dvar = (t_p - min_tp_h2o) / dtp_h2o
            itp = min(max(int(aint(dvar,r8)) + 1, 1), n_tp - 1)
            itp1 = itp + 1
            wtp = dvar - floor(dvar)
            wtp1 = 1.0 - wtp
            
            t_e = min(max(tplnka(i,k2)-t_p, min_te_h2o), max_te_h2o)
            dvar = (t_e - min_te_h2o) / dte_h2o
            ite = min(max(int(aint(dvar,r8)) + 1, 1), n_te - 1)
            ite1 = ite + 1
            wte = dvar - floor(dvar)
            wte1 = 1.0 - wte
            
            rh_path = min(max(q_path / qsx, min_rh_h2o), max_rh_h2o)
            dvar = (rh_path - min_rh_h2o) / drh_h2o
            irh = min(max(int(aint(dvar,r8)) + 1, 1), n_rh - 1)
            irh1 = irh + 1
            wrh = dvar - floor(dvar)
            wrh1 = 1.0 - wrh

            w_0_0_ = wtp  * wte
            w_0_1_ = wtp  * wte1
            w_1_0_ = wtp1 * wte 
            w_1_1_ = wtp1 * wte1
            
            w_0_00 = w_0_0_ * wrh
            w_0_01 = w_0_0_ * wrh1
            w_0_10 = w_0_1_ * wrh
            w_0_11 = w_0_1_ * wrh1
            w_1_00 = w_1_0_ * wrh
            w_1_01 = w_1_0_ * wrh1
            w_1_10 = w_1_1_ * wrh
            w_1_11 = w_1_1_ * wrh1

!
! h2o continuum path for 0-800 and 1200-2200 cm^-1
!
!    assume foreign continuum dominates total h2o continuum in these bands
!    per clough et al, jgr, v. 97, no. d14 (oct 20, 1992), p. 15776
!    then the effective h2o path is just 
!         u_c = integral[ f(p) dw ]
!    where 
!           w = water-vapor mass and 
!        f(p) = dependence of foreign continuum on pressure 
!             = p / sslp
!    then 
!         u_c = u (the same effective h2o path as for lines)
!
!
! continuum terms for 800-1200 cm^-1
!
!    assume self continuum dominates total h2o continuum for this band
!    per clough et al, jgr, v. 97, no. d14 (oct 20, 1992), p. 15776
!    then the effective h2o self-continuum path is 
!         u_c = integral[ h(e,t) dw ]                        (*eq. 1*)
!    where 
!           w = water-vapor mass and 
!           e = partial pressure of h2o along path
!           t = temperature along path
!      h(e,t) = dependence of foreign continuum on e,t
!             = e / sslp * f(t)
!
!    replacing
!           e =~ q * p / epsilo
!           q = mixing ratio of h2o
!     epsilo = 0.622
!
!    and using the definition
!           u = integral [ (p / sslp) dw ]
!             = (p / sslp) w                                 (homogeneous path)
!
!    the effective path length for the self continuum is
!         u_c = (q / epsilo) f(t) u                         (*eq. 2*)
!
!    once values of t, u, and q have been calculated for the inhomogeneous
!        path, this sets u_c for the corresponding
!        homogeneous atmosphere.  however, this need not equal the
!        value of u_c' defined by eq. 1 for the actual inhomogeneous atmosphere
!        under consideration.
!
!    solution: hold t and q constant, solve for u' that gives u_c' by
!        inverting eq. (2):
!
!        u' = (u_c * epsilo) / (q * f(t))
!
            fch2o = fh2oself(t_p) 
            uch2o = (pch2o * epsilo) / (q_path * fch2o)

!
! band-dependent indices for non-window
!
            ib = 1

            uvar = ub(ib) * fdif
            log_u  = min(log10(max(uvar, min_u_h2o)), max_lu_h2o)
            dvar = (log_u - min_lu_h2o) / dlu_h2o
            iu = min(max(int(aint(dvar,r8)) + 1, 1), n_u - 1)
            iu1 = iu + 1
            wu = dvar - floor(dvar)
            wu1 = 1.0 - wu
            
            log_p  = min(log10(max(pnewb(ib), min_p_h2o)), max_lp_h2o)
            dvar = (log_p - min_lp_h2o) / dlp_h2o
            ip = min(max(int(aint(dvar,r8)) + 1, 1), n_p - 1)
            ip1 = ip + 1
            wp = dvar - floor(dvar)
            wp1 = 1.0 - wp
         
            w00_00 = wp  * w_0_00 
            w00_01 = wp  * w_0_01 
            w00_10 = wp  * w_0_10 
            w00_11 = wp  * w_0_11 
            w01_00 = wp  * w_1_00 
            w01_01 = wp  * w_1_01 
            w01_10 = wp  * w_1_10 
            w01_11 = wp  * w_1_11 
            w10_00 = wp1 * w_0_00 
            w10_01 = wp1 * w_0_01 
            w10_10 = wp1 * w_0_10 
            w10_11 = wp1 * w_0_11 
            w11_00 = wp1 * w_1_00 
            w11_01 = wp1 * w_1_01 
            w11_10 = wp1 * w_1_10 
            w11_11 = wp1 * w_1_11 
!
! asymptotic value of absorptivity as u->infinity
!
            fa = fat(1,ib) + &
                 fat(2,ib) * te1 + &
                 fat(3,ib) * te2 + &
                 fat(4,ib) * te3 + &
                 fat(5,ib) * te4 + &
                 fat(6,ib) * te5

            a_star = &
                 ah2onw(ip , itp , iu , ite , irh ) * w11_11 * wu1 + &
                 ah2onw(ip , itp , iu , ite , irh1) * w11_10 * wu1 + &
                 ah2onw(ip , itp , iu , ite1, irh ) * w11_01 * wu1 + &
                 ah2onw(ip , itp , iu , ite1, irh1) * w11_00 * wu1 + &
                 ah2onw(ip , itp , iu1, ite , irh ) * w11_11 * wu  + &
                 ah2onw(ip , itp , iu1, ite , irh1) * w11_10 * wu  + &
                 ah2onw(ip , itp , iu1, ite1, irh ) * w11_01 * wu  + &
                 ah2onw(ip , itp , iu1, ite1, irh1) * w11_00 * wu  + &
                 ah2onw(ip , itp1, iu , ite , irh ) * w10_11 * wu1 + &
                 ah2onw(ip , itp1, iu , ite , irh1) * w10_10 * wu1 + &
                 ah2onw(ip , itp1, iu , ite1, irh ) * w10_01 * wu1 + &
                 ah2onw(ip , itp1, iu , ite1, irh1) * w10_00 * wu1 + &
                 ah2onw(ip , itp1, iu1, ite , irh ) * w10_11 * wu  + &
                 ah2onw(ip , itp1, iu1, ite , irh1) * w10_10 * wu  + &
                 ah2onw(ip , itp1, iu1, ite1, irh ) * w10_01 * wu  + &
                 ah2onw(ip , itp1, iu1, ite1, irh1) * w10_00 * wu  + &
                 ah2onw(ip1, itp , iu , ite , irh ) * w01_11 * wu1 + &
                 ah2onw(ip1, itp , iu , ite , irh1) * w01_10 * wu1 + &
                 ah2onw(ip1, itp , iu , ite1, irh ) * w01_01 * wu1 + &
                 ah2onw(ip1, itp , iu , ite1, irh1) * w01_00 * wu1 + &
                 ah2onw(ip1, itp , iu1, ite , irh ) * w01_11 * wu  + &
                 ah2onw(ip1, itp , iu1, ite , irh1) * w01_10 * wu  + &
                 ah2onw(ip1, itp , iu1, ite1, irh ) * w01_01 * wu  + &
                 ah2onw(ip1, itp , iu1, ite1, irh1) * w01_00 * wu  + &
                 ah2onw(ip1, itp1, iu , ite , irh ) * w00_11 * wu1 + &
                 ah2onw(ip1, itp1, iu , ite , irh1) * w00_10 * wu1 + &
                 ah2onw(ip1, itp1, iu , ite1, irh ) * w00_01 * wu1 + &
                 ah2onw(ip1, itp1, iu , ite1, irh1) * w00_00 * wu1 + &
                 ah2onw(ip1, itp1, iu1, ite , irh ) * w00_11 * wu  + &
                 ah2onw(ip1, itp1, iu1, ite , irh1) * w00_10 * wu  + &
                 ah2onw(ip1, itp1, iu1, ite1, irh ) * w00_01 * wu  + &
                 ah2onw(ip1, itp1, iu1, ite1, irh1) * w00_00 * wu 
            abso(i,ib) = min(max(fa * (1.0 - (1.0 - a_star) * &
                                 aer_trn_ttl(i,k1,k2,ib)), &
                             0.0_r8), 1.0_r8)
!
! invoke linear limit for scaling wrt u below min_u_h2o
!
            if (uvar < min_u_h2o) then
               uscl = uvar / min_u_h2o
               abso(i,ib) = abso(i,ib) * uscl
            endif
                         
!
! band-dependent indices for window
!
            ib = 2

            uvar = ub(ib) * fdif
            log_u  = min(log10(max(uvar, min_u_h2o)), max_lu_h2o)
            dvar = (log_u - min_lu_h2o) / dlu_h2o
            iu = min(max(int(aint(dvar,r8)) + 1, 1), n_u - 1)
            iu1 = iu + 1
            wu = dvar - floor(dvar)
            wu1 = 1.0 - wu
            
            log_p  = min(log10(max(pnewb(ib), min_p_h2o)), max_lp_h2o)
            dvar = (log_p - min_lp_h2o) / dlp_h2o
            ip = min(max(int(aint(dvar,r8)) + 1, 1), n_p - 1)
            ip1 = ip + 1
            wp = dvar - floor(dvar)
            wp1 = 1.0 - wp
         
            w00_00 = wp  * w_0_00 
            w00_01 = wp  * w_0_01 
            w00_10 = wp  * w_0_10 
            w00_11 = wp  * w_0_11 
            w01_00 = wp  * w_1_00 
            w01_01 = wp  * w_1_01 
            w01_10 = wp  * w_1_10 
            w01_11 = wp  * w_1_11 
            w10_00 = wp1 * w_0_00 
            w10_01 = wp1 * w_0_01 
            w10_10 = wp1 * w_0_10 
            w10_11 = wp1 * w_0_11 
            w11_00 = wp1 * w_1_00 
            w11_01 = wp1 * w_1_01 
            w11_10 = wp1 * w_1_10 
            w11_11 = wp1 * w_1_11 

            log_uc  = min(log10(max(uch2o * fdif, min_u_h2o)), max_lu_h2o)
            dvar = (log_uc - min_lu_h2o) / dlu_h2o
            iuc = min(max(int(aint(dvar,r8)) + 1, 1), n_u - 1)
            iuc1 = iuc + 1
            wuc = dvar - floor(dvar)
            wuc1 = 1.0 - wuc
!
! asymptotic value of absorptivity as u->infinity
!
            fa = fat(1,ib) + &
                 fat(2,ib) * te1 + &
                 fat(3,ib) * te2 + &
                 fat(4,ib) * te3 + &
                 fat(5,ib) * te4 + &
                 fat(6,ib) * te5

            l_star = &
                 ln_ah2ow(ip , itp , iu , ite , irh ) * w11_11 * wu1 + &
                 ln_ah2ow(ip , itp , iu , ite , irh1) * w11_10 * wu1 + &
                 ln_ah2ow(ip , itp , iu , ite1, irh ) * w11_01 * wu1 + &
                 ln_ah2ow(ip , itp , iu , ite1, irh1) * w11_00 * wu1 + &
                 ln_ah2ow(ip , itp , iu1, ite , irh ) * w11_11 * wu  + &
                 ln_ah2ow(ip , itp , iu1, ite , irh1) * w11_10 * wu  + &
                 ln_ah2ow(ip , itp , iu1, ite1, irh ) * w11_01 * wu  + &
                 ln_ah2ow(ip , itp , iu1, ite1, irh1) * w11_00 * wu  + &
                 ln_ah2ow(ip , itp1, iu , ite , irh ) * w10_11 * wu1 + &
                 ln_ah2ow(ip , itp1, iu , ite , irh1) * w10_10 * wu1 + &
                 ln_ah2ow(ip , itp1, iu , ite1, irh ) * w10_01 * wu1 + &
                 ln_ah2ow(ip , itp1, iu , ite1, irh1) * w10_00 * wu1 + &
                 ln_ah2ow(ip , itp1, iu1, ite , irh ) * w10_11 * wu  + &
                 ln_ah2ow(ip , itp1, iu1, ite , irh1) * w10_10 * wu  + &
                 ln_ah2ow(ip , itp1, iu1, ite1, irh ) * w10_01 * wu  + &
                 ln_ah2ow(ip , itp1, iu1, ite1, irh1) * w10_00 * wu  + &
                 ln_ah2ow(ip1, itp , iu , ite , irh ) * w01_11 * wu1 + &
                 ln_ah2ow(ip1, itp , iu , ite , irh1) * w01_10 * wu1 + &
                 ln_ah2ow(ip1, itp , iu , ite1, irh ) * w01_01 * wu1 + &
                 ln_ah2ow(ip1, itp , iu , ite1, irh1) * w01_00 * wu1 + &
                 ln_ah2ow(ip1, itp , iu1, ite , irh ) * w01_11 * wu  + &
                 ln_ah2ow(ip1, itp , iu1, ite , irh1) * w01_10 * wu  + &
                 ln_ah2ow(ip1, itp , iu1, ite1, irh ) * w01_01 * wu  + &
                 ln_ah2ow(ip1, itp , iu1, ite1, irh1) * w01_00 * wu  + &
                 ln_ah2ow(ip1, itp1, iu , ite , irh ) * w00_11 * wu1 + &
                 ln_ah2ow(ip1, itp1, iu , ite , irh1) * w00_10 * wu1 + &
                 ln_ah2ow(ip1, itp1, iu , ite1, irh ) * w00_01 * wu1 + &
                 ln_ah2ow(ip1, itp1, iu , ite1, irh1) * w00_00 * wu1 + &
                 ln_ah2ow(ip1, itp1, iu1, ite , irh ) * w00_11 * wu  + &
                 ln_ah2ow(ip1, itp1, iu1, ite , irh1) * w00_10 * wu  + &
                 ln_ah2ow(ip1, itp1, iu1, ite1, irh ) * w00_01 * wu  + &
                 ln_ah2ow(ip1, itp1, iu1, ite1, irh1) * w00_00 * wu 

            c_star = &
                 cn_ah2ow(ip , itp , iuc , ite , irh ) * w11_11 * wuc1 + &
                 cn_ah2ow(ip , itp , iuc , ite , irh1) * w11_10 * wuc1 + &
                 cn_ah2ow(ip , itp , iuc , ite1, irh ) * w11_01 * wuc1 + &
                 cn_ah2ow(ip , itp , iuc , ite1, irh1) * w11_00 * wuc1 + &
                 cn_ah2ow(ip , itp , iuc1, ite , irh ) * w11_11 * wuc  + &
                 cn_ah2ow(ip , itp , iuc1, ite , irh1) * w11_10 * wuc  + &
                 cn_ah2ow(ip , itp , iuc1, ite1, irh ) * w11_01 * wuc  + &
                 cn_ah2ow(ip , itp , iuc1, ite1, irh1) * w11_00 * wuc  + &
                 cn_ah2ow(ip , itp1, iuc , ite , irh ) * w10_11 * wuc1 + &
                 cn_ah2ow(ip , itp1, iuc , ite , irh1) * w10_10 * wuc1 + &
                 cn_ah2ow(ip , itp1, iuc , ite1, irh ) * w10_01 * wuc1 + &
                 cn_ah2ow(ip , itp1, iuc , ite1, irh1) * w10_00 * wuc1 + &
                 cn_ah2ow(ip , itp1, iuc1, ite , irh ) * w10_11 * wuc  + &
                 cn_ah2ow(ip , itp1, iuc1, ite , irh1) * w10_10 * wuc  + &
                 cn_ah2ow(ip , itp1, iuc1, ite1, irh ) * w10_01 * wuc  + &
                 cn_ah2ow(ip , itp1, iuc1, ite1, irh1) * w10_00 * wuc  + &
                 cn_ah2ow(ip1, itp , iuc , ite , irh ) * w01_11 * wuc1 + &
                 cn_ah2ow(ip1, itp , iuc , ite , irh1) * w01_10 * wuc1 + &
                 cn_ah2ow(ip1, itp , iuc , ite1, irh ) * w01_01 * wuc1 + &
                 cn_ah2ow(ip1, itp , iuc , ite1, irh1) * w01_00 * wuc1 + &
                 cn_ah2ow(ip1, itp , iuc1, ite , irh ) * w01_11 * wuc  + &
                 cn_ah2ow(ip1, itp , iuc1, ite , irh1) * w01_10 * wuc  + &
                 cn_ah2ow(ip1, itp , iuc1, ite1, irh ) * w01_01 * wuc  + &
                 cn_ah2ow(ip1, itp , iuc1, ite1, irh1) * w01_00 * wuc  + &
                 cn_ah2ow(ip1, itp1, iuc , ite , irh ) * w00_11 * wuc1 + &
                 cn_ah2ow(ip1, itp1, iuc , ite , irh1) * w00_10 * wuc1 + &
                 cn_ah2ow(ip1, itp1, iuc , ite1, irh ) * w00_01 * wuc1 + &
                 cn_ah2ow(ip1, itp1, iuc , ite1, irh1) * w00_00 * wuc1 + &
                 cn_ah2ow(ip1, itp1, iuc1, ite , irh ) * w00_11 * wuc  + &
                 cn_ah2ow(ip1, itp1, iuc1, ite , irh1) * w00_10 * wuc  + &
                 cn_ah2ow(ip1, itp1, iuc1, ite1, irh ) * w00_01 * wuc  + &
                 cn_ah2ow(ip1, itp1, iuc1, ite1, irh1) * w00_00 * wuc 
            abso(i,ib) = min(max(fa * (1.0 - l_star * c_star * &
                                 aer_trn_ttl(i,k1,k2,ib)), &
                             0.0_r8), 1.0_r8) 
!
! invoke linear limit for scaling wrt u below min_u_h2o
!
            if (uvar < min_u_h2o) then
               uscl = uvar / min_u_h2o
               abso(i,ib) = abso(i,ib) * uscl
            endif

         end do
!
! line transmission in 800-1000 and 1000-1200 cm-1 intervals
!
         do i=1,ncol
            term7(i,1) = coefj(1,1) + coefj(2,1)*dty(i)*(1. + c16*dty(i))
            term8(i,1) = coefk(1,1) + coefk(2,1)*dty(i)*(1. + c17*dty(i))
            term7(i,2) = coefj(1,2) + coefj(2,2)*dty(i)*(1. + c26*dty(i))
            term8(i,2) = coefk(1,2) + coefk(2,2)*dty(i)*(1. + c27*dty(i))
         end do
!
! 500 -  800 cm-1   h2o rotation band overlap with co2
!
         do i=1,ncol
            k21    = term7(i,1) + term8(i,1)/ &
               (1. + (c30 + c31*(dty(i)-10.)*(dty(i)-10.))*sqrtu(i))
            k22    = term7(i,2) + term8(i,2)/ &
               (1. + (c28 + c29*(dty(i)-10.))*sqrtu(i))
            tr1    = exp(-(k21*(sqrtu(i) + fc1*fwku(i))))
            tr2    = exp(-(k22*(sqrtu(i) + fc1*fwku(i))))
            tr1=tr1*aer_trn_ttl(i,k1,k2,idx_lw_0650_0800) 
!                                          ! h2o line+straer trn 650--800 cm-1
            tr2=tr2*aer_trn_ttl(i,k1,k2,idx_lw_0500_0650)
!                                          ! h2o line+straer trn 500--650 cm-1
            tr5    = exp(-((coefh(1,3) + coefh(2,3)*dtx(i))*uc1(i)))
            tr6    = exp(-((coefh(1,4) + coefh(2,4)*dtx(i))*uc1(i)))
            tr9(i)   = tr1*tr5
            tr10(i)  = tr2*tr6
            th2o(i) = tr10(i)
            trab2(i) = 0.65*tr9(i) + 0.35*tr10(i)
         end do
         if (k2 < k1) then
            do i=1,ncol
               to3h2o(i) = h2otr(i,k1)/h2otr(i,k2)
            end do
         else
            do i=1,ncol
               to3h2o(i) = h2otr(i,k2)/h2otr(i,k1)
            end do
         end if
!
! abso(i,3)   o3  9.6 micrometer band (nu3 and nu1 bands)
!
         do i=1,ncol
            dpnm(i)  = pnm(i,k1) - pnm(i,k2)
            to3co2(i) = (pnm(i,k1)*co2t(i,k1) - pnm(i,k2)*co2t(i,k2))/dpnm(i)
            te       = (to3co2(i)*r293)**.7
            dplos    = plos(i,k1) - plos(i,k2)
            dplol    = plol(i,k1) - plol(i,k2)
            u1       = 18.29*abs(dplos)/te
            u2       = .5649*abs(dplos)/te
            rphat    = dplol/dplos
            tlocal   = tint(i,k2)
            tcrfac   = sqrt(tlocal*r250)*te
            beta     = r3205*(rphat + dpfo3*tcrfac)
            realnu   = te/beta
            tmp1     = u1/sqrt(4. + u1*(1. + realnu))
            tmp2     = u2/sqrt(4. + u2*(1. + realnu))
            o3bndi    = 74.*te*log(1. + tmp1 + tmp2)
            abso(i,3) = o3bndi*to3h2o(i)*dbvtit(i,k2)
            to3(i)   = 1.0/(1. + 0.1*tmp1 + 0.1*tmp2)
         end do
!
! abso(i,4)      co2 15  micrometer band system
!
         do i=1,ncol
            sqwp      = sqrt(abs(plco2(i,k1) - plco2(i,k2)))
            et        = exp(-480./to3co2(i))
            sqti(i)   = sqrt(to3co2(i))
            rsqti     = 1./sqti(i)
            et2       = et*et
            et4       = et2*et2
            omet      = 1. - 1.5*et2
            f1co2     = 899.70*omet*(1. + 1.94774*et + 4.73486*et2)*rsqti
            f1sqwp(i) = f1co2*sqwp
            t1co2(i)  = 1./(1. + (245.18*omet*sqwp*rsqti))
            oneme     = 1. - et2
            alphat    = oneme**3*rsqti
            pi        = abs(dpnm(i))
            wco2      =  2.5221*co2vmr*pi*rga
            u7(i)     =  4.9411e4*alphat*et2*wco2
            u8        =  3.9744e4*alphat*et4*wco2
            u9        =  1.0447e5*alphat*et4*et2*wco2
            u13       = 2.8388e3*alphat*et4*wco2
            tpath     = to3co2(i)
            tlocal    = tint(i,k2)
            tcrfac    = sqrt(tlocal*r250*tpath*r300)
            posqt     = ((pnm(i,k2) + pnm(i,k1))*r2sslp + dpfco2*tcrfac)*rsqti
            rbeta7(i) = 1./(5.3228*posqt)
            rbeta8    = 1./(10.6576*posqt)
            rbeta9    = rbeta7(i)
            rbeta13   = rbeta9
            f2co2(i)  = (u7(i)/sqrt(4. + u7(i)*(1. + rbeta7(i)))) + &
               (u8   /sqrt(4. + u8*(1. + rbeta8))) + &
               (u9   /sqrt(4. + u9*(1. + rbeta9)))
            f3co2(i)  = u13/sqrt(4. + u13*(1. + rbeta13))
         end do
         if (k2 >= k1) then
            do i=1,ncol
               sqti(i) = sqrt(tlayr(i,k2))
            end do
         end if
!
         do i=1,ncol
            tmp1      = log(1. + f1sqwp(i))
            tmp2      = log(1. + f2co2(i))
            tmp3      = log(1. + f3co2(i))
            absbnd    = (tmp1 + 2.*t1co2(i)*tmp2 + 2.*tmp3)*sqti(i)
            abso(i,4) = trab2(i)*co2em(i,k2)*absbnd
            tco2(i)   = 1./(1.0+10.0*(u7(i)/sqrt(4. + u7(i)*(1. + rbeta7(i)))))
         end do
!
! calculate absorptivity due to trace gases, abstrc
!
         call trcab( lchnk   ,ncol    ,pcols, pverp,                   &
            k1      ,k2      ,ucfc11  ,ucfc12  ,un2o0   , &
            un2o1   ,uch4    ,uco211  ,uco212  ,uco213  , &
            uco221  ,uco222  ,uco223  ,bn2o0   ,bn2o1   , &
            bch4    ,to3co2  ,pnm     ,dw      ,pnew    , &
            s2c     ,uptype  ,u       ,abplnk1 ,tco2    , &
            th2o    ,to3     ,abstrc  , &
            aer_trn_ttl)
!
! sum total absorptivity
!
         do i=1,ncol
            abstot(i,k1,k2) = abso(i,1) + abso(i,2) + &
               abso(i,3) + abso(i,4) + abstrc(i)
         end do
      end do ! do k2 = 
   end do ! do k1 = 
!
! adjacent layer absorptivity:
!
! abso(i,1)     0 -  800 cm-1   h2o rotation band
! abso(i,1)  1200 - 2200 cm-1   h2o vibration-rotation band
! abso(i,2)   800 - 1200 cm-1   h2o window
!
! separation between rotation and vibration-rotation dropped, so
!                only 2 slots needed for h2o absorptivity
!
! 500-800 cm^-1 h2o continuum/line overlap already included
!                in abso(i,1).  this used to be in abso(i,4)
!
! abso(i,3)   o3  9.6 micrometer band (nu3 and nu1 bands)
! abso(i,4)   co2 15  micrometer band system
!
! nearest layer level loop
!
   do k2=pver,ntoplw,-1
      do i=1,ncol
         tbar(i,1)   = 0.5*(tint(i,k2+1) + tlayr(i,k2+1))
         emm(i,1)    = 0.5*(co2em(i,k2+1) + co2eml(i,k2))
         tbar(i,2)   = 0.5*(tlayr(i,k2+1) + tint(i,k2))
         emm(i,2)    = 0.5*(co2em(i,k2) + co2eml(i,k2))
         tbar(i,3)   = 0.5*(tbar(i,2) + tbar(i,1))
         emm(i,3)    = emm(i,1)
         tbar(i,4)   = tbar(i,3)
         emm(i,4)    = emm(i,2)
         o3emm(i,1)  = 0.5*(dbvtit(i,k2+1) + dbvtly(i,k2))
         o3emm(i,2)  = 0.5*(dbvtit(i,k2) + dbvtly(i,k2))
         o3emm(i,3)  = o3emm(i,1)
         o3emm(i,4)  = o3emm(i,2)
         temh2o(i,1) = tbar(i,1)
         temh2o(i,2) = tbar(i,2)
         temh2o(i,3) = tbar(i,1)
         temh2o(i,4) = tbar(i,2)
         dpnm(i)     = pnm(i,k2+1) - pnm(i,k2)
      end do
!
!  weighted planck functions for trace gases
!
      do wvl = 1,14
         do i = 1,ncol
            bplnk(wvl,i,1) = 0.5*(abplnk1(wvl,i,k2+1) + abplnk2(wvl,i,k2))
            bplnk(wvl,i,2) = 0.5*(abplnk1(wvl,i,k2) + abplnk2(wvl,i,k2))
            bplnk(wvl,i,3) = bplnk(wvl,i,1)
            bplnk(wvl,i,4) = bplnk(wvl,i,2)
         end do
      end do
      
      do i=1,ncol
         rdpnmsq    = 1./(pnmsq(i,k2+1) - pnmsq(i,k2))
         rdpnm      = 1./dpnm(i)
         p1         = .5*(pbr(i,k2) + pnm(i,k2+1))
         p2         = .5*(pbr(i,k2) + pnm(i,k2  ))
         uinpl(i,1) =  (pnmsq(i,k2+1) - p1**2)*rdpnmsq
         uinpl(i,2) = -(pnmsq(i,k2  ) - p2**2)*rdpnmsq
         uinpl(i,3) = -(pnmsq(i,k2  ) - p1**2)*rdpnmsq
         uinpl(i,4) =  (pnmsq(i,k2+1) - p2**2)*rdpnmsq
         winpl(i,1) = (.5*( pnm(i,k2+1) - pbr(i,k2)))*rdpnm
         winpl(i,2) = (.5*(-pnm(i,k2  ) + pbr(i,k2)))*rdpnm
         winpl(i,3) = (.5*( pnm(i,k2+1) + pbr(i,k2)) - pnm(i,k2  ))*rdpnm
         winpl(i,4) = (.5*(-pnm(i,k2  ) - pbr(i,k2)) + pnm(i,k2+1))*rdpnm
         tmp1       = 1./(piln(i,k2+1) - piln(i,k2))
         tmp2       = piln(i,k2+1) - pmln(i,k2)
         tmp3       = piln(i,k2  ) - pmln(i,k2)
         zinpl(i,1) = (.5*tmp2          )*tmp1
         zinpl(i,2) = (        - .5*tmp3)*tmp1
         zinpl(i,3) = (.5*tmp2 -    tmp3)*tmp1
         zinpl(i,4) = (   tmp2 - .5*tmp3)*tmp1
         pinpl(i,1) = 0.5*(p1 + pnm(i,k2+1))
         pinpl(i,2) = 0.5*(p2 + pnm(i,k2  ))
         pinpl(i,3) = 0.5*(p1 + pnm(i,k2  ))
         pinpl(i,4) = 0.5*(p2 + pnm(i,k2+1))
         if(strat_volcanic) then
           aer_pth_ngh(i) = abs(aer_mpp(i,k2)-aer_mpp(i,k2+1))
         endif
      end do
      do kn=1,4
         do i=1,ncol
            u(i)     = uinpl(i,kn)*abs(plh2o(i,k2) - plh2o(i,k2+1))
            sqrtu(i) = sqrt(u(i))
            dw(i)    = abs(w(i,k2) - w(i,k2+1))
            pnew(i)  = u(i)/(winpl(i,kn)*dw(i))
            pnew_mks  = pnew(i) * sslp_mks
            t_p = min(max(tbar(i,kn), min_tp_h2o), max_tp_h2o)
            iest = floor(t_p) - min_tp_h2o
            esx = estblh2o(iest) + (estblh2o(iest+1)-estblh2o(iest)) * &
                 (t_p - min_tp_h2o - iest)
            qsx = epsilo * esx / (pnew_mks - omeps * esx)
            q_path = dw(i) / abs(dpnm(i)) / rga
            
            ds2c     = abs(s2c(i,k2) - s2c(i,k2+1))
            uc1(i)   = uinpl(i,kn)*ds2c
            pch2o    = uc1(i)
            uc1(i)   = (uc1(i) + 1.7e-3*u(i))*(1. +  2.*uc1(i))/(1. + 15.*uc1(i))
            dtx(i)      = temh2o(i,kn) - 250.
            dty(i)      = tbar(i,kn) - 250.
            
            fwk(i)    = fwcoef + fwc1/(1. + fwc2*u(i))
            fwku(i)   = fwk(i)*u(i)

            if(strat_volcanic) then
              aer_pth_dlt=uinpl(i,kn)*aer_pth_ngh(i)
  
              do bnd_idx=1,bnd_nbr_lw
                 odap_aer_ttl=abs_cff_mss_aer(bnd_idx) * aer_pth_dlt 
                 aer_trn_ngh(i,bnd_idx)=exp(-fdif * odap_aer_ttl)
              end do
            else
              aer_trn_ngh(i,:) = 1.0
            endif

!
! define variables for c/h/e (now c/lt/e) fit
!
! abso(i,1)     0 -  800 cm-1   h2o rotation band
! abso(i,1)  1200 - 2200 cm-1   h2o vibration-rotation band
! abso(i,2)   800 - 1200 cm-1   h2o window
!
! separation between rotation and vibration-rotation dropped, so
!                only 2 slots needed for h2o absorptivity
!
! notation:
! u   = integral (p/p_0 dw)  
! p   = atmospheric pressure
! p_0 = reference atmospheric pressure
! w   = precipitable water path
! t_e = emission temperature
! t_p = path temperature
! rh  = path relative humidity
!
!
! terms for asymptotic value of emissivity
!
            te1  = temh2o(i,kn)
            te2  = te1 * te1
            te3  = te2 * te1
            te4  = te3 * te1
            te5  = te4 * te1

!
! indices for lines and continuum tables 
! note: because we are dealing with the nearest layer,
!       the hulst-curtis-godson corrections
!       for inhomogeneous paths are not applied.
!
            uvar = u(i)*fdif
            log_u  = min(log10(max(uvar, min_u_h2o)), max_lu_h2o)
            dvar = (log_u - min_lu_h2o) / dlu_h2o
            iu = min(max(int(aint(dvar,r8)) + 1, 1), n_u - 1)
            iu1 = iu + 1
            wu = dvar - floor(dvar)
            wu1 = 1.0 - wu
            
            log_p  = min(log10(max(pnew(i), min_p_h2o)), max_lp_h2o)
            dvar = (log_p - min_lp_h2o) / dlp_h2o
            ip = min(max(int(aint(dvar,r8)) + 1, 1), n_p - 1)
            ip1 = ip + 1
            wp = dvar - floor(dvar)
            wp1 = 1.0 - wp
            
            dvar = (t_p - min_tp_h2o) / dtp_h2o
            itp = min(max(int(aint(dvar,r8)) + 1, 1), n_tp - 1)
            itp1 = itp + 1
            wtp = dvar - floor(dvar)
            wtp1 = 1.0 - wtp
            
            t_e = min(max(temh2o(i,kn)-t_p,min_te_h2o),max_te_h2o)
            dvar = (t_e - min_te_h2o) / dte_h2o
            ite = min(max(int(aint(dvar,r8)) + 1, 1), n_te - 1)
            ite1 = ite + 1
            wte = dvar - floor(dvar)
            wte1 = 1.0 - wte
            
            rh_path = min(max(q_path / qsx, min_rh_h2o), max_rh_h2o)
            dvar = (rh_path - min_rh_h2o) / drh_h2o
            irh = min(max(int(aint(dvar,r8)) + 1, 1), n_rh - 1)
            irh1 = irh + 1
            wrh = dvar - floor(dvar)
            wrh1 = 1.0 - wrh
            
            w_0_0_ = wtp  * wte
            w_0_1_ = wtp  * wte1
            w_1_0_ = wtp1 * wte 
            w_1_1_ = wtp1 * wte1
            
            w_0_00 = w_0_0_ * wrh
            w_0_01 = w_0_0_ * wrh1
            w_0_10 = w_0_1_ * wrh
            w_0_11 = w_0_1_ * wrh1
            w_1_00 = w_1_0_ * wrh
            w_1_01 = w_1_0_ * wrh1
            w_1_10 = w_1_1_ * wrh
            w_1_11 = w_1_1_ * wrh1
            
            w00_00 = wp  * w_0_00 
            w00_01 = wp  * w_0_01 
            w00_10 = wp  * w_0_10 
            w00_11 = wp  * w_0_11 
            w01_00 = wp  * w_1_00 
            w01_01 = wp  * w_1_01 
            w01_10 = wp  * w_1_10 
            w01_11 = wp  * w_1_11 
            w10_00 = wp1 * w_0_00 
            w10_01 = wp1 * w_0_01 
            w10_10 = wp1 * w_0_10 
            w10_11 = wp1 * w_0_11 
            w11_00 = wp1 * w_1_00 
            w11_01 = wp1 * w_1_01 
            w11_10 = wp1 * w_1_10 
            w11_11 = wp1 * w_1_11 

!
! non-window absorptivity
!
            ib = 1
            
            fa = fat(1,ib) + &
                 fat(2,ib) * te1 + &
                 fat(3,ib) * te2 + &
                 fat(4,ib) * te3 + &
                 fat(5,ib) * te4 + &
                 fat(6,ib) * te5
            
            a_star = &
                 ah2onw(ip , itp , iu , ite , irh ) * w11_11 * wu1 + &
                 ah2onw(ip , itp , iu , ite , irh1) * w11_10 * wu1 + &
                 ah2onw(ip , itp , iu , ite1, irh ) * w11_01 * wu1 + &
                 ah2onw(ip , itp , iu , ite1, irh1) * w11_00 * wu1 + &
                 ah2onw(ip , itp , iu1, ite , irh ) * w11_11 * wu  + &
                 ah2onw(ip , itp , iu1, ite , irh1) * w11_10 * wu  + &
                 ah2onw(ip , itp , iu1, ite1, irh ) * w11_01 * wu  + &
                 ah2onw(ip , itp , iu1, ite1, irh1) * w11_00 * wu  + &
                 ah2onw(ip , itp1, iu , ite , irh ) * w10_11 * wu1 + &
                 ah2onw(ip , itp1, iu , ite , irh1) * w10_10 * wu1 + &
                 ah2onw(ip , itp1, iu , ite1, irh ) * w10_01 * wu1 + &
                 ah2onw(ip , itp1, iu , ite1, irh1) * w10_00 * wu1 + &
                 ah2onw(ip , itp1, iu1, ite , irh ) * w10_11 * wu  + &
                 ah2onw(ip , itp1, iu1, ite , irh1) * w10_10 * wu  + &
                 ah2onw(ip , itp1, iu1, ite1, irh ) * w10_01 * wu  + &
                 ah2onw(ip , itp1, iu1, ite1, irh1) * w10_00 * wu  + &
                 ah2onw(ip1, itp , iu , ite , irh ) * w01_11 * wu1 + &
                 ah2onw(ip1, itp , iu , ite , irh1) * w01_10 * wu1 + &
                 ah2onw(ip1, itp , iu , ite1, irh ) * w01_01 * wu1 + &
                 ah2onw(ip1, itp , iu , ite1, irh1) * w01_00 * wu1 + &
                 ah2onw(ip1, itp , iu1, ite , irh ) * w01_11 * wu  + &
                 ah2onw(ip1, itp , iu1, ite , irh1) * w01_10 * wu  + &
                 ah2onw(ip1, itp , iu1, ite1, irh ) * w01_01 * wu  + &
                 ah2onw(ip1, itp , iu1, ite1, irh1) * w01_00 * wu  + &
                 ah2onw(ip1, itp1, iu , ite , irh ) * w00_11 * wu1 + &
                 ah2onw(ip1, itp1, iu , ite , irh1) * w00_10 * wu1 + &
                 ah2onw(ip1, itp1, iu , ite1, irh ) * w00_01 * wu1 + &
                 ah2onw(ip1, itp1, iu , ite1, irh1) * w00_00 * wu1 + &
                 ah2onw(ip1, itp1, iu1, ite , irh ) * w00_11 * wu  + &
                 ah2onw(ip1, itp1, iu1, ite , irh1) * w00_10 * wu  + &
                 ah2onw(ip1, itp1, iu1, ite1, irh ) * w00_01 * wu  + &
                 ah2onw(ip1, itp1, iu1, ite1, irh1) * w00_00 * wu
            
            abso(i,ib) = min(max(fa * (1.0 - (1.0 - a_star) * &
                                 aer_trn_ngh(i,ib)), &
                             0.0_r8), 1.0_r8)

!
! invoke linear limit for scaling wrt u below min_u_h2o
!
            if (uvar < min_u_h2o) then
               uscl = uvar / min_u_h2o
               abso(i,ib) = abso(i,ib) * uscl
            endif
            
!
! window absorptivity
!
            ib = 2
            
            fa = fat(1,ib) + &
                 fat(2,ib) * te1 + &
                 fat(3,ib) * te2 + &
                 fat(4,ib) * te3 + &
                 fat(5,ib) * te4 + &
                 fat(6,ib) * te5
            
            a_star = &
                 ah2ow(ip , itp , iu , ite , irh ) * w11_11 * wu1 + &
                 ah2ow(ip , itp , iu , ite , irh1) * w11_10 * wu1 + &
                 ah2ow(ip , itp , iu , ite1, irh ) * w11_01 * wu1 + &
                 ah2ow(ip , itp , iu , ite1, irh1) * w11_00 * wu1 + &
                 ah2ow(ip , itp , iu1, ite , irh ) * w11_11 * wu  + &
                 ah2ow(ip , itp , iu1, ite , irh1) * w11_10 * wu  + &
                 ah2ow(ip , itp , iu1, ite1, irh ) * w11_01 * wu  + &
                 ah2ow(ip , itp , iu1, ite1, irh1) * w11_00 * wu  + &
                 ah2ow(ip , itp1, iu , ite , irh ) * w10_11 * wu1 + &
                 ah2ow(ip , itp1, iu , ite , irh1) * w10_10 * wu1 + &
                 ah2ow(ip , itp1, iu , ite1, irh ) * w10_01 * wu1 + &
                 ah2ow(ip , itp1, iu , ite1, irh1) * w10_00 * wu1 + &
                 ah2ow(ip , itp1, iu1, ite , irh ) * w10_11 * wu  + &
                 ah2ow(ip , itp1, iu1, ite , irh1) * w10_10 * wu  + &
                 ah2ow(ip , itp1, iu1, ite1, irh ) * w10_01 * wu  + &
                 ah2ow(ip , itp1, iu1, ite1, irh1) * w10_00 * wu  + &
                 ah2ow(ip1, itp , iu , ite , irh ) * w01_11 * wu1 + &
                 ah2ow(ip1, itp , iu , ite , irh1) * w01_10 * wu1 + &
                 ah2ow(ip1, itp , iu , ite1, irh ) * w01_01 * wu1 + &
                 ah2ow(ip1, itp , iu , ite1, irh1) * w01_00 * wu1 + &
                 ah2ow(ip1, itp , iu1, ite , irh ) * w01_11 * wu  + &
                 ah2ow(ip1, itp , iu1, ite , irh1) * w01_10 * wu  + &
                 ah2ow(ip1, itp , iu1, ite1, irh ) * w01_01 * wu  + &
                 ah2ow(ip1, itp , iu1, ite1, irh1) * w01_00 * wu  + &
                 ah2ow(ip1, itp1, iu , ite , irh ) * w00_11 * wu1 + &
                 ah2ow(ip1, itp1, iu , ite , irh1) * w00_10 * wu1 + &
                 ah2ow(ip1, itp1, iu , ite1, irh ) * w00_01 * wu1 + &
                 ah2ow(ip1, itp1, iu , ite1, irh1) * w00_00 * wu1 + &
                 ah2ow(ip1, itp1, iu1, ite , irh ) * w00_11 * wu  + &
                 ah2ow(ip1, itp1, iu1, ite , irh1) * w00_10 * wu  + &
                 ah2ow(ip1, itp1, iu1, ite1, irh ) * w00_01 * wu  + &
                 ah2ow(ip1, itp1, iu1, ite1, irh1) * w00_00 * wu
            
            abso(i,ib) = min(max(fa * (1.0 - (1.0 - a_star) * &
                                 aer_trn_ngh(i,ib)), &
                             0.0_r8), 1.0_r8)

!
! invoke linear limit for scaling wrt u below min_u_h2o
!
            if (uvar < min_u_h2o) then
               uscl = uvar / min_u_h2o
               abso(i,ib) = abso(i,ib) * uscl
            endif
            
         end do
!
! line transmission in 800-1000 and 1000-1200 cm-1 intervals
!
         do i=1,ncol
            term7(i,1) = coefj(1,1) + coefj(2,1)*dty(i)*(1. + c16*dty(i))
            term8(i,1) = coefk(1,1) + coefk(2,1)*dty(i)*(1. + c17*dty(i))
            term7(i,2) = coefj(1,2) + coefj(2,2)*dty(i)*(1. + c26*dty(i))
            term8(i,2) = coefk(1,2) + coefk(2,2)*dty(i)*(1. + c27*dty(i))
         end do
!
! 500 -  800 cm-1   h2o rotation band overlap with co2
!
         do i=1,ncol
            dtym10     = dty(i) - 10.
            denom      = 1. + (c30 + c31*dtym10*dtym10)*sqrtu(i)
            k21        = term7(i,1) + term8(i,1)/denom
            denom      = 1. + (c28 + c29*dtym10       )*sqrtu(i)
            k22        = term7(i,2) + term8(i,2)/denom
            tr1     = exp(-(k21*(sqrtu(i) + fc1*fwku(i))))
            tr2     = exp(-(k22*(sqrtu(i) + fc1*fwku(i))))
            tr1=tr1*aer_trn_ngh(i,idx_lw_0650_0800) 
!                                         ! h2o line+straer trn 650--800 cm-1
            tr2=tr2*aer_trn_ngh(i,idx_lw_0500_0650) 
!                                         ! h2o line+straer trn 500--650 cm-1
            tr5     = exp(-((coefh(1,3) + coefh(2,3)*dtx(i))*uc1(i)))
            tr6     = exp(-((coefh(1,4) + coefh(2,4)*dtx(i))*uc1(i)))
            tr9(i)  = tr1*tr5
            tr10(i) = tr2*tr6
            trab2(i)= 0.65*tr9(i) + 0.35*tr10(i)
            th2o(i) = tr10(i)
         end do
!
! abso(i,3)  o3  9.6 micrometer (nu3 and nu1 bands)
!
         do i=1,ncol
            te        = (tbar(i,kn)*r293)**.7
            dplos     = abs(plos(i,k2+1) - plos(i,k2))
            u1        = zinpl(i,kn)*18.29*dplos/te
            u2        = zinpl(i,kn)*.5649*dplos/te
            tlocal    = tbar(i,kn)
            tcrfac    = sqrt(tlocal*r250)*te
            beta      = r3205*(pinpl(i,kn)*rsslp + dpfo3*tcrfac)
            realnu    = te/beta
            tmp1      = u1/sqrt(4. + u1*(1. + realnu))
            tmp2      = u2/sqrt(4. + u2*(1. + realnu))
            o3bndi    = 74.*te*log(1. + tmp1 + tmp2)
            abso(i,3) = o3bndi*o3emm(i,kn)*(h2otr(i,k2+1)/h2otr(i,k2))
            to3(i)    = 1.0/(1. + 0.1*tmp1 + 0.1*tmp2)
         end do
!
! abso(i,4)   co2 15  micrometer band system
!
         do i=1,ncol
            dplco2   = plco2(i,k2+1) - plco2(i,k2)
            sqwp     = sqrt(uinpl(i,kn)*dplco2)
            et       = exp(-480./tbar(i,kn))
            sqti(i)  = sqrt(tbar(i,kn))
            rsqti    = 1./sqti(i)
            et2      = et*et
            et4      = et2*et2
            omet     = (1. - 1.5*et2)
            f1co2    = 899.70*omet*(1. + 1.94774*et + 4.73486*et2)*rsqti
            f1sqwp(i)= f1co2*sqwp
            t1co2(i) = 1./(1. + (245.18*omet*sqwp*rsqti))
            oneme    = 1. - et2
            alphat   = oneme**3*rsqti
            pi       = abs(dpnm(i))*winpl(i,kn)
            wco2     = 2.5221*co2vmr*pi*rga
            u7(i)    = 4.9411e4*alphat*et2*wco2
            u8       = 3.9744e4*alphat*et4*wco2
            u9       = 1.0447e5*alphat*et4*et2*wco2
            u13      = 2.8388e3*alphat*et4*wco2
            tpath    = tbar(i,kn)
            tlocal   = tbar(i,kn)
            tcrfac   = sqrt((tlocal*r250)*(tpath*r300))
            posqt    = (pinpl(i,kn)*rsslp + dpfco2*tcrfac)*rsqti
            rbeta7(i)= 1./(5.3228*posqt)
            rbeta8   = 1./(10.6576*posqt)
            rbeta9   = rbeta7(i)
            rbeta13  = rbeta9
            f2co2(i) = u7(i)/sqrt(4. + u7(i)*(1. + rbeta7(i))) + &
                 u8   /sqrt(4. + u8*(1. + rbeta8)) + &
                 u9   /sqrt(4. + u9*(1. + rbeta9))
            f3co2(i) = u13/sqrt(4. + u13*(1. + rbeta13))
            tmp1     = log(1. + f1sqwp(i))
            tmp2     = log(1. + f2co2(i))
            tmp3     = log(1. + f3co2(i))
            absbnd   = (tmp1 + 2.*t1co2(i)*tmp2 + 2.*tmp3)*sqti(i)
            abso(i,4)= trab2(i)*emm(i,kn)*absbnd
            tco2(i)  = 1.0/(1.0+ 10.0*u7(i)/sqrt(4. + u7(i)*(1. + rbeta7(i))))
         end do ! do i =
!
! calculate trace gas absorptivity for nearest layer, abstrc
!
         call trcabn(lchnk   ,ncol    ,pcols, pverp,                   &
              k2      ,kn      ,ucfc11  ,ucfc12  ,un2o0   , &
              un2o1   ,uch4    ,uco211  ,uco212  ,uco213  , &
              uco221  ,uco222  ,uco223  ,tbar    ,bplnk   , &
              winpl   ,pinpl   ,tco2    ,th2o    ,to3     , &
              uptype  ,dw      ,s2c     ,u       ,pnew    , &
              abstrc  ,uinpl   , &
              aer_trn_ngh)
!
! total next layer absorptivity:
!
         do i=1,ncol
            absnxt(i,k2,kn) = abso(i,1) + abso(i,2) + &
                 abso(i,3) + abso(i,4) + abstrc(i)
         end do
      end do ! do kn =
   end do ! do k2 =

   return
end subroutine radabs



subroutine radems(lchnk   ,ncol    ,pcols, pver, pverp,         &
                  s2c     ,tcg     ,w       ,tplnke  ,plh2o   , &
                  pnm     ,plco2   ,tint    ,tint4   ,tlayr   , &
                  tlayr4  ,plol    ,plos    ,ucfc11  ,ucfc12  , &
                  un2o0   ,un2o1   ,uch4    ,uco211 ,uco212   , &
                  uco213  ,uco221  ,uco222  ,uco223  ,uptype  , &
                  bn2o0   ,bn2o1   ,bch4    ,co2em   ,co2eml  , &
                  co2t    ,h2otr   ,abplnk1 ,abplnk2 ,emstot  , &
                  plh2ob  ,wb      , &
                  aer_trn_ttl)
!----------------------------------------------------------------------- 
! 
! purpose: 
! compute emissivity for h2o, co2, o3, ch4, n2o, cfc11 and cfc12
! 
! method: 
! h2o  ....  uses nonisothermal emissivity method for water vapor from
!            ramanathan, v. and  p.downey, 1986: a nonisothermal
!            emissivity and absorptivity formulation for water vapor
!            jouranl of geophysical research, vol. 91., d8, pp 8649-8666
!
!            implementation updated by collins,hackney, and edwards 2001
!               using line-by-line calculations based upon hitran 1996 and
!               ckd 2.1 for absorptivity and emissivity
!
!            implementation updated by collins, lee-taylor, and edwards (2003)
!               using line-by-line calculations based upon hitran 2000 and
!               ckd 2.4 for absorptivity and emissivity
!
! co2  ....  uses absorptance parameterization of the 15 micro-meter
!            (500 - 800 cm-1) band system of carbon dioxide, from
!            kiehl, j.t. and b.p.briegleb, 1991: a new parameterization
!            of the absorptance due to the 15 micro-meter band system
!            of carbon dioxide jouranl of geophysical research,
!            vol. 96., d5, pp 9013-9019. also includes the effects
!            of the 9.4 and 10.4 micron bands of co2.
!
! o3   ....  uses absorptance parameterization of the 9.6 micro-meter
!            band system of ozone, from ramanathan, v. and r. dickinson,
!            1979: the role of stratospheric ozone in the zonal and
!            seasonal radiative energy balance of the earth-troposphere
!            system. journal of the atmospheric sciences, vol. 36,
!            pp 1084-1104
!
! ch4  ....  uses a broad band model for the 7.7 micron band of methane.
!
! n20  ....  uses a broad band model for the 7.8, 8.6 and 17.0 micron
!            bands of nitrous oxide
!
! cfc11 ...  uses a quasi-linear model for the 9.2, 10.7, 11.8 and 12.5
!            micron bands of cfc11
!
! cfc12 ...  uses a quasi-linear model for the 8.6, 9.1, 10.8 and 11.2
!            micron bands of cfc12
!
!
! computes individual emissivities, accounting for band overlap, and
! sums to obtain the total.
!
! author: w. collins (h2o emissivity) and j. kiehl
! 
!-----------------------------------------------------------------------
!------------------------------arguments--------------------------------
!
! input arguments
!
   integer, intent(in) :: lchnk                    ! chunk identifier
   integer, intent(in) :: ncol                     ! number of atmospheric columns
   integer, intent(in) :: pcols, pver, pverp

   real(r8), intent(in) :: s2c(pcols,pverp)        ! h2o continuum path length
   real(r8), intent(in) :: tcg(pcols,pverp)        ! h2o-mass-wgted temp. (curtis-godson approx.)
   real(r8), intent(in) :: w(pcols,pverp)          ! h2o path length
   real(r8), intent(in) :: tplnke(pcols)           ! layer planck temperature
   real(r8), intent(in) :: plh2o(pcols,pverp)      ! h2o prs wghted path length
   real(r8), intent(in) :: pnm(pcols,pverp)        ! model interface pressure
   real(r8), intent(in) :: plco2(pcols,pverp)      ! prs wghted path of co2
   real(r8), intent(in) :: tint(pcols,pverp)       ! model interface temperatures
   real(r8), intent(in) :: tint4(pcols,pverp)      ! tint to the 4th power
   real(r8), intent(in) :: tlayr(pcols,pverp)      ! k-1 model layer temperature
   real(r8), intent(in) :: tlayr4(pcols,pverp)     ! tlayr to the 4th power
   real(r8), intent(in) :: plol(pcols,pverp)       ! pressure wghtd ozone path
   real(r8), intent(in) :: plos(pcols,pverp)       ! ozone path
   real(r8), intent(in) :: plh2ob(nbands,pcols,pverp) ! pressure weighted h2o path with 
                                                      !    hulst-curtis-godson temp. factor 
                                                      !    for h2o bands 
   real(r8), intent(in) :: wb(nbands,pcols,pverp)     ! h2o path length with 
                                                      !    hulst-curtis-godson temp. factor 
                                                      !    for h2o bands 

   real(r8), intent(in) :: aer_trn_ttl(pcols,pverp,pverp,bnd_nbr_lw) 
!                               ! [fraction] total strat. aerosol
!                               ! transmission between interfaces k1 and k2  

!
! trace gas variables
!
   real(r8), intent(in) :: ucfc11(pcols,pverp)     ! cfc11 path length
   real(r8), intent(in) :: ucfc12(pcols,pverp)     ! cfc12 path length
   real(r8), intent(in) :: un2o0(pcols,pverp)      ! n2o path length
   real(r8), intent(in) :: un2o1(pcols,pverp)      ! n2o path length (hot band)
   real(r8), intent(in) :: uch4(pcols,pverp)       ! ch4 path length
   real(r8), intent(in) :: uco211(pcols,pverp)     ! co2 9.4 micron band path length
   real(r8), intent(in) :: uco212(pcols,pverp)     ! co2 9.4 micron band path length
   real(r8), intent(in) :: uco213(pcols,pverp)     ! co2 9.4 micron band path length
   real(r8), intent(in) :: uco221(pcols,pverp)     ! co2 10.4 micron band path length
   real(r8), intent(in) :: uco222(pcols,pverp)     ! co2 10.4 micron band path length
   real(r8), intent(in) :: uco223(pcols,pverp)     ! co2 10.4 micron band path length
   real(r8), intent(in) :: bn2o0(pcols,pverp)      ! pressure factor for n2o
   real(r8), intent(in) :: bn2o1(pcols,pverp)      ! pressure factor for n2o
   real(r8), intent(in) :: bch4(pcols,pverp)       ! pressure factor for ch4
   real(r8), intent(in) :: uptype(pcols,pverp)     ! p-type continuum path length
!
! output arguments
!
   real(r8), intent(out) :: emstot(pcols,pverp)     ! total emissivity
   real(r8), intent(out) :: co2em(pcols,pverp)      ! layer co2 normalzd plnck funct drvtv
   real(r8), intent(out) :: co2eml(pcols,pver)      ! intrfc co2 normalzd plnck func drvtv
   real(r8), intent(out) :: co2t(pcols,pverp)       ! tmp and prs weighted path length
   real(r8), intent(out) :: h2otr(pcols,pverp)      ! h2o transmission over o3 band
   real(r8), intent(out) :: abplnk1(14,pcols,pverp) ! non-nearest layer plack factor
   real(r8), intent(out) :: abplnk2(14,pcols,pverp) ! nearest layer factor

!
!---------------------------local variables-----------------------------
!
   integer i                    ! longitude index
   integer k                    ! level index]
   integer k1                   ! level index
!
! local variables for h2o:
!
   real(r8) h2oems(pcols,pverp)     ! h2o emissivity
   real(r8) tpathe                  ! used to compute h2o emissivity
   real(r8) dtx(pcols)              ! planck temperature minus 250 k
   real(r8) dty(pcols)              ! path temperature minus 250 k
!
! the 500-800 cm^-1 emission in emis(i,4) has been combined
!              into the 0-800 cm^-1 emission in emis(i,1)
!
   real(r8) emis(pcols,2)           ! h2o emissivity 
!
!
!
   real(r8) term7(pcols,2)          ! kl_inf(i) in eq(r8) of table a3a of r&d
   real(r8) term8(pcols,2)          ! delta kl_inf(i) in eq(r8)
   real(r8) tr1(pcols)              ! equation(6) in table a2 for 650-800
   real(r8) tr2(pcols)              ! equation(6) in table a2 for 500-650
   real(r8) tr3(pcols)              ! equation(4) in table a2 for 650-800
   real(r8) tr4(pcols)              ! equation(4),table a2 of r&d for 500-650
   real(r8) tr7(pcols)              ! equation (6) times eq(4) in table a2
!                                      of r&d for 650-800 cm-1 region
   real(r8) tr8(pcols)              ! equation (6) times eq(4) in table a2
!                                      of r&d for 500-650 cm-1 region
   real(r8) k21(pcols)              ! exponential coefficient used to calc
!                                     rot band transmissivity in the 650-800
!                                     cm-1 region (tr1)
   real(r8) k22(pcols)              ! exponential coefficient used to calc
!                                     rot band transmissivity in the 500-650
!                                     cm-1 region (tr2)
   real(r8) u(pcols)                ! pressure weighted h2o path length
   real(r8) ub(nbands)              ! pressure weighted h2o path length with
                                    !  hulst-curtis-godson correction for
                                    !  each band
   real(r8) pnew                    ! effective pressure for h2o linewidth
   real(r8) pnewb(nbands)           ! effective pressure for h2o linewidth w/
                                    !  hulst-curtis-godson correction for
                                    !  each band
   real(r8) uc1(pcols)              ! h2o continuum pathlength 500-800 cm-1
   real(r8) fwk                     ! equation(33) in r&d far wing correction
   real(r8) troco2(pcols,pverp)     ! h2o overlap factor for co2 absorption
   real(r8) emplnk(14,pcols)        ! emissivity planck factor
   real(r8) emstrc(pcols,pverp)     ! total trace gas emissivity
!
! local variables for co2:
!
   real(r8) co2ems(pcols,pverp)      ! co2 emissivity
   real(r8) co2plk(pcols)            ! used to compute co2 emissivity
   real(r8) sum(pcols)               ! used to calculate path temperature
   real(r8) t1i                      ! co2 hot band temperature factor
   real(r8) sqti                     ! sqrt of temperature
   real(r8) pi                       ! pressure used in co2 mean line width
   real(r8) et                       ! co2 hot band factor
   real(r8) et2                      ! co2 hot band factor
   real(r8) et4                      ! co2 hot band factor
   real(r8) omet                     ! co2 stimulated emission term
   real(r8) ex                       ! part of co2 planck function
   real(r8) f1co2                    ! co2 weak band factor
   real(r8) f2co2                    ! co2 weak band factor
   real(r8) f3co2                    ! co2 weak band factor
   real(r8) t1co2                    ! overlap factor weak bands strong band
   real(r8) sqwp                     ! sqrt of co2 pathlength
   real(r8) f1sqwp                   ! main co2 band factor
   real(r8) oneme                    ! co2 stimulated emission term
   real(r8) alphat                   ! part of the co2 stimulated emiss term
   real(r8) wco2                     ! consts used to define co2 pathlength
   real(r8) posqt                    ! effective pressure for co2 line width
   real(r8) rbeta7                   ! inverse of co2 hot band line width par
   real(r8) rbeta8                   ! inverse of co2 hot band line width par
   real(r8) rbeta9                   ! inverse of co2 hot band line width par
   real(r8) rbeta13                  ! inverse of co2 hot band line width par
   real(r8) tpath                    ! path temp used in co2 band model
   real(r8) tmp1                     ! co2 band factor
   real(r8) tmp2                     ! co2 band factor
   real(r8) tmp3                     ! co2 band factor
   real(r8) tlayr5                   ! temperature factor in co2 planck func
   real(r8) rsqti                    ! reciprocal of sqrt of temperature
   real(r8) exm1sq                   ! part of co2 planck function
   real(r8) u7                       ! absorber amt for various co2 band systems
   real(r8) u8                       ! absorber amt for various co2 band systems
   real(r8) u9                       ! absorber amt for various co2 band systems
   real(r8) u13                      ! absorber amt for various co2 band systems
   real(r8) r250                     ! inverse 250k
   real(r8) r300                     ! inverse 300k
   real(r8) rsslp                    ! inverse standard sea-level pressure
!
! local variables for o3:
!
   real(r8) o3ems(pcols,pverp)       ! ozone emissivity
   real(r8) dbvtt(pcols)             ! tmp drvtv of planck fctn for tplnke
   real(r8) dbvt,fo3,t,ux,vx
   real(r8) te                       ! temperature factor
   real(r8) u1                       ! path length factor
   real(r8) u2                       ! path length factor
   real(r8) phat                     ! effecitive path length pressure
   real(r8) tlocal                   ! local planck function temperature
   real(r8) tcrfac                   ! scaled temperature factor
   real(r8) beta                     ! absorption funct factor voigt effect
   real(r8) realnu                   ! absorption function factor
   real(r8) o3bndi                   ! band absorption factor
!
! transmission terms for various spectral intervals:
!
   real(r8) absbnd                   ! proportional to co2 band absorptance
   real(r8) tco2(pcols)              ! co2 overlap factor
   real(r8) th2o(pcols)              ! h2o overlap factor
   real(r8) to3(pcols)               ! o3 overlap factor
!
! variables for new h2o parameterization
!
! notation:
! u   = integral (p/p_0 dw)  eq. 15 in ramanathan/downey 1986
! p   = atmospheric pressure
! p_0 = reference atmospheric pressure
! w   = precipitable water path
! t_e = emission temperature
! t_p = path temperature
! rh  = path relative humidity
!
   real(r8) fe               ! asymptotic value of emis. as u->infinity
   real(r8) e_star           ! normalized non-window emissivity
   real(r8) l_star           ! interpolated line transmission
   real(r8) c_star           ! interpolated continuum transmission

   real(r8) te1              ! emission temperature
   real(r8) te2              ! te^2
   real(r8) te3              ! te^3
   real(r8) te4              ! te^4
   real(r8) te5              ! te^5

   real(r8) log_u            ! log base 10 of u 
   real(r8) log_uc           ! log base 10 of h2o continuum path
   real(r8) log_p            ! log base 10 of p
   real(r8) t_p              ! t_p
   real(r8) t_e              ! t_e (offset by t_p)

   integer iu                ! index for log10(u)
   integer iu1               ! iu + 1
   integer iuc               ! index for log10(h2o continuum path)
   integer iuc1              ! iuc + 1
   integer ip                ! index for log10(p)
   integer ip1               ! ip + 1
   integer itp               ! index for t_p
   integer itp1              ! itp + 1
   integer ite               ! index for t_e
   integer ite1              ! ite + 1
   integer irh               ! index for rh
   integer irh1              ! irh + 1

   real(r8) dvar             ! normalized variation in t_p/t_e/p/u
   real(r8) uvar             ! u * diffusivity factor
   real(r8) uscl             ! factor for lineary scaling as u->0

   real(r8) wu               ! weight for u
   real(r8) wu1              ! 1 - wu
   real(r8) wuc              ! weight for h2o continuum path
   real(r8) wuc1             ! 1 - wuc
   real(r8) wp               ! weight for p
   real(r8) wp1              ! 1 - wp
   real(r8) wtp              ! weight for t_p
   real(r8) wtp1             ! 1 - wtp
   real(r8) wte              ! weight for t_e
   real(r8) wte1             ! 1 - wte
   real(r8) wrh              ! weight for rh
   real(r8) wrh1             ! 1 - wrh

   real(r8) w_0_0_           ! weight for tp/te combination
   real(r8) w_0_1_           ! weight for tp/te combination
   real(r8) w_1_0_           ! weight for tp/te combination
   real(r8) w_1_1_           ! weight for tp/te combination

   real(r8) w_0_00           ! weight for tp/te/rh combination
   real(r8) w_0_01           ! weight for tp/te/rh combination
   real(r8) w_0_10           ! weight for tp/te/rh combination
   real(r8) w_0_11           ! weight for tp/te/rh combination
   real(r8) w_1_00           ! weight for tp/te/rh combination
   real(r8) w_1_01           ! weight for tp/te/rh combination
   real(r8) w_1_10           ! weight for tp/te/rh combination
   real(r8) w_1_11           ! weight for tp/te/rh combination

   real(r8) w00_00           ! weight for p/tp/te/rh combination
   real(r8) w00_01           ! weight for p/tp/te/rh combination
   real(r8) w00_10           ! weight for p/tp/te/rh combination
   real(r8) w00_11           ! weight for p/tp/te/rh combination
   real(r8) w01_00           ! weight for p/tp/te/rh combination
   real(r8) w01_01           ! weight for p/tp/te/rh combination
   real(r8) w01_10           ! weight for p/tp/te/rh combination
   real(r8) w01_11           ! weight for p/tp/te/rh combination
   real(r8) w10_00           ! weight for p/tp/te/rh combination
   real(r8) w10_01           ! weight for p/tp/te/rh combination
   real(r8) w10_10           ! weight for p/tp/te/rh combination
   real(r8) w10_11           ! weight for p/tp/te/rh combination
   real(r8) w11_00           ! weight for p/tp/te/rh combination
   real(r8) w11_01           ! weight for p/tp/te/rh combination
   real(r8) w11_10           ! weight for p/tp/te/rh combination
   real(r8) w11_11           ! weight for p/tp/te/rh combination

   integer ib                ! spectral interval:
                             !   1 = 0-800 cm^-1 and 1200-2200 cm^-1
                             !   2 = 800-1200 cm^-1

   real(r8) pch2o            ! h2o continuum path
   real(r8) fch2o            ! temp. factor for continuum
   real(r8) uch2o            ! u corresponding to h2o cont. path (window)

   real(r8) fdif             ! secant(zenith angle) for diffusivity approx.

   real(r8) sslp_mks         ! sea-level pressure in mks units
   real(r8) esx              ! saturation vapor pressure returned by vqsatd
   real(r8) qsx              ! saturation mixing ratio returned by vqsatd
   real(r8) pnew_mks         ! pnew in mks units
   real(r8) q_path           ! effective specific humidity along path
   real(r8) rh_path          ! effective relative humidity along path
   real(r8) omeps            ! 1 - epsilo

   integer  iest             ! index in estblh2o

!
!---------------------------statement functions-------------------------
!
! derivative of planck function at 9.6 micro-meter wavelength, and
! an absorption function factor:
!
!
   dbvt(t)=(-2.8911366682e-4+(2.3771251896e-6+1.1305188929e-10*t)*t)/ &
           (1.0+(-6.1364820707e-3+1.5550319767e-5*t)*t)
!
   fo3(ux,vx)=ux/sqrt(4.+ux*(1.+vx))
!
!
!
!-----------------------------------------------------------------------
!
! initialize
!
   r250  = 1./250.
   r300  = 1./300.
   rsslp = 1./sslp
!
! constants for computing u corresponding to h2o cont. path
!
   fdif       = 1.66
   sslp_mks   = sslp / 10.0
   omeps      = 1.0 - epsilo
!
! planck function for co2
!
   do i=1,ncol
      ex             = exp(960./tplnke(i))
      co2plk(i)      = 5.e8/((tplnke(i)**4)*(ex - 1.))
      co2t(i,ntoplw) = tplnke(i)
      sum(i)         = co2t(i,ntoplw)*pnm(i,ntoplw)
   end do
   k = ntoplw
   do k1=pverp,ntoplw+1,-1
      k = k + 1
      do i=1,ncol
         sum(i)         = sum(i) + tlayr(i,k)*(pnm(i,k)-pnm(i,k-1))
         ex             = exp(960./tlayr(i,k1))
         tlayr5         = tlayr(i,k1)*tlayr4(i,k1)
         co2eml(i,k1-1) = 1.2e11*ex/(tlayr5*(ex - 1.)**2)
         co2t(i,k)      = sum(i)/pnm(i,k)
      end do
   end do
!
! initialize planck function derivative for o3
!
   do i=1,ncol
      dbvtt(i) = dbvt(tplnke(i))
   end do
!
! calculate trace gas planck functions
!
   call trcplk(lchnk   ,ncol    ,pcols, pver, pverp,         &
               tint    ,tlayr   ,tplnke  ,emplnk  ,abplnk1 , &
               abplnk2 )
!
! interface loop
!
   do k1=ntoplw,pverp
!
! h2o emissivity
!
! emis(i,1)     0 -  800 cm-1   h2o rotation band
! emis(i,1)  1200 - 2200 cm-1   h2o vibration-rotation band
! emis(i,2)   800 - 1200 cm-1   h2o window
!
! separation between rotation and vibration-rotation dropped, so
!                only 2 slots needed for h2o emissivity
!
!      emis(i,3)   = 0.0
!
! for the p type continuum
!
      do i=1,ncol
         u(i)        = plh2o(i,k1)
         pnew        = u(i)/w(i,k1)
         pnew_mks    = pnew * sslp_mks
!
! apply scaling factor for 500-800 continuum
!
         uc1(i)      = (s2c(i,k1) + 1.7e-3*plh2o(i,k1))*(1. + 2.*s2c(i,k1))/ &
                       (1. + 15.*s2c(i,k1))
         pch2o       = s2c(i,k1)
!
! changed effective path temperature to std. curtis-godson form
!
         tpathe   = tcg(i,k1)/w(i,k1)
         t_p = min(max(tpathe, min_tp_h2o), max_tp_h2o)
         iest = floor(t_p) - min_tp_h2o
         esx = estblh2o(iest) + (estblh2o(iest+1)-estblh2o(iest)) * &
               (t_p - min_tp_h2o - iest)
         qsx = epsilo * esx / (pnew_mks - omeps * esx)
!
! compute effective rh along path
!
         q_path = w(i,k1) / pnm(i,k1) / rga
!
! calculate effective u, pnew for each band using
!        hulst-curtis-godson approximation:
! formulae: goody and yung, atmospheric radiation: theoretical basis, 
!           2nd edition, oxford university press, 1989.
! effective h2o path (w)
!      eq. 6.24, p. 228
! effective h2o path pressure (pnew = u/w):
!      eq. 6.29, p. 228
!
         ub(1) = plh2ob(1,i,k1) / psi(t_p,1)
         ub(2) = plh2ob(2,i,k1) / psi(t_p,2)

         pnewb(1) = ub(1) / wb(1,i,k1) * phi(t_p,1)
         pnewb(2) = ub(2) / wb(2,i,k1) * phi(t_p,2)
!
!
!
         dtx(i) = tplnke(i) - 250.
         dty(i) = tpathe - 250.
!
! define variables for c/h/e (now c/lt/e) fit
!
! emis(i,1)     0 -  800 cm-1   h2o rotation band
! emis(i,1)  1200 - 2200 cm-1   h2o vibration-rotation band
! emis(i,2)   800 - 1200 cm-1   h2o window
!
! separation between rotation and vibration-rotation dropped, so
!                only 2 slots needed for h2o emissivity
!
! emis(i,3)   = 0.0
!
! notation:
! u   = integral (p/p_0 dw)  
! p   = atmospheric pressure
! p_0 = reference atmospheric pressure
! w   = precipitable water path
! t_e = emission temperature
! t_p = path temperature
! rh  = path relative humidity
!
! terms for asymptotic value of emissivity
!
         te1  = tplnke(i)
         te2  = te1 * te1
         te3  = te2 * te1
         te4  = te3 * te1
         te5  = te4 * te1
!
! band-independent indices for lines and continuum tables
!
         dvar = (t_p - min_tp_h2o) / dtp_h2o
         itp = min(max(int(aint(dvar,r8)) + 1, 1), n_tp - 1)
         itp1 = itp + 1
         wtp = dvar - floor(dvar)
         wtp1 = 1.0 - wtp

         t_e = min(max(tplnke(i) - t_p, min_te_h2o), max_te_h2o)
         dvar = (t_e - min_te_h2o) / dte_h2o
         ite = min(max(int(aint(dvar,r8)) + 1, 1), n_te - 1)
         ite1 = ite + 1
         wte = dvar - floor(dvar)
         wte1 = 1.0 - wte

         rh_path = min(max(q_path / qsx, min_rh_h2o), max_rh_h2o)
         dvar = (rh_path - min_rh_h2o) / drh_h2o
         irh = min(max(int(aint(dvar,r8)) + 1, 1), n_rh - 1)
         irh1 = irh + 1
         wrh = dvar - floor(dvar)
         wrh1 = 1.0 - wrh

         w_0_0_ = wtp  * wte
         w_0_1_ = wtp  * wte1
         w_1_0_ = wtp1 * wte 
         w_1_1_ = wtp1 * wte1

         w_0_00 = w_0_0_ * wrh
         w_0_01 = w_0_0_ * wrh1
         w_0_10 = w_0_1_ * wrh
         w_0_11 = w_0_1_ * wrh1
         w_1_00 = w_1_0_ * wrh
         w_1_01 = w_1_0_ * wrh1
         w_1_10 = w_1_1_ * wrh
         w_1_11 = w_1_1_ * wrh1
!
! h2o continuum path for 0-800 and 1200-2200 cm^-1
!
!    assume foreign continuum dominates total h2o continuum in these bands
!    per clough et al, jgr, v. 97, no. d14 (oct 20, 1992), p. 15776
!    then the effective h2o path is just 
!         u_c = integral[ f(p) dw ]
!    where 
!           w = water-vapor mass and 
!        f(p) = dependence of foreign continuum on pressure 
!             = p / sslp
!    then 
!         u_c = u (the same effective h2o path as for lines)
!
!
! continuum terms for 800-1200 cm^-1
!
!    assume self continuum dominates total h2o continuum for this band
!    per clough et al, jgr, v. 97, no. d14 (oct 20, 1992), p. 15776
!    then the effective h2o self-continuum path is 
!         u_c = integral[ h(e,t) dw ]                        (*eq. 1*)
!    where 
!           w = water-vapor mass and 
!           e = partial pressure of h2o along path
!           t = temperature along path
!      h(e,t) = dependence of foreign continuum on e,t
!             = e / sslp * f(t)
!
!    replacing
!           e =~ q * p / epsilo
!           q = mixing ratio of h2o
!     epsilo = 0.622
!
!    and using the definition
!           u = integral [ (p / sslp) dw ]
!             = (p / sslp) w                                 (homogeneous path)
!
!    the effective path length for the self continuum is
!         u_c = (q / epsilo) f(t) u                         (*eq. 2*)
!
!    once values of t, u, and q have been calculated for the inhomogeneous
!        path, this sets u_c for the corresponding
!        homogeneous atmosphere.  however, this need not equal the
!        value of u_c' defined by eq. 1 for the actual inhomogeneous atmosphere
!        under consideration.
!
!    solution: hold t and q constant, solve for u' that gives u_c' by
!        inverting eq. (2):
!
!        u' = (u_c * epsilo) / (q * f(t))
!
         fch2o = fh2oself(t_p)
         uch2o = (pch2o * epsilo) / (q_path * fch2o)

!
! band-dependent indices for non-window
!
         ib = 1

         uvar = ub(ib) * fdif
         log_u  = min(log10(max(uvar, min_u_h2o)), max_lu_h2o)
         dvar = (log_u - min_lu_h2o) / dlu_h2o
         iu = min(max(int(aint(dvar,r8)) + 1, 1), n_u - 1)
         iu1 = iu + 1
         wu = dvar - floor(dvar)
         wu1 = 1.0 - wu
         
         log_p  = min(log10(max(pnewb(ib), min_p_h2o)), max_lp_h2o)
         dvar = (log_p - min_lp_h2o) / dlp_h2o
         ip = min(max(int(aint(dvar,r8)) + 1, 1), n_p - 1)
         ip1 = ip + 1
         wp = dvar - floor(dvar)
         wp1 = 1.0 - wp

         w00_00 = wp  * w_0_00 
         w00_01 = wp  * w_0_01 
         w00_10 = wp  * w_0_10 
         w00_11 = wp  * w_0_11 
         w01_00 = wp  * w_1_00 
         w01_01 = wp  * w_1_01 
         w01_10 = wp  * w_1_10 
         w01_11 = wp  * w_1_11 
         w10_00 = wp1 * w_0_00 
         w10_01 = wp1 * w_0_01 
         w10_10 = wp1 * w_0_10 
         w10_11 = wp1 * w_0_11 
         w11_00 = wp1 * w_1_00 
         w11_01 = wp1 * w_1_01 
         w11_10 = wp1 * w_1_10 
         w11_11 = wp1 * w_1_11 

!
! asymptotic value of emissivity as u->infinity
!
         fe = fet(1,ib) + &
              fet(2,ib) * te1 + &
              fet(3,ib) * te2 + &
              fet(4,ib) * te3 + &
              fet(5,ib) * te4 + &
              fet(6,ib) * te5

         e_star = &
              eh2onw(ip , itp , iu , ite , irh ) * w11_11 * wu1 + &
              eh2onw(ip , itp , iu , ite , irh1) * w11_10 * wu1 + &
              eh2onw(ip , itp , iu , ite1, irh ) * w11_01 * wu1 + &
              eh2onw(ip , itp , iu , ite1, irh1) * w11_00 * wu1 + &
              eh2onw(ip , itp , iu1, ite , irh ) * w11_11 * wu  + &
              eh2onw(ip , itp , iu1, ite , irh1) * w11_10 * wu  + &
              eh2onw(ip , itp , iu1, ite1, irh ) * w11_01 * wu  + &
              eh2onw(ip , itp , iu1, ite1, irh1) * w11_00 * wu  + &
              eh2onw(ip , itp1, iu , ite , irh ) * w10_11 * wu1 + &
              eh2onw(ip , itp1, iu , ite , irh1) * w10_10 * wu1 + &
              eh2onw(ip , itp1, iu , ite1, irh ) * w10_01 * wu1 + &
              eh2onw(ip , itp1, iu , ite1, irh1) * w10_00 * wu1 + &
              eh2onw(ip , itp1, iu1, ite , irh ) * w10_11 * wu  + &
              eh2onw(ip , itp1, iu1, ite , irh1) * w10_10 * wu  + &
              eh2onw(ip , itp1, iu1, ite1, irh ) * w10_01 * wu  + &
              eh2onw(ip , itp1, iu1, ite1, irh1) * w10_00 * wu  + &
              eh2onw(ip1, itp , iu , ite , irh ) * w01_11 * wu1 + &
              eh2onw(ip1, itp , iu , ite , irh1) * w01_10 * wu1 + &
              eh2onw(ip1, itp , iu , ite1, irh ) * w01_01 * wu1 + &
              eh2onw(ip1, itp , iu , ite1, irh1) * w01_00 * wu1 + &
              eh2onw(ip1, itp , iu1, ite , irh ) * w01_11 * wu  + &
              eh2onw(ip1, itp , iu1, ite , irh1) * w01_10 * wu  + &
              eh2onw(ip1, itp , iu1, ite1, irh ) * w01_01 * wu  + &
              eh2onw(ip1, itp , iu1, ite1, irh1) * w01_00 * wu  + &
              eh2onw(ip1, itp1, iu , ite , irh ) * w00_11 * wu1 + &
              eh2onw(ip1, itp1, iu , ite , irh1) * w00_10 * wu1 + &
              eh2onw(ip1, itp1, iu , ite1, irh ) * w00_01 * wu1 + &
              eh2onw(ip1, itp1, iu , ite1, irh1) * w00_00 * wu1 + &
              eh2onw(ip1, itp1, iu1, ite , irh ) * w00_11 * wu  + &
              eh2onw(ip1, itp1, iu1, ite , irh1) * w00_10 * wu  + &
              eh2onw(ip1, itp1, iu1, ite1, irh ) * w00_01 * wu  + &
              eh2onw(ip1, itp1, iu1, ite1, irh1) * w00_00 * wu 
         emis(i,ib) = min(max(fe * (1.0 - (1.0 - e_star) * &
                              aer_trn_ttl(i,k1,1,ib)), &
                          0.0_r8), 1.0_r8)
!
! invoke linear limit for scaling wrt u below min_u_h2o
!
         if (uvar < min_u_h2o) then
            uscl = uvar / min_u_h2o
            emis(i,ib) = emis(i,ib) * uscl
         endif

                      

!
! band-dependent indices for window
!
         ib = 2

         uvar = ub(ib) * fdif
         log_u  = min(log10(max(uvar, min_u_h2o)), max_lu_h2o)
         dvar = (log_u - min_lu_h2o) / dlu_h2o
         iu = min(max(int(aint(dvar,r8)) + 1, 1), n_u - 1)
         iu1 = iu + 1
         wu = dvar - floor(dvar)
         wu1 = 1.0 - wu
         
         log_p  = min(log10(max(pnewb(ib), min_p_h2o)), max_lp_h2o)
         dvar = (log_p - min_lp_h2o) / dlp_h2o
         ip = min(max(int(aint(dvar,r8)) + 1, 1), n_p - 1)
         ip1 = ip + 1
         wp = dvar - floor(dvar)
         wp1 = 1.0 - wp

         w00_00 = wp  * w_0_00 
         w00_01 = wp  * w_0_01 
         w00_10 = wp  * w_0_10 
         w00_11 = wp  * w_0_11 
         w01_00 = wp  * w_1_00 
         w01_01 = wp  * w_1_01 
         w01_10 = wp  * w_1_10 
         w01_11 = wp  * w_1_11 
         w10_00 = wp1 * w_0_00 
         w10_01 = wp1 * w_0_01 
         w10_10 = wp1 * w_0_10 
         w10_11 = wp1 * w_0_11 
         w11_00 = wp1 * w_1_00 
         w11_01 = wp1 * w_1_01 
         w11_10 = wp1 * w_1_10 
         w11_11 = wp1 * w_1_11 

         log_uc  = min(log10(max(uch2o * fdif, min_u_h2o)), max_lu_h2o)
         dvar = (log_uc - min_lu_h2o) / dlu_h2o
         iuc = min(max(int(aint(dvar,r8)) + 1, 1), n_u - 1)
         iuc1 = iuc + 1
         wuc = dvar - floor(dvar)
         wuc1 = 1.0 - wuc
!
! asymptotic value of emissivity as u->infinity
!
         fe = fet(1,ib) + &
              fet(2,ib) * te1 + &
              fet(3,ib) * te2 + &
              fet(4,ib) * te3 + &
              fet(5,ib) * te4 + &
              fet(6,ib) * te5

         l_star = &
              ln_eh2ow(ip , itp , iu , ite , irh ) * w11_11 * wu1 + &
              ln_eh2ow(ip , itp , iu , ite , irh1) * w11_10 * wu1 + &
              ln_eh2ow(ip , itp , iu , ite1, irh ) * w11_01 * wu1 + &
              ln_eh2ow(ip , itp , iu , ite1, irh1) * w11_00 * wu1 + &
              ln_eh2ow(ip , itp , iu1, ite , irh ) * w11_11 * wu  + &
              ln_eh2ow(ip , itp , iu1, ite , irh1) * w11_10 * wu  + &
              ln_eh2ow(ip , itp , iu1, ite1, irh ) * w11_01 * wu  + &
              ln_eh2ow(ip , itp , iu1, ite1, irh1) * w11_00 * wu  + &
              ln_eh2ow(ip , itp1, iu , ite , irh ) * w10_11 * wu1 + &
              ln_eh2ow(ip , itp1, iu , ite , irh1) * w10_10 * wu1 + &
              ln_eh2ow(ip , itp1, iu , ite1, irh ) * w10_01 * wu1 + &
              ln_eh2ow(ip , itp1, iu , ite1, irh1) * w10_00 * wu1 + &
              ln_eh2ow(ip , itp1, iu1, ite , irh ) * w10_11 * wu  + &
              ln_eh2ow(ip , itp1, iu1, ite , irh1) * w10_10 * wu  + &
              ln_eh2ow(ip , itp1, iu1, ite1, irh ) * w10_01 * wu  + &
              ln_eh2ow(ip , itp1, iu1, ite1, irh1) * w10_00 * wu  + &
              ln_eh2ow(ip1, itp , iu , ite , irh ) * w01_11 * wu1 + &
              ln_eh2ow(ip1, itp , iu , ite , irh1) * w01_10 * wu1 + &
              ln_eh2ow(ip1, itp , iu , ite1, irh ) * w01_01 * wu1 + &
              ln_eh2ow(ip1, itp , iu , ite1, irh1) * w01_00 * wu1 + &
              ln_eh2ow(ip1, itp , iu1, ite , irh ) * w01_11 * wu  + &
              ln_eh2ow(ip1, itp , iu1, ite , irh1) * w01_10 * wu  + &
              ln_eh2ow(ip1, itp , iu1, ite1, irh ) * w01_01 * wu  + &
              ln_eh2ow(ip1, itp , iu1, ite1, irh1) * w01_00 * wu  + &
              ln_eh2ow(ip1, itp1, iu , ite , irh ) * w00_11 * wu1 + &
              ln_eh2ow(ip1, itp1, iu , ite , irh1) * w00_10 * wu1 + &
              ln_eh2ow(ip1, itp1, iu , ite1, irh ) * w00_01 * wu1 + &
              ln_eh2ow(ip1, itp1, iu , ite1, irh1) * w00_00 * wu1 + &
              ln_eh2ow(ip1, itp1, iu1, ite , irh ) * w00_11 * wu  + &
              ln_eh2ow(ip1, itp1, iu1, ite , irh1) * w00_10 * wu  + &
              ln_eh2ow(ip1, itp1, iu1, ite1, irh ) * w00_01 * wu  + &
              ln_eh2ow(ip1, itp1, iu1, ite1, irh1) * w00_00 * wu 

         c_star = &
              cn_eh2ow(ip , itp , iuc , ite , irh ) * w11_11 * wuc1 + &
              cn_eh2ow(ip , itp , iuc , ite , irh1) * w11_10 * wuc1 + &
              cn_eh2ow(ip , itp , iuc , ite1, irh ) * w11_01 * wuc1 + &
              cn_eh2ow(ip , itp , iuc , ite1, irh1) * w11_00 * wuc1 + &
              cn_eh2ow(ip , itp , iuc1, ite , irh ) * w11_11 * wuc  + &
              cn_eh2ow(ip , itp , iuc1, ite , irh1) * w11_10 * wuc  + &
              cn_eh2ow(ip , itp , iuc1, ite1, irh ) * w11_01 * wuc  + &
              cn_eh2ow(ip , itp , iuc1, ite1, irh1) * w11_00 * wuc  + &
              cn_eh2ow(ip , itp1, iuc , ite , irh ) * w10_11 * wuc1 + &
              cn_eh2ow(ip , itp1, iuc , ite , irh1) * w10_10 * wuc1 + &
              cn_eh2ow(ip , itp1, iuc , ite1, irh ) * w10_01 * wuc1 + &
              cn_eh2ow(ip , itp1, iuc , ite1, irh1) * w10_00 * wuc1 + &
              cn_eh2ow(ip , itp1, iuc1, ite , irh ) * w10_11 * wuc  + &
              cn_eh2ow(ip , itp1, iuc1, ite , irh1) * w10_10 * wuc  + &
              cn_eh2ow(ip , itp1, iuc1, ite1, irh ) * w10_01 * wuc  + &
              cn_eh2ow(ip , itp1, iuc1, ite1, irh1) * w10_00 * wuc  + &
              cn_eh2ow(ip1, itp , iuc , ite , irh ) * w01_11 * wuc1 + &
              cn_eh2ow(ip1, itp , iuc , ite , irh1) * w01_10 * wuc1 + &
              cn_eh2ow(ip1, itp , iuc , ite1, irh ) * w01_01 * wuc1 + &
              cn_eh2ow(ip1, itp , iuc , ite1, irh1) * w01_00 * wuc1 + &
              cn_eh2ow(ip1, itp , iuc1, ite , irh ) * w01_11 * wuc  + &
              cn_eh2ow(ip1, itp , iuc1, ite , irh1) * w01_10 * wuc  + &
              cn_eh2ow(ip1, itp , iuc1, ite1, irh ) * w01_01 * wuc  + &
              cn_eh2ow(ip1, itp , iuc1, ite1, irh1) * w01_00 * wuc  + &
              cn_eh2ow(ip1, itp1, iuc , ite , irh ) * w00_11 * wuc1 + &
              cn_eh2ow(ip1, itp1, iuc , ite , irh1) * w00_10 * wuc1 + &
              cn_eh2ow(ip1, itp1, iuc , ite1, irh ) * w00_01 * wuc1 + &
              cn_eh2ow(ip1, itp1, iuc , ite1, irh1) * w00_00 * wuc1 + &
              cn_eh2ow(ip1, itp1, iuc1, ite , irh ) * w00_11 * wuc  + &
              cn_eh2ow(ip1, itp1, iuc1, ite , irh1) * w00_10 * wuc  + &
              cn_eh2ow(ip1, itp1, iuc1, ite1, irh ) * w00_01 * wuc  + &
              cn_eh2ow(ip1, itp1, iuc1, ite1, irh1) * w00_00 * wuc 
         emis(i,ib) = min(max(fe * (1.0 - l_star * c_star * &
                              aer_trn_ttl(i,k1,1,ib)), &
                          0.0_r8), 1.0_r8) 
!
! invoke linear limit for scaling wrt u below min_u_h2o
!
         if (uvar < min_u_h2o) then
            uscl = uvar / min_u_h2o
            emis(i,ib) = emis(i,ib) * uscl
         endif

                      
!
! compute total emissivity for h2o
!
         h2oems(i,k1) = emis(i,1)+emis(i,2)

      end do
!
!
!

      do i=1,ncol
         term7(i,1) = coefj(1,1) + coefj(2,1)*dty(i)*(1.+c16*dty(i))
         term8(i,1) = coefk(1,1) + coefk(2,1)*dty(i)*(1.+c17*dty(i))
         term7(i,2) = coefj(1,2) + coefj(2,2)*dty(i)*(1.+c26*dty(i))
         term8(i,2) = coefk(1,2) + coefk(2,2)*dty(i)*(1.+c27*dty(i))
      end do
      do i=1,ncol
!
! 500 -  800 cm-1   rotation band overlap with co2
!
         k21(i) = term7(i,1) + term8(i,1)/ &
                 (1. + (c30 + c31*(dty(i)-10.)*(dty(i)-10.))*sqrt(u(i)))
         k22(i) = term7(i,2) + term8(i,2)/ &
                 (1. + (c28 + c29*(dty(i)-10.))*sqrt(u(i)))
         fwk    = fwcoef + fwc1/(1.+fwc2*u(i))
         tr1(i) = exp(-(k21(i)*(sqrt(u(i)) + fc1*fwk*u(i))))
         tr2(i) = exp(-(k22(i)*(sqrt(u(i)) + fc1*fwk*u(i))))
         tr1(i)=tr1(i)*aer_trn_ttl(i,k1,1,idx_lw_0650_0800) 
!                                            ! h2o line+aer trn 650--800 cm-1
         tr2(i)=tr2(i)*aer_trn_ttl(i,k1,1,idx_lw_0500_0650) 
!                                            ! h2o line+aer trn 500--650 cm-1
         tr3(i) = exp(-((coefh(1,1) + coefh(2,1)*dtx(i))*uc1(i)))
         tr4(i) = exp(-((coefh(1,2) + coefh(2,2)*dtx(i))*uc1(i)))
         tr7(i) = tr1(i)*tr3(i)
         tr8(i) = tr2(i)*tr4(i)
         troco2(i,k1) = 0.65*tr7(i) + 0.35*tr8(i)
         th2o(i) = tr8(i)
      end do
!
! co2 emissivity for 15 micron band system
!
      do i=1,ncol
         t1i    = exp(-480./co2t(i,k1))
         sqti   = sqrt(co2t(i,k1))
         rsqti  = 1./sqti
         et     = t1i
         et2    = et*et
         et4    = et2*et2
         omet   = 1. - 1.5*et2
         f1co2  = 899.70*omet*(1. + 1.94774*et + 4.73486*et2)*rsqti
         sqwp   = sqrt(plco2(i,k1))
         f1sqwp = f1co2*sqwp
         t1co2  = 1./(1. + 245.18*omet*sqwp*rsqti)
         oneme  = 1. - et2
         alphat = oneme**3*rsqti
         wco2   = 2.5221*co2vmr*pnm(i,k1)*rga
         u7     = 4.9411e4*alphat*et2*wco2
         u8     = 3.9744e4*alphat*et4*wco2
         u9     = 1.0447e5*alphat*et4*et2*wco2
         u13    = 2.8388e3*alphat*et4*wco2
!
         tpath  = co2t(i,k1)
         tlocal = tplnke(i)
         tcrfac = sqrt((tlocal*r250)*(tpath*r300))
         pi     = pnm(i,k1)*rsslp + 2.*dpfco2*tcrfac
         posqt  = pi/(2.*sqti)
         rbeta7 =  1./( 5.3288*posqt)
         rbeta8 = 1./ (10.6576*posqt)
         rbeta9 = rbeta7
         rbeta13= rbeta9
         f2co2  = (u7/sqrt(4. + u7*(1. + rbeta7))) + &
                  (u8/sqrt(4. + u8*(1. + rbeta8))) + &
                  (u9/sqrt(4. + u9*(1. + rbeta9)))
         f3co2  = u13/sqrt(4. + u13*(1. + rbeta13))
         tmp1   = log(1. + f1sqwp)
         tmp2   = log(1. +  f2co2)
         tmp3   = log(1. +  f3co2)
         absbnd = (tmp1 + 2.*t1co2*tmp2 + 2.*tmp3)*sqti
         tco2(i)=1.0/(1.0+10.0*(u7/sqrt(4. + u7*(1. + rbeta7))))
         co2ems(i,k1)  = troco2(i,k1)*absbnd*co2plk(i)
         ex     = exp(960./tint(i,k1))
         exm1sq = (ex - 1.)**2
         co2em(i,k1) = 1.2e11*ex/(tint(i,k1)*tint4(i,k1)*exm1sq)
      end do
!
! o3 emissivity
!
      do i=1,ncol
         h2otr(i,k1) = exp(-12.*s2c(i,k1))
          h2otr(i,k1)=h2otr(i,k1)*aer_trn_ttl(i,k1,1,idx_lw_1000_1200)
         te          = (co2t(i,k1)/293.)**.7
         u1          = 18.29*plos(i,k1)/te
         u2          = .5649*plos(i,k1)/te
         phat        = plos(i,k1)/plol(i,k1)
         tlocal      = tplnke(i)
         tcrfac      = sqrt(tlocal*r250)*te
         beta        = (1./.3205)*((1./phat) + (dpfo3*tcrfac))
         realnu      = (1./beta)*te
         o3bndi      = 74.*te*(tplnke(i)/375.)*log(1. + fo3(u1,realnu) + fo3(u2,realnu))
         o3ems(i,k1) = dbvtt(i)*h2otr(i,k1)*o3bndi
         to3(i)=1.0/(1. + 0.1*fo3(u1,realnu) + 0.1*fo3(u2,realnu))
      end do
!
!   calculate trace gas emissivities
!
      call trcems(lchnk   ,ncol    ,pcols, pverp,               &
                  k1      ,co2t    ,pnm     ,ucfc11  ,ucfc12  , &
                  un2o0   ,un2o1   ,bn2o0   ,bn2o1   ,uch4    , &
                  bch4    ,uco211  ,uco212  ,uco213  ,uco221  , &
                  uco222  ,uco223  ,uptype  ,w       ,s2c     , &
                  u       ,emplnk  ,th2o    ,tco2    ,to3     , &
                  emstrc  , &
                  aer_trn_ttl)
!
! total emissivity:
!
      do i=1,ncol
         emstot(i,k1) = h2oems(i,k1) + co2ems(i,k1) + o3ems(i,k1)  &
                        + emstrc(i,k1)
      end do
   end do ! end of interface loop

   return
end subroutine radems

subroutine radtpl(lchnk   ,ncol    ,pcols, pver, pverp,                 &
                  tnm     ,lwupcgs ,qnm     ,pnm     ,plco2   ,plh2o   , &
                  tplnka  ,s2c     ,tcg     ,w       ,tplnke  , &
                  tint    ,tint4   ,tlayr   ,tlayr4  ,pmln    , &
                  piln    ,plh2ob  ,wb      )
!--------------------------------------------------------------------
!
! purpose:
! compute temperatures and path lengths for longwave radiation
!
! method:
! <describe the algorithm(s) used in the routine.>
! <also include any applicable external references.>
!
! author: ccm1
!
!--------------------------------------------------------------------

!------------------------------arguments-----------------------------
!
! input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns
   integer, intent(in) :: pcols, pver, pverp

   real(r8), intent(in) :: tnm(pcols,pver)      ! model level temperatures
   real(r8), intent(in) :: lwupcgs(pcols)       ! surface longwave up flux
   real(r8), intent(in) :: qnm(pcols,pver)      ! model level specific humidity
   real(r8), intent(in) :: pnm(pcols,pverp)     ! pressure at model interfaces (dynes/cm2)
   real(r8), intent(in) :: pmln(pcols,pver)     ! ln(pmidm1)
   real(r8), intent(in) :: piln(pcols,pverp)    ! ln(pintm1)
!
! output arguments
!
   real(r8), intent(out) :: plco2(pcols,pverp)   ! pressure weighted co2 path
   real(r8), intent(out) :: plh2o(pcols,pverp)   ! pressure weighted h2o path
   real(r8), intent(out) :: tplnka(pcols,pverp)  ! level temperature from interface temperatures
   real(r8), intent(out) :: s2c(pcols,pverp)     ! h2o continuum path length
   real(r8), intent(out) :: tcg(pcols,pverp)     ! h2o-mass-wgted temp. (curtis-godson approx.)
   real(r8), intent(out) :: w(pcols,pverp)       ! h2o path length
   real(r8), intent(out) :: tplnke(pcols)        ! equal to tplnka
   real(r8), intent(out) :: tint(pcols,pverp)    ! layer interface temperature
   real(r8), intent(out) :: tint4(pcols,pverp)   ! tint to the 4th power
   real(r8), intent(out) :: tlayr(pcols,pverp)   ! k-1 level temperature
   real(r8), intent(out) :: tlayr4(pcols,pverp)  ! tlayr to the 4th power
   real(r8), intent(out) :: plh2ob(nbands,pcols,pverp)! pressure weighted h2o path with 
                                                      !    hulst-curtis-godson temp. factor 
                                                      !    for h2o bands 
   real(r8), intent(out) :: wb(nbands,pcols,pverp)    ! h2o path length with 
                                                      !    hulst-curtis-godson temp. factor 
                                                      !    for h2o bands 

!
!---------------------------local variables--------------------------
!
   integer i                 ! longitude index
   integer k                 ! level index
   integer kp1               ! level index + 1

   real(r8) repsil               ! inver ratio mol weight h2o to dry air
   real(r8) dy                   ! thickness of layer for tmp interp
   real(r8) dpnm                 ! pressure thickness of layer
   real(r8) dpnmsq               ! prs squared difference across layer
   real(r8) dw                   ! increment in h2o path length
   real(r8) dplh2o               ! increment in plh2o
   real(r8) cpwpl                ! const in co2 mix ratio to path length conversn

!--------------------------------------------------------------------
!
   repsil = 1./epsilo
!
! compute co2 and h2o paths
!
   cpwpl = amco2/amd * 0.5/(gravit*p0)
   do i=1,ncol
      plh2o(i,ntoplw)  = rgsslp*qnm(i,ntoplw)*pnm(i,ntoplw)*pnm(i,ntoplw)
      plco2(i,ntoplw)  = co2vmr*cpwpl*pnm(i,ntoplw)*pnm(i,ntoplw)
   end do
   do k=ntoplw,pver
      do i=1,ncol
         plh2o(i,k+1)  = plh2o(i,k) + rgsslp* &
                         (pnm(i,k+1)**2 - pnm(i,k)**2)*qnm(i,k)
         plco2(i,k+1)  = co2vmr*cpwpl*pnm(i,k+1)**2
      end do
   end do
!
! set the top and bottom intermediate level temperatures,
! top level planck temperature and top layer temp**4.
!
! tint is lower interface temperature
! (not available for bottom layer, so use ground temperature)
!
   do i=1,ncol
      tint4(i,pverp)   = lwupcgs(i)/stebol
      tint(i,pverp)    = sqrt(sqrt(tint4(i,pverp)))
      tplnka(i,ntoplw) = tnm(i,ntoplw)
      tint(i,ntoplw)   = tplnka(i,ntoplw)
      tlayr4(i,ntoplw) = tplnka(i,ntoplw)**4
      tint4(i,ntoplw)  = tlayr4(i,ntoplw)
   end do
!
! intermediate level temperatures are computed using temperature
! at the full level below less dy*delta t,between the full level
!
   do k=ntoplw+1,pver
      do i=1,ncol
         dy = (piln(i,k) - pmln(i,k))/(pmln(i,k-1) - pmln(i,k))
         tint(i,k)  = tnm(i,k) - dy*(tnm(i,k)-tnm(i,k-1))
         tint4(i,k) = tint(i,k)**4
      end do
   end do
!
! now set the layer temp=full level temperatures and establish a
! planck temperature for absorption (tplnka) which is the average
! the intermediate level temperatures.  note that tplnka is not
! equal to the full level temperatures.
!
   do k=ntoplw+1,pverp
      do i=1,ncol
         tlayr(i,k)  = tnm(i,k-1)
         tlayr4(i,k) = tlayr(i,k)**4
         tplnka(i,k) = .5*(tint(i,k) + tint(i,k-1))
      end do
   end do
!
! calculate tplank for emissivity calculation.
! assume isothermal tplnke i.e. all levels=ttop.
!
   do i=1,ncol
      tplnke(i)       = tplnka(i,ntoplw)
      tlayr(i,ntoplw) = tint(i,ntoplw)
   end do
!
! now compute h2o path fields:
!
   do i=1,ncol
!
! changed effective path temperature to std. curtis-godson form
!
      tcg(i,ntoplw) = rga*qnm(i,ntoplw)*pnm(i,ntoplw)*tnm(i,ntoplw)
      w(i,ntoplw)   = sslp * (plh2o(i,ntoplw)*2.) / pnm(i,ntoplw)
!
! hulst-curtis-godson scaling for h2o path
!
      wb(1,i,ntoplw) = w(i,ntoplw) * phi(tnm(i,ntoplw),1)
      wb(2,i,ntoplw) = w(i,ntoplw) * phi(tnm(i,ntoplw),2)
!
! hulst-curtis-godson scaling for effective pressure along h2o path
!
      plh2ob(1,i,ntoplw) = plh2o(i,ntoplw) * psi(tnm(i,ntoplw),1)
      plh2ob(2,i,ntoplw) = plh2o(i,ntoplw) * psi(tnm(i,ntoplw),2)

      s2c(i,ntoplw) = plh2o(i,ntoplw)*fh2oself(tnm(i,ntoplw))*qnm(i,ntoplw)*repsil
   end do

   do k=ntoplw,pver
      do i=1,ncol
         dpnm       = pnm(i,k+1) - pnm(i,k)
         dpnmsq     = pnm(i,k+1)**2 - pnm(i,k)**2
         dw         = rga*qnm(i,k)*dpnm
         kp1        = k+1
         w(i,kp1)   = w(i,k) + dw
!
! hulst-curtis-godson scaling for h2o path
!
         wb(1,i,kp1) = wb(1,i,k) + dw * phi(tnm(i,k),1)
         wb(2,i,kp1) = wb(2,i,k) + dw * phi(tnm(i,k),2)
!
! hulst-curtis-godson scaling for effective pressure along h2o path
!
         dplh2o = plh2o(i,kp1) - plh2o(i,k)

         plh2ob(1,i,kp1) = plh2ob(1,i,k) + dplh2o * psi(tnm(i,k),1)
         plh2ob(2,i,kp1) = plh2ob(2,i,k) + dplh2o * psi(tnm(i,k),2)
!
! changed effective path temperature to std. curtis-godson form
!
         tcg(i,kp1) = tcg(i,k) + dw*tnm(i,k)
         s2c(i,kp1) = s2c(i,k) + rgsslp*dpnmsq*qnm(i,k)* &
                      fh2oself(tnm(i,k))*qnm(i,k)*repsil
      end do
   end do
!
   return
end subroutine radtpl


subroutine radclwmx(lchnk   ,ncol    ,pcols, pver, pverp,         &
                    lwupcgs ,tnm     ,qnm     ,o3vmr   , &
                    pmid    ,pint    ,pmln    ,piln    ,          &
                             n2o     ,ch4     ,cfc11   ,cfc12   , &
                    cld     ,emis    ,pmxrgn  ,nmxrgn  ,qrl     , &
                    doabsems, abstot, absnxt, emstot,             &
                    flns    ,flnt    ,flnsc   ,flntc   ,flwds   , &
                    flut    ,flutc   , &
                    flup    ,flupc   ,fldn    ,fldnc   ,          &
                    aer_mass)
!----------------------------------------------------------------------- 
! 
! purpose: 
! compute longwave radiation heating rates and boundary fluxes
! 
! method: 
! uses broad band absorptivity/emissivity method to compute clear sky;
! assumes randomly overlapped clouds with variable cloud emissivity to
! include effects of clouds.
!
! computes clear sky absorptivity/emissivity at lower frequency (in
! general) than the model radiation frequency; uses previously computed
! and stored values for efficiency
!
! note: this subroutine contains vertical indexing which proceeds
!       from bottom to top rather than the top to bottom indexing
!       used in the rest of the model.
! 
! author: b. collins
! 
!-----------------------------------------------------------------------
!  use shr_kind_mod, only: r8 => shr_kind_r8
!  use ppgrid
!  use radae, only: nbands, radems, radabs, radtpl, abstot_3d, absnxt_3d, emstot_3d
!  use volcrad

   implicit none

   integer pverp2,pverp3,pverp4
!  parameter (pverp2=pver+2,pverp3=pver+3,pverp4=pver+4)

   real(r8) cldmin
   parameter (cldmin = 1.0d-80)
!------------------------------commons----------------------------------
!-----------------------------------------------------------------------
!------------------------------arguments--------------------------------
!
! input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: pcols, pver, pverp
   integer, intent(in) :: ncol                  ! number of atmospheric columns
!    maximally overlapped region.
!    0->pmxrgn(i,1) is range of pmid for
!    1st region, pmxrgn(i,1)->pmxrgn(i,2) for
!    2nd region, etc
   integer, intent(in) :: nmxrgn(pcols)         ! number of maximally overlapped regions
   logical, intent(in) :: doabsems

   real(r8), intent(in) :: pmxrgn(pcols,pverp)  ! maximum values of pmid for each
   real(r8), intent(in) :: lwupcgs(pcols)       ! longwave up flux in cgs units
!
! input arguments which are only passed to other routines
!
   real(r8), intent(in) :: tnm(pcols,pver)      ! level temperature
   real(r8), intent(in) :: qnm(pcols,pver)      ! level moisture field
   real(r8), intent(in) :: o3vmr(pcols,pver)    ! ozone volume mixing ratio
   real(r8), intent(in) :: pmid(pcols,pver)     ! level pressure
   real(r8), intent(in) :: pint(pcols,pverp)    ! model interface pressure
   real(r8), intent(in) :: pmln(pcols,pver)     ! ln(pmid)
   real(r8), intent(in) :: piln(pcols,pverp)    ! ln(pint)
   real(r8), intent(in) :: n2o(pcols,pver)      ! nitrous oxide mass mixing ratio
   real(r8), intent(in) :: ch4(pcols,pver)      ! methane mass mixing ratio
   real(r8), intent(in) :: cfc11(pcols,pver)    ! cfc11 mass mixing ratio
   real(r8), intent(in) :: cfc12(pcols,pver)    ! cfc12 mass mixing ratio
   real(r8), intent(in) :: cld(pcols,pver)      ! cloud cover
   real(r8), intent(in) :: emis(pcols,pver)     ! cloud emissivity
   real(r8), intent(in) :: aer_mass(pcols,pver) ! straer mass in layer

!
! output arguments
!
   real(r8), intent(out) :: qrl(pcols,pver)      ! longwave heating rate
   real(r8), intent(out) :: flns(pcols)          ! surface cooling flux
   real(r8), intent(out) :: flnt(pcols)          ! net outgoing flux
   real(r8), intent(out) :: flut(pcols)          ! upward flux at top of model
   real(r8), intent(out) :: flnsc(pcols)         ! clear sky surface cooing
   real(r8), intent(out) :: flntc(pcols)         ! net clear sky outgoing flux
   real(r8), intent(out) :: flutc(pcols)         ! upward clear-sky flux at top of model
   real(r8), intent(out) :: flwds(pcols)         ! down longwave flux at surface
! added downward/upward total and clear sky fluxes
   real(r8), intent(out) :: flup(pcols,pverp)      ! total sky upward longwave flux 
   real(r8), intent(out) :: flupc(pcols,pverp)     ! clear sky upward longwave flux 
   real(r8), intent(out) :: fldn(pcols,pverp)      ! total sky downward longwave flux 
   real(r8), intent(out) :: fldnc(pcols,pverp)     ! clear sky downward longwave flux
!
   real(r8), intent(inout) :: abstot(pcols,pverp,pverp) ! total absorptivity
   real(r8), intent(inout) :: absnxt(pcols,pver,4)      ! total nearest layer absorptivity
   real(r8), intent(inout) :: emstot(pcols,pverp)     ! total emissivity

!---------------------------local variables-----------------------------
!
   integer i                 ! longitude index
   integer ilon              ! longitude index
   integer ii                ! longitude index
   integer iimx              ! longitude index (max overlap)
   integer k                 ! level index
   integer k1                ! level index
   integer k2                ! level index
   integer k3                ! level index
   integer km                ! level index
   integer km1               ! level index
   integer km3               ! level index
   integer km4               ! level index
   integer irgn              ! index for max-overlap regions
   integer l                 ! index for clouds to overlap
   integer l1                ! index for clouds to overlap
   integer n                 ! counter

!
   real(r8) :: plco2(pcols,pverp)   ! path length co2
   real(r8) :: plh2o(pcols,pverp)   ! path length h2o
   real(r8) tmp(pcols)           ! temporary workspace
   real(r8) tmp2(pcols)          ! temporary workspace
   real(r8) absbt(pcols)         ! downward emission at model top
   real(r8) plol(pcols,pverp)    ! o3 pressure wghted path length
   real(r8) plos(pcols,pverp)    ! o3 path length
   real(r8) aer_mpp(pcols,pverp) ! straer path above kth interface level
   real(r8) co2em(pcols,pverp)   ! layer co2 normalized planck funct. derivative
   real(r8) co2eml(pcols,pver)   ! interface co2 normalized planck funct. deriv.
   real(r8) delt(pcols)          ! diff t**4 mid layer to top interface
   real(r8) delt1(pcols)         ! diff t**4 lower intrfc to mid layer
   real(r8) bk1(pcols)           ! absrptvty for vertical quadrature
   real(r8) bk2(pcols)           ! absrptvty for vertical quadrature
   real(r8) cldp(pcols,pverp)    ! cloud cover with extra layer
   real(r8) ful(pcols,pverp)     ! total upwards longwave flux
   real(r8) fsul(pcols,pverp)    ! clear sky upwards longwave flux
   real(r8) fdl(pcols,pverp)     ! total downwards longwave flux
   real(r8) fsdl(pcols,pverp)    ! clear sky downwards longwv flux
   real(r8) fclb4(pcols,-1:pver)    ! sig t**4 for cld bottom interfc
   real(r8) fclt4(pcols,0:pver)    ! sig t**4 for cloud top interfc
   real(r8) s(pcols,pverp,pverp) ! flx integral sum
   real(r8) tplnka(pcols,pverp)  ! planck fnctn temperature
   real(r8) s2c(pcols,pverp)     ! h2o cont amount
   real(r8) tcg(pcols,pverp)     ! h2o-mass-wgted temp. (curtis-godson approx.)
   real(r8) w(pcols,pverp)       ! h2o path
   real(r8) tplnke(pcols)        ! planck fnctn temperature
   real(r8) h2otr(pcols,pverp)   ! h2o trnmsn for o3 overlap
   real(r8) co2t(pcols,pverp)    ! prs wghted temperature path
   real(r8) tint(pcols,pverp)    ! interface temperature
   real(r8) tint4(pcols,pverp)   ! interface temperature**4
   real(r8) tlayr(pcols,pverp)   ! level temperature
   real(r8) tlayr4(pcols,pverp)  ! level temperature**4
   real(r8) plh2ob(nbands,pcols,pverp)! pressure weighted h2o path with 
                                      !    hulst-curtis-godson temp. factor 
                                      !    for h2o bands 
   real(r8) wb(nbands,pcols,pverp)    ! h2o path length with 
                                      !    hulst-curtis-godson temp. factor 
                                      !    for h2o bands 

   real(r8) cld0                 ! previous cloud amt (for max overlap)
   real(r8) cld1                 ! next cloud amt (for max overlap)
   real(r8) emx(0:pverp)         ! emissivity factors (max overlap)
   real(r8) emx0                 ! emissivity factors for bcs (max overlap)
   real(r8) trans                ! 1 - emis
   real(r8) asort(pver)          ! 1 - cloud amounts to be sorted for max ovrlp.
   real(r8) atmp                 ! temporary storage for sort when nxs = 2
   real(r8) maxcld(pcols)        ! maximum cloud at any layer

   integer indx(pcols)       ! index vector of gathered array values
!!$   integer indxmx(pcols+1,pverp)! index vector of gathered array values
   integer indxmx(pcols,pverp)! index vector of gathered array values
!    (max overlap)
   integer nrgn(pcols)       ! number of max overlap regions at longitude
   integer npts              ! number of values satisfying some criterion
   integer ncolmx(pverp)     ! number of columns with clds in region
   integer kx1(pcols,pverp)  ! level index for top of max-overlap region
   integer kx2(pcols,0:pverp)! level index for bottom of max-overlap region
   integer kxs(0:pverp,pcols,pverp)! level indices for cld layers sorted by cld()
!    in descending order
   integer nxs(pcols,pverp)  ! number of cloudy layers between kx1 and kx2
   integer nxsk              ! number of cloudy layers between (kx1/kx2)&k
   integer ksort(0:pverp)    ! level indices of cloud amounts to be sorted
!    for max ovrlp. calculation
   integer ktmp              ! temporary storage for sort when nxs = 2

!  real aer_trn_ttl(pcols,pverp,pverp,bnd_nbr_lw) ! [fraction] total
  real(r8) aer_trn_ttl(pcols,pverp,pverp,bnd_nbr_lw) ! [fraction] total
!                               ! transmission between interfaces k1 and k2  
!
! pointer variables to 3d structures
!
!  real(r8), pointer :: abstot(:,:,:)
!  real(r8), pointer :: absnxt(:,:,:)
!  real(r8), pointer :: emstot(:,:)

!
! trace gas variables
!
   real(r8) ucfc11(pcols,pverp)  ! cfc11 path length
   real(r8) ucfc12(pcols,pverp)  ! cfc12 path length
   real(r8) un2o0(pcols,pverp)   ! n2o path length
   real(r8) un2o1(pcols,pverp)   ! n2o path length (hot band)
   real(r8) uch4(pcols,pverp)    ! ch4 path length
   real(r8) uco211(pcols,pverp)  ! co2 9.4 micron band path length
   real(r8) uco212(pcols,pverp)  ! co2 9.4 micron band path length
   real(r8) uco213(pcols,pverp)  ! co2 9.4 micron band path length
   real(r8) uco221(pcols,pverp)  ! co2 10.4 micron band path length
   real(r8) uco222(pcols,pverp)  ! co2 10.4 micron band path length
   real(r8) uco223(pcols,pverp)  ! co2 10.4 micron band path length
   real(r8) bn2o0(pcols,pverp)   ! pressure factor for n2o
   real(r8) bn2o1(pcols,pverp)   ! pressure factor for n2o
   real(r8) bch4(pcols,pverp)    ! pressure factor for ch4
   real(r8) uptype(pcols,pverp)  ! p-type continuum path length
   real(r8) abplnk1(14,pcols,pverp)  ! non-nearest layer plack factor
   real(r8) abplnk2(14,pcols,pverp)  ! nearest layer factor
!
!
!-----------------------------------------------------------------------
!
!
   pverp2=pver+2
   pverp3=pver+3
   pverp4=pver+4
!
! set pointer variables
!
!  abstot => abstot_3d(:,:,:,lchnk)
!  absnxt => absnxt_3d(:,:,:,lchnk)
!  emstot => emstot_3d(:,:,lchnk)
!
! accumulate mass path from top of atmosphere
!
  call aer_pth(aer_mass, aer_mpp, ncol, pcols, pver, pverp)

!
! calculate some temperatures needed to derive absorptivity and
! emissivity, as well as some h2o path lengths
!
   call radtpl(lchnk   ,ncol    ,pcols, pver, pverp,                  &
               tnm     ,lwupcgs ,qnm     ,pint    ,plco2   ,plh2o   , &
               tplnka  ,s2c     ,tcg     ,w       ,tplnke  , &
               tint    ,tint4   ,tlayr   ,tlayr4  ,pmln    , &
               piln    ,plh2ob  ,wb      )
   if (doabsems) then
!
! compute ozone path lengths at frequency of a/e calculation.
!
      call radoz2(lchnk, ncol, pcols, pver, pverp, o3vmr   ,pint    ,plol    ,plos, ntoplw    )
!
! compute trace gas path lengths
!
      call trcpth(lchnk   ,ncol    ,pcols, pver, pverp,         &
                  tnm     ,pint    ,cfc11   ,cfc12   ,n2o     , &
                  ch4     ,qnm     ,ucfc11  ,ucfc12  ,un2o0   , &
                  un2o1   ,uch4    ,uco211  ,uco212  ,uco213  , &
                  uco221  ,uco222  ,uco223  ,bn2o0   ,bn2o1   , &
                  bch4    ,uptype  )

!     compute transmission through straer absorption continuum
      call aer_trn(aer_mpp, aer_trn_ttl, pcols, pver, pverp)

!
!
! compute total emissivity:
!
      call radems(lchnk   ,ncol    ,pcols, pver, pverp,         &
                  s2c     ,tcg     ,w       ,tplnke  ,plh2o   , &
                  pint    ,plco2   ,tint    ,tint4   ,tlayr   , &
                  tlayr4  ,plol    ,plos    ,ucfc11  ,ucfc12  , &
                  un2o0   ,un2o1   ,uch4    ,uco211  ,uco212  , &
                  uco213  ,uco221  ,uco222  ,uco223  ,uptype  , &
                  bn2o0   ,bn2o1   ,bch4    ,co2em   ,co2eml  , &
                  co2t    ,h2otr   ,abplnk1 ,abplnk2 ,emstot  , &
                  plh2ob  ,wb      , &
                  aer_trn_ttl)
!
! compute total absorptivity:
!
      call radabs(lchnk   ,ncol    ,pcols, pver, pverp,         &
                  pmid    ,pint    ,co2em   ,co2eml  ,tplnka  , &
                  s2c     ,tcg     ,w       ,h2otr   ,plco2   , &
                  plh2o   ,co2t    ,tint    ,tlayr   ,plol    , &
                  plos    ,pmln    ,piln    ,ucfc11  ,ucfc12  , &
                  un2o0   ,un2o1   ,uch4    ,uco211  ,uco212  , &
                  uco213  ,uco221  ,uco222  ,uco223  ,uptype  , &
                  bn2o0   ,bn2o1   ,bch4    ,abplnk1 ,abplnk2 , &
                  abstot  ,absnxt  ,plh2ob  ,wb      , &
                  aer_mpp ,aer_trn_ttl)
   end if
!
! compute sums used in integrals (all longitude points)
!
! definition of bk1 & bk2 depends on finite differencing.  for
! trapezoidal rule bk1=bk2. trapezoidal rule applied for nonadjacent
! layers only.
!
! delt=t**4 in layer above current sigma level km.
! delt1=t**4 in layer below current sigma level km.
!
   do i=1,ncol
      delt(i) = tint4(i,pver) - tlayr4(i,pverp)
      delt1(i) = tlayr4(i,pverp) - tint4(i,pverp)
      s(i,pverp,pverp) = stebol*(delt1(i)*absnxt(i,pver,1) + delt (i)*absnxt(i,pver,4))
      s(i,pver,pverp)  = stebol*(delt (i)*absnxt(i,pver,2) + delt1(i)*absnxt(i,pver,3))
   end do
   do k=ntoplw,pver-1
      do i=1,ncol
         bk2(i) = (abstot(i,k,pver) + abstot(i,k,pverp))*0.5
         bk1(i) = bk2(i)
         s(i,k,pverp) = stebol*(bk2(i)*delt(i) + bk1(i)*delt1(i))
      end do
   end do
!
! all k, km>1
!
   do km=pver,ntoplw+1,-1
      do i=1,ncol
         delt(i)  = tint4(i,km-1) - tlayr4(i,km)
         delt1(i) = tlayr4(i,km) - tint4(i,km)
      end do
      do k=pverp,ntoplw,-1
         if (k == km) then
            do i=1,ncol
               bk2(i) = absnxt(i,km-1,4)
               bk1(i) = absnxt(i,km-1,1)
            end do
         else if (k == km-1) then
            do i=1,ncol
               bk2(i) = absnxt(i,km-1,2)
               bk1(i) = absnxt(i,km-1,3)
            end do
         else
            do i=1,ncol
               bk2(i) = (abstot(i,k,km-1) + abstot(i,k,km))*0.5
               bk1(i) = bk2(i)
            end do
         end if
         do i=1,ncol
            s(i,k,km) = s(i,k,km+1) + stebol*(bk2(i)*delt(i) + bk1(i)*delt1(i))
         end do
      end do
   end do
!
! computation of clear sky fluxes always set first level of fsul
!
   do i=1,ncol
      fsul(i,pverp) = lwupcgs(i)
   end do
!
! downward clear sky fluxes store intermediate quantities in down flux
! initialize fluxes to clear sky values.
!
   do i=1,ncol
      tmp(i) = fsul(i,pverp) - stebol*tint4(i,pverp)
      fsul(i,ntoplw) = fsul(i,pverp) - abstot(i,ntoplw,pverp)*tmp(i) + s(i,ntoplw,ntoplw+1)
      fsdl(i,ntoplw) = stebol*(tplnke(i)**4)*emstot(i,ntoplw)
   end do
!
! fsdl(i,pverp) assumes isothermal layer
!
   do k=ntoplw+1,pver
      do i=1,ncol
         fsul(i,k) = fsul(i,pverp) - abstot(i,k,pverp)*tmp(i) + s(i,k,k+1)
         fsdl(i,k) = stebol*(tplnke(i)**4)*emstot(i,k) - (s(i,k,ntoplw+1) - s(i,k,k+1))
      end do
   end do
!
! store the downward emission from level 1 = total gas emission * sigma
! t**4.  fsdl does not yet include all terms
!
   do i=1,ncol
      absbt(i) = stebol*(tplnke(i)**4)*emstot(i,pverp)
      fsdl(i,pverp) = absbt(i) - s(i,pverp,ntoplw+1)
   end do
!
!----------------------------------------------------------------------
! modifications for clouds -- max/random overlap assumption
!
! the column is divided into sets of adjacent layers, called regions,
!   in which the clouds are maximally overlapped.  the clouds are
!   randomly overlapped between different regions.  the number of
!   regions in a column is set by nmxrgn, and the range of pressures
!   included in each region is set by pmxrgn.  the max/random overlap
!   can be written in terms of the solutions of random overlap with
!   cloud amounts = 1.  the random overlap assumption is equivalent to
!   setting the flux boundary conditions (bcs) at the edges of each region
!   equal to the mean all-sky flux at those boundaries.  since the
!   emissivity array for propogating bcs is only computed for the
!   toa bc, the flux bcs elsewhere in the atmosphere have to be formulated
!   in terms of solutions to the random overlap equations.  this is done
!   by writing the flux bcs as the sum of a clear-sky flux and emission
!   from a cloud outside the region weighted by an emissivity.  this
!   emissivity is determined from the location of the cloud and the
!   flux bc.
!
! copy cloud amounts to buffer with extra layer (needed for overlap logic)
!
   cldp(:ncol,ntoplw:pver) = cld(:ncol,ntoplw:pver)
   cldp(:ncol,pverp) = 0.0
!
!
! select only those locations where there are no clouds
!    (maximum cloud fraction <= 1.e-3 treated as clear)
!    set all-sky fluxes to clear-sky values.
!
   maxcld(1:ncol) = maxval(cldp(1:ncol,ntoplw:pver),dim=2)

   npts = 0
   do i=1,ncol
      if (maxcld(i) < cldmin) then
         npts = npts + 1
         indx(npts) = i
      end if
   end do

   do ii = 1, npts
      i = indx(ii)
      do k = ntoplw, pverp
         fdl(i,k) = fsdl(i,k)
         ful(i,k) = fsul(i,k)
      end do
   end do
!
! select only those locations where there are clouds
!
   npts = 0
   do i=1,ncol
      if (maxcld(i) >= cldmin) then
         npts = npts + 1
         indx(npts) = i
      end if
   end do

!
! initialize all-sky fluxes. fdl(i,1) & ful(i,pverp) are boundary conditions
!
   do ii = 1, npts
      i = indx(ii)
      fdl(i,ntoplw) = fsdl(i,ntoplw)
      fdl(i,pverp)  = 0.0
      ful(i,ntoplw) = 0.0
      ful(i,pverp)  = fsul(i,pverp)
      do k = ntoplw+1, pver
         fdl(i,k) = 0.0
         ful(i,k) = 0.0
      end do
!
! initialize planck emission from layer boundaries
!
      do k = ntoplw, pver
         fclt4(i,k-1) = stebol*tint4(i,k)
         fclb4(i,k-1) = stebol*tint4(i,k+1)
      enddo
      fclb4(i,ntoplw-2) =  stebol*tint4(i,ntoplw)
      fclt4(i,pver)     = stebol*tint4(i,pverp)
!
! initialize indices for layers to be max-overlapped
!
      do irgn = 0, nmxrgn(i)
         kx2(i,irgn) = ntoplw-1
      end do
      nrgn(i) = 0
   end do

!----------------------------------------------------------------------
! index calculations for max overlap

   do ii = 1, npts
      ilon = indx(ii)

!
! outermost loop over regions (sets of adjacent layers) to be max overlapped
!
      do irgn = 1, nmxrgn(ilon)
!
! calculate min/max layer indices inside region.
!
         n = 0
         if (kx2(ilon,irgn-1) < pver) then
            nrgn(ilon) = irgn
            k1 = kx2(ilon,irgn-1)+1
            kx1(ilon,irgn) = k1
            kx2(ilon,irgn) = 0
            do k2 = pver, k1, -1
               if (pmid(ilon,k2) <= pmxrgn(ilon,irgn)) then
                  kx2(ilon,irgn) = k2
                  exit
               end if
            end do
!
! identify columns with clouds in the given region.
!
            do k = k1, k2
               if (cldp(ilon,k) >= cldmin) then
                  n = n+1
                  indxmx(n,irgn) = ilon
                  exit
               endif
            end do
         endif
         ncolmx(irgn) = n
!
! dummy value for handling clear-sky regions
!
!!$         indxmx(ncolmx(irgn)+1,irgn) = ncol+1
!
! outer loop over columns with clouds in the max-overlap region
!
         do iimx = 1, ncolmx(irgn)
            i = indxmx(iimx,irgn)
!
! sort cloud areas and corresponding level indices.
!
            n = 0
            do k = kx1(i,irgn),kx2(i,irgn)
               if (cldp(i,k) >= cldmin) then
                  n = n+1
                  ksort(n) = k
!
! we need indices for clouds in order of largest to smallest, so
!    sort 1-cld in ascending order
!
                  asort(n) = 1.0-cldp(i,k)
               end if
            end do
            nxs(i,irgn) = n
!
! if nxs(i,irgn) eq 1, no need to sort.
! if nxs(i,irgn) eq 2, sort by swapping if necessary
! if nxs(i,irgn) ge 3, sort using local sort routine
!
            if (nxs(i,irgn) == 2) then
               if (asort(2) < asort(1)) then
                  ktmp = ksort(1)
                  ksort(1) = ksort(2)
                  ksort(2) = ktmp

                  atmp = asort(1)
                  asort(1) = asort(2)
                  asort(2) = atmp
               endif
            else if (nxs(i,irgn) >= 3) then
               call sortarray(nxs(i,irgn),asort,ksort(1:))
            endif

            do l = 1, nxs(i,irgn)
               kxs(l,i,irgn) = ksort(l)
            end do
!
! end loop over longitude i for fluxes
!
         end do
!
! end loop over regions irgn for max-overlap
!
      end do
!
!----------------------------------------------------------------------
! downward fluxes:
! outermost loop over regions (sets of adjacent layers) to be max overlapped
!
      do irgn = 1, nmxrgn(ilon)
!
! compute clear-sky fluxes for regions without clouds
!
         iimx = 1
         if (ilon < indxmx(iimx,irgn) .and. irgn <= nrgn(ilon)) then
!
! calculate emissivity so that downward flux at upper boundary of region
!    can be cast in form of solution for downward flux from cloud above
!    that boundary.  then solutions for fluxes at other levels take form of
!    random overlap expressions.  try to locate "cloud" as close as possible
!    to toa such that the "cloud" pseudo-emissivity is between 0 and 1.
!
            k1 = kx1(ilon,irgn)
            do km1 = ntoplw-2, k1-2
               km4 = km1+3
               k2 = k1
               k3 = k2+1
               tmp(ilon) = s(ilon,k2,min(k3,pverp))*min(1,pverp2-k3)
               emx0 = (fdl(ilon,k1)-fsdl(ilon,k1))/ &
                      ((fclb4(ilon,km1)-s(ilon,k2,km4)+tmp(ilon))- fsdl(ilon,k1))
               if (emx0 >= 0.0 .and. emx0 <= 1.0) exit
            end do
            km1 = min(km1,k1-2)
            do k2 = kx1(ilon,irgn)+1, kx2(ilon,irgn)+1
               k3 = k2+1
               tmp(ilon) = s(ilon,k2,min(k3,pverp))*min(1,pverp2-k3)
               fdl(ilon,k2) = (1.0-emx0)*fsdl(ilon,k2) + &
                               emx0*(fclb4(ilon,km1)-s(ilon,k2,km4)+tmp(ilon))
            end do
         else if (ilon==indxmx(iimx,irgn) .and. iimx<=ncolmx(irgn)) then
            iimx = iimx+1
         end if
!
! outer loop over columns with clouds in the max-overlap region
!
         do iimx = 1, ncolmx(irgn)
            i = indxmx(iimx,irgn)

!
! calculate emissivity so that downward flux at upper boundary of region
!    can be cast in form of solution for downward flux from cloud above that
!    boundary.  then solutions for fluxes at other levels take form of
!    random overlap expressions.  try to locate "cloud" as close as possible
!    to toa such that the "cloud" pseudo-emissivity is between 0 and 1.
!
            k1 = kx1(i,irgn)
            do km1 = ntoplw-2,k1-2
               km4 = km1+3
               k2 = k1
               k3 = k2 + 1
               tmp(i) = s(i,k2,min(k3,pverp))*min(1,pverp2-k3)
               tmp2(i) = s(i,k2,min(km4,pverp))*min(1,pverp2-km4)
               emx0 = (fdl(i,k1)-fsdl(i,k1))/((fclb4(i,km1)-tmp2(i)+tmp(i))-fsdl(i,k1))
               if (emx0 >= 0.0 .and. emx0 <= 1.0) exit
            end do
            km1 = min(km1,k1-2)
            ksort(0) = km1 + 1
!
! loop to calculate fluxes at level k
!
            nxsk = 0
            do k = kx1(i,irgn), kx2(i,irgn)
!
! identify clouds (largest to smallest area) between kx1 and k
!    since nxsk will increase with increasing k up to nxs(i,irgn), once
!    nxsk == nxs(i,irgn) then use the list constructed for previous k
!
               if (nxsk < nxs(i,irgn)) then
                  nxsk = 0
                  do l = 1, nxs(i,irgn)
                     k1 = kxs(l,i,irgn)
                     if (k >= k1) then
                        nxsk = nxsk + 1
                        ksort(nxsk) = k1
                     endif
                  end do
               endif
!
! dummy value of index to insure computation of cloud amt is valid for l=nxsk+1
!
               ksort(nxsk+1) = pverp
!
! initialize iterated emissivity factors
!
               do l = 1, nxsk
                  emx(l) = emis(i,ksort(l))
               end do
!
! initialize iterated emissivity factor for bnd. condition at upper interface
!
               emx(0) = emx0
!
! initialize previous cloud amounts
!
               cld0 = 1.0
!
! indices for flux calculations
!
               k2 = k+1
               k3 = k2+1
               tmp(i) = s(i,k2,min(k3,pverp))*min(1,pverp2-k3)
!
! loop over number of cloud levels inside region (biggest to smallest cld area)
!
               do l = 1, nxsk+1
!
! calculate downward fluxes
!
                  cld1 = cldp(i,ksort(l))*min(1,nxsk+1-l)
                  if (cld0 /= cld1) then
                     fdl(i,k2) = fdl(i,k2)+(cld0-cld1)*fsdl(i,k2)
                     do l1 = 0, l - 1
                        km1 = ksort(l1)-1
                        km4 = km1+3
                        tmp2(i) = s(i,k2,min(km4,pverp))* min(1,pverp2-km4)
                        fdl(i,k2) = fdl(i,k2)+(cld0-cld1)*emx(l1)*(fclb4(i,km1)-tmp2(i)+tmp(i)- &
                                    fsdl(i,k2))
                     end do
                  endif
                  cld0 = cld1
!
! multiply emissivity factors by current cloud transmissivity
!
                  if (l <= nxsk) then
                     k1 = ksort(l)
                     trans = 1.0-emis(i,k1)
!
! ideally the upper bound on l1 would be l-1, but the sort routine
!    scrambles the order of layers with identical cloud amounts
!
                     do l1 = 0, nxsk
                        if (ksort(l1) < k1) then
                           emx(l1) = emx(l1)*trans
                        endif
                     end do
                  end if
!
! end loop over number l of cloud levels
!
               end do
!
! end loop over level k for fluxes
!
            end do
!
! end loop over longitude i for fluxes
!
         end do
!
! end loop over regions irgn for max-overlap
!
      end do

!
!----------------------------------------------------------------------
! upward fluxes:
! outermost loop over regions (sets of adjacent layers) to be max overlapped
!
      do irgn = nmxrgn(ilon), 1, -1
!
! compute clear-sky fluxes for regions without clouds
!
         iimx = 1
         if (ilon < indxmx(iimx,irgn) .and. irgn <= nrgn(ilon)) then
!
! calculate emissivity so that upward flux at lower boundary of region
!    can be cast in form of solution for upward flux from cloud below that
!    boundary.  then solutions for fluxes at other levels take form of
!    random overlap expressions.  try to locate "cloud" as close as possible
!    to surface such that the "cloud" pseudo-emissivity is between 0 and 1.
! include allowance for surface emissivity (both numerator and denominator
!    equal 1)
!
            k1 = kx2(ilon,irgn)+1
            if (k1 < pverp) then
               do km1 = pver-1,kx2(ilon,irgn),-1
                  km3 = km1+2
                  k2 = k1
                  k3 = k2+1
                  tmp(ilon) = s(ilon,k2,min(km3,pverp))* min(1,pverp2-km3)
                  emx0 = (ful(ilon,k1)-fsul(ilon,k1))/ &
                         ((fclt4(ilon,km1)+s(ilon,k2,k3)-tmp(ilon))- fsul(ilon,k1))
                  if (emx0 >= 0.0 .and. emx0 <= 1.0) exit
               end do
               km1 = max(km1,kx2(ilon,irgn))
            else
               km1 = k1-1
               km3 = km1+2
               emx0 = 1.0
            endif

            do k2 = kx1(ilon,irgn), kx2(ilon,irgn)
               k3 = k2+1
!
! if km3 == pver+2, one of the s integrals = 0 (integration limits both = p_s)
!
               tmp(ilon) = s(ilon,k2,min(km3,pverp))* min(1,pverp2-km3)
               ful(ilon,k2) =(1.0-emx0)*fsul(ilon,k2) + emx0* &
                             (fclt4(ilon,km1)+s(ilon,k2,k3)-tmp(ilon))
            end do
         else if (ilon==indxmx(iimx,irgn) .and. iimx<=ncolmx(irgn)) then
            iimx = iimx+1
         end if
!
! outer loop over columns with clouds in the max-overlap region
!
         do iimx = 1, ncolmx(irgn)
            i = indxmx(iimx,irgn)

!
! calculate emissivity so that upward flux at lower boundary of region
!    can be cast in form of solution for upward flux from cloud at that
!    boundary.  then solutions for fluxes at other levels take form of
!    random overlap expressions.  try to locate "cloud" as close as possible
!    to surface such that the "cloud" pseudo-emissivity is between 0 and 1.
! include allowance for surface emissivity (both numerator and denominator
!    equal 1)
!
            k1 = kx2(i,irgn)+1
            if (k1 < pverp) then
               do km1 = pver-1,kx2(i,irgn),-1
                  km3 = km1+2
                  k2 = k1
                  k3 = k2+1
                  tmp(i) = s(i,k2,min(km3,pverp))*min(1,pverp2-km3)
                  emx0 = (ful(i,k1)-fsul(i,k1))/((fclt4(i,km1)+s(i,k2,k3)-tmp(i))-fsul(i,k1))
                  if (emx0 >= 0.0 .and. emx0 <= 1.0) exit
               end do
               km1 = max(km1,kx2(i,irgn))
            else
               emx0 = 1.0
               km1 = k1-1
            endif
            ksort(0) = km1 + 1

!
! loop to calculate fluxes at level k
!
            nxsk = 0
            do k = kx2(i,irgn), kx1(i,irgn), -1
!
! identify clouds (largest to smallest area) between k and kx2
!    since nxsk will increase with decreasing k up to nxs(i,irgn), once
!    nxsk == nxs(i,irgn) then use the list constructed for previous k
!
               if (nxsk < nxs(i,irgn)) then
                  nxsk = 0
                  do l = 1, nxs(i,irgn)
                     k1 = kxs(l,i,irgn)
                     if (k <= k1) then
                        nxsk = nxsk + 1
                        ksort(nxsk) = k1
                     endif
                  end do
               endif
!
! dummy value of index to insure computation of cloud amt is valid for l=nxsk+1
!
               ksort(nxsk+1) = pverp
!
! initialize iterated emissivity factors
!
               do l = 1, nxsk
                  emx(l) = emis(i,ksort(l))
               end do
!
! initialize iterated emissivity factor for bnd. condition at lower interface
!
               emx(0) = emx0
!
! initialize previous cloud amounts
!
               cld0 = 1.0
!
! indices for flux calculations
!
               k2 = k
               k3 = k2+1
!
! loop over number of cloud levels inside region (biggest to smallest cld area)
!
               do l = 1, nxsk+1
!
! calculate upward fluxes
!
                  cld1 = cldp(i,ksort(l))*min(1,nxsk+1-l)
                  if (cld0 /= cld1) then
                     ful(i,k2) = ful(i,k2)+(cld0-cld1)*fsul(i,k2)
                     do l1 = 0, l - 1
                        km1 = ksort(l1)-1
                        km3 = km1+2
!
! if km3 == pver+2, one of the s integrals = 0 (integration limits both = p_s)
!
                        tmp(i) = s(i,k2,min(km3,pverp))* min(1,pverp2-km3)
                        ful(i,k2) = ful(i,k2)+(cld0-cld1)*emx(l1)* &
                                   (fclt4(i,km1)+s(i,k2,k3)-tmp(i)- fsul(i,k2))
                     end do
                  endif
                  cld0 = cld1
!
! multiply emissivity factors by current cloud transmissivity
!
                  if (l <= nxsk) then
                     k1 = ksort(l)
                     trans = 1.0-emis(i,k1)
!
! ideally the upper bound on l1 would be l-1, but the sort routine
!    scrambles the order of layers with identical cloud amounts
!
                     do l1 = 0, nxsk
                        if (ksort(l1) > k1) then
                           emx(l1) = emx(l1)*trans
                        endif
                     end do
                  end if
!
! end loop over number l of cloud levels
!
               end do
!
! end loop over level k for fluxes
!
            end do
!
! end loop over longitude i for fluxes
!
         end do
!
! end loop over regions irgn for max-overlap
!
      end do
!
! end outermost longitude loop
!
   end do
!
! end cloud modification loops
!
!----------------------------------------------------------------------
! all longitudes: store history tape quantities
!
   do i=1,ncol
      flwds(i) = fdl (i,pverp )
      flns(i)  = ful (i,pverp ) - fdl (i,pverp )
      flnsc(i) = fsul(i,pverp ) - fsdl(i,pverp )
      flnt(i)  = ful (i,ntoplw) - fdl (i,ntoplw)
      flntc(i) = fsul(i,ntoplw) - fsdl(i,ntoplw)
      flut(i)  = ful (i,ntoplw)
      flutc(i) = fsul(i,ntoplw)
   end do
!
! computation of longwave heating (j/kg/s)
!
   do k=ntoplw,pver
      do i=1,ncol
         qrl(i,k) = (ful(i,k) - fdl(i,k) - ful(i,k+1) + fdl(i,k+1))* &
                     1.e-4*gravit/((pint(i,k) - pint(i,k+1)))
      end do
   end do
! return 0 above solution domain
   if ( ntoplw > 1 )then
      qrl(:ncol,:ntoplw-1) = 0.
   end if

! added downward/upward total and clear sky fluxes
!
   do k=ntoplw,pverp
      do i=1,ncol
        flup(i,k)  = ful(i,k)
        flupc(i,k) = fsul(i,k)
        fldn(i,k)  = fdl(i,k)
        fldnc(i,k) = fsdl(i,k)
      end do
   end do
! return 0 above solution domain
   if ( ntoplw > 1 )then
      flup(:ncol,:ntoplw-1) = 0.
      flupc(:ncol,:ntoplw-1) = 0.
      fldn(:ncol,:ntoplw-1) = 0.
      fldnc(:ncol,:ntoplw-1) = 0.
   end if
!
   return
end subroutine radclwmx

subroutine radcswmx(jj, lchnk   ,ncol    ,pcols, pver, pverp,         &
                    pint    ,pmid    ,h2ommr  ,rh      ,o3mmr   , &
                    aermmr  ,cld     ,cicewp  ,cliqwp  ,rel     , &
!                   rei     ,eccf    ,coszrs  ,scon    ,solin   ,solcon,  &
                    rei     ,tauxcl  ,tauxci  ,eccf    ,coszrs  ,scon    ,solin   ,solcon,  &
                    asdir   ,asdif   ,aldir   ,aldif   ,nmxrgn  , &
                    pmxrgn  ,qrs     ,fsnt    ,fsntc   ,fsntoa  , &
                    fsntoac ,fsnirtoa,fsnrtoac,fsnrtoaq,fsns    , &
                    fsnsc   ,fsdsc   ,fsds    ,sols    ,soll    , &
                    solsd   ,solld   ,frc_day ,                   &
                    fsup    ,fsupc   ,fsdn    ,fsdnc   ,          &
                    aertau  ,aerssa  ,aerasm  ,aerfwd             )
!-----------------------------------------------------------------------
! 
! purpose: 
! solar radiation code
! 
! method: 
! basic method is delta-eddington as described in:
! 
! briegleb, bruce p., 1992: delta-eddington
! approximation for solar radiation in the ncar community climate model,
! journal of geophysical research, vol 97, d7, pp7603-7612).
! 
! five changes to the basic method described above are:
! (1) addition of sulfate aerosols (kiehl and briegleb, 1993)
! (2) the distinction between liquid and ice particle clouds 
! (kiehl et al, 1996);
! (3) provision for calculating toa fluxes with spectral response to
! match nimbus-7 visible/near-ir radiometers (collins, 1998);
! (4) max-random overlap (collins, 2001)
! (5) the near-ir absorption by h2o was updated in 2003 by collins, 
!     lee-taylor, and edwards for consistency with the new line data in
!     hitran 2000 and the h2o continuum version ckd 2.4.  modifications
!     were optimized by reducing rms errors in heating rates relative
!     to a series of benchmark calculations for the 5 standard afgl 
!     atmospheres.  the benchmarks were performed using disort2 combined
!     with genln3.  the near-ir scattering optical depths for rayleigh
!     scattering were also adjusted, as well as the correction for
!     stratospheric heating by h2o.
!
! the treatment of maximum-random overlap is described in the
! comment block "index calculations for max overlap".
! 
! divides solar spectrum into 19 intervals from 0.2-5.0 micro-meters.
! solar flux fractions specified for each interval. allows for
! seasonally and diurnally varying solar input.  includes molecular,
! cloud, aerosol, and surface scattering, along with h2o,o3,co2,o2,cloud, 
! and surface absorption. computes delta-eddington reflections and
! transmissions assuming homogeneously mixed layers. adds the layers 
! assuming scattering between layers to be isotropic, and distinguishes 
! direct solar beam from scattered radiation.
! 
! longitude loops are broken into 1 or 2 sections, so that only daylight
! (i.e. coszrs > 0) computations are done.
! 
! note that an extra layer above the model top layer is added.
! 
! cgs units are used.
! 
! special diagnostic calculation of the clear sky surface and total column
! absorbed flux is also done for cloud forcing diagnostics.
! 
!-----------------------------------------------------------------------
!  use shr_kind_mod, only: r8 => shr_kind_r8
!  use ppgrid
!  use ghg_surfvals, only: co2mmr
!  use prescribed_aerosols, only: idxbg, idxsul, idxsslt, idxocpho, idxbcpho, idxocphi, idxbcphi, &
!    idxdustfirst, numdust, idxvolc, naer_all
!  use aer_optics, only: nrh, ndstsz, ksul, wsul, gsul, &
!    ksslt, wsslt, gsslt, kcphil, wcphil, gcphil, kcphob, wcphob, gcphob, &
!    kcb, wcb, gcb, kdst, wdst, gdst, kbg, wbg, gbg, kvolc, wvolc, gvolc
!  use abortutils, only: endrun

   implicit none

   integer nspint            ! num of spctrl intervals across solar spectrum
   integer naer_groups       ! num of aerosol groups for optical diagnostics

   parameter ( nspint = 19 )
   parameter ( naer_groups = 7 )    ! current groupings are sul, sslt, all carbons, all dust, and all aerosols
!-----------------------constants for new band (640-700 nm)-------------

!-------------parameters for accelerating max-random solution-------------
! 
! the solution time scales like prod(j:1->n) (1 + n_j) where 
! n   = number of max-overlap regions (nmxrgn)
! n_j = number of unique cloud amounts in region j
! 
! therefore the solution cost can be reduced by decreasing n_j.
! cldmin reduces n_j by treating cloud amounts < cldmin as clear sky.
! cldeps reduces n_j by treating cloud amounts identical to log(1/cldeps)
! decimal places as identical
! 
! areamin reduces the cost by dropping configurations that occupy
! a surface area < areamin of the model grid box.  the surface area
! for a configuration c(j,k_j), where j is the region number and k_j is the
! index for a unique cloud amount (in descending order from biggest to
! smallest clouds) in region j, is
! 
! a = prod(j:1->n) [c(j,k_j) - c(j,k_j+1)]
! 
! where c(j,0) = 1.0 and c(j,n_j+1) = 0.0.
! 
! nconfgmax reduces the cost and improves load balancing by setting an upper
! bound on the number of cloud configurations in the solution.  if the number
! of configurations exceeds nconfgmax, the nconfgmax configurations with the
! largest area are retained, and the fluxes are normalized by the total area
! of these nconfgmax configurations.  for the current max/random overlap 
! assumption (see subroutine cldovrlap), 30 levels, and cloud-amount 
! parameterization, the mean and rms number of configurations are 
! both roughly 5.  nconfgmax has been set to the mean+2*rms number, or 15.
! 
! minimum cloud amount (as a fraction of the grid-box area) to 
! distinguish from clear sky
! 
   real(r8) cldmin
   parameter (cldmin = 1.0e-80_r8)
! 
! minimimum horizontal area (as a fraction of the grid-box area) to retain 
! for a unique cloud configuration in the max-random solution
! 
   real(r8) areamin
   parameter (areamin = 0.01_r8)
! 
! decimal precision of cloud amount (0 -> preserve full resolution;
! 10^-n -> preserve n digits of cloud amount)
! 
   real(r8) cldeps
   parameter (cldeps = 0.0_r8)
! 
! maximum number of configurations to include in solution
! 
   integer nconfgmax
   parameter (nconfgmax = 15)
!------------------------------commons----------------------------------
! 
! input arguments
! 
   integer, intent(in) :: lchnk,jj             ! chunk identifier
   integer, intent(in) :: pcols, pver, pverp
   integer, intent(in) :: ncol              ! number of atmospheric columns

   real(r8), intent(in) :: pmid(pcols,pver) ! level pressure
   real(r8), intent(in) :: pint(pcols,pverp) ! interface pressure
   real(r8), intent(in) :: h2ommr(pcols,pver) ! specific humidity (h2o mass mix ratio)
   real(r8), intent(in) :: o3mmr(pcols,pver) ! ozone mass mixing ratio
   real(r8), intent(in) :: aermmr(pcols,pver,naer_all) ! aerosol mass mixing ratio
   real(r8), intent(in) :: rh(pcols,pver)   ! relative humidity (fraction)
! 
   real(r8), intent(in) :: cld(pcols,pver)  ! fractional cloud cover
   real(r8), intent(in) :: cicewp(pcols,pver) ! in-cloud cloud ice water path
   real(r8), intent(in) :: cliqwp(pcols,pver) ! in-cloud cloud liquid water path
   real(r8), intent(in) :: rel(pcols,pver)  ! liquid effective drop size (microns)
   real(r8), intent(in) :: rei(pcols,pver)  ! ice effective drop size (microns)
! 
   real(r8), intent(in) :: eccf             ! eccentricity factor (1./earth-sun dist^2)
   real, intent(in) :: solcon           ! solar constant with eccentricity factor
   real(r8), intent(in) :: coszrs(pcols)    ! cosine solar zenith angle
   real(r8), intent(in) :: asdir(pcols)     ! 0.2-0.7 micro-meter srfc alb: direct rad
   real(r8), intent(in) :: aldir(pcols)     ! 0.7-5.0 micro-meter srfc alb: direct rad
   real(r8), intent(in) :: asdif(pcols)     ! 0.2-0.7 micro-meter srfc alb: diffuse rad
   real(r8), intent(in) :: aldif(pcols)     ! 0.7-5.0 micro-meter srfc alb: diffuse rad

   real(r8), intent(in) :: scon             ! solar constant 
! 
! in/out arguments
! 
   real(r8), intent(inout) :: pmxrgn(pcols,pverp) ! maximum values of pressure for each
!                                                 !    maximally overlapped region. 
!                                                 !    0->pmxrgn(i,1) is range of pressure for
!                                                 !    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
!                                                 !    2nd region, etc
   integer, intent(inout) ::  nmxrgn(pcols)    ! number of maximally overlapped regions
! 
! output arguments
! 

   real(r8), intent(out) :: solin(pcols)     ! incident solar flux
   real(r8), intent(out) :: qrs(pcols,pver)  ! solar heating rate
   real(r8), intent(out) :: fsns(pcols)      ! surface absorbed solar flux
   real(r8), intent(out) :: fsnt(pcols)      ! total column absorbed solar flux
   real(r8), intent(out) :: fsntoa(pcols)    ! net solar flux at toa
   real(r8), intent(out) :: fsds(pcols)      ! flux shortwave downwelling surface
! 
   real(r8), intent(out) :: fsnsc(pcols)     ! clear sky surface absorbed solar flux
   real(r8), intent(out) :: fsdsc(pcols)     ! clear sky surface downwelling solar flux
   real(r8), intent(out) :: fsntc(pcols)     ! clear sky total column absorbed solar flx
   real(r8), intent(out) :: fsntoac(pcols)   ! clear sky net solar flx at toa
   real(r8), intent(out) :: sols(pcols)      ! direct solar rad on surface (< 0.7)
   real(r8), intent(out) :: soll(pcols)      ! direct solar rad on surface (>= 0.7)
   real(r8), intent(out) :: solsd(pcols)     ! diffuse solar rad on surface (< 0.7)
   real(r8), intent(out) :: solld(pcols)     ! diffuse solar rad on surface (>= 0.7)
   real(r8), intent(out) :: fsnirtoa(pcols)  ! near-ir flux absorbed at toa
   real(r8), intent(out) :: fsnrtoac(pcols)  ! clear sky near-ir flux absorbed at toa
   real(r8), intent(out) :: fsnrtoaq(pcols)  ! net near-ir flux at toa >= 0.7 microns
   real(r8), intent(out) :: tauxcl(pcols,0:pver) ! water cloud extinction optical depth
   real(r8), intent(out) :: tauxci(pcols,0:pver) ! ice cloud extinction optical depth

! added downward/upward total and clear sky fluxes
   real(r8), intent(out) :: fsup(pcols,pverp)      ! total sky upward solar flux (spectrally summed)
   real(r8), intent(out) :: fsupc(pcols,pverp)     ! clear sky upward solar flux (spectrally summed)
   real(r8), intent(out) :: fsdn(pcols,pverp)      ! total sky downward solar flux (spectrally summed)
   real(r8), intent(out) :: fsdnc(pcols,pverp)     ! clear sky downward solar flux (spectrally summed)
!
   real(r8) , intent(out) :: frc_day(pcols) ! = 1 for daylight, =0 for night columns
   real(r8) :: aertau(pcols,nspint,naer_groups) ! aerosol column optical depth
   real(r8) :: aerssa(pcols,nspint,naer_groups) ! aerosol column averaged single scattering albedo
   real(r8) :: aerasm(pcols,nspint,naer_groups) ! aerosol column averaged asymmetry parameter
   real(r8) :: aerfwd(pcols,nspint,naer_groups) ! aerosol column averaged forward scattering
!  real(r8), intent(out) :: aertau(pcols,nspint,naer_groups) ! aerosol column optical depth
!  real(r8), intent(out) :: aerssa(pcols,nspint,naer_groups) ! aerosol column averaged single scattering albedo
!  real(r8), intent(out) :: aerasm(pcols,nspint,naer_groups) ! aerosol column averaged asymmetry parameter
!  real(r8), intent(out) :: aerfwd(pcols,nspint,naer_groups) ! aerosol column averaged forward scattering
! 
!---------------------------local variables-----------------------------
! 
! max/random overlap variables
! 
   real(r8) asort(pverp)     ! 1 - cloud amounts to be sorted for max ovrlp.
   real(r8) atmp             ! temporary storage for sort when nxs = 2
   real(r8) cld0             ! 1 - (cld amt) used to make wstr, cstr, nstr
   real(r8) totwgt           ! total of xwgts = total fractional area of 
!   grid-box covered by cloud configurations
!   included in solution to fluxes

   real(r8) wgtv(nconfgmax)  ! weights for fluxes
!   1st index is configuration number
   real(r8) wstr(pverp,pverp) ! area weighting factors for streams
!   1st index is for stream #, 
!   2nd index is for region #

   real(r8) xexpt            ! solar direct beam trans. for layer above
   real(r8) xrdnd            ! diffuse reflectivity for layer above
   real(r8) xrupd            ! diffuse reflectivity for layer below
   real(r8) xrups            ! direct-beam reflectivity for layer below
   real(r8) xtdnt            ! total trans for layers above

   real(r8) xwgt             ! product of cloud amounts

   real(r8) yexpt            ! solar direct beam trans. for layer above
   real(r8) yrdnd            ! diffuse reflectivity for layer above
   real(r8) yrupd            ! diffuse reflectivity for layer below
   real(r8) ytdnd            ! dif-beam transmission for layers above
   real(r8) ytupd            ! dif-beam transmission for layers below

   real(r8) zexpt            ! solar direct beam trans. for layer above
   real(r8) zrdnd            ! diffuse reflectivity for layer above
   real(r8) zrupd            ! diffuse reflectivity for layer below
   real(r8) zrups            ! direct-beam reflectivity for layer below
   real(r8) ztdnt            ! total trans for layers above

   logical new_term          ! flag for configurations to include in fluxes
   logical region_found      ! flag for identifying regions

   integer ccon(0:pverp,nconfgmax)                                
! flags for presence of clouds
!   1st index is for level # (including 
!    layer above top of model and at surface)
!   2nd index is for configuration #
   integer cstr(0:pverp,pverp)                                
! flags for presence of clouds
!   1st index is for level # (including 
!    layer above top of model and at surface)
!   2nd index is for stream #
   integer icond(0:pverp,nconfgmax)
! indices for copying rad. properties from
!     one identical downward cld config.
!     to another in adding method (step 2)
!   1st index is for interface # (including 
!     layer above top of model and at surface)
!   2nd index is for configuration # range
   integer iconu(0:pverp,nconfgmax)
! indices for copying rad. properties from
!     one identical upward configuration
!     to another in adding method (step 2)
!   1st index is for interface # (including 
!     layer above top of model and at surface)
!   2nd index is for configuration # range
   integer iconfig           ! counter for random-ovrlap configurations
   integer irgn              ! index for max-overlap regions
   integer is0               ! lower end of stream index range
   integer is1               ! upper end of stream index range
   integer isn               ! stream index
   integer istr(pverp+1)     ! index for stream #s during flux calculation
   integer istrtd(0:pverp,0:nconfgmax+1)
! indices into icond 
!   1st index is for interface # (including 
!     layer above top of model and at surface)
!   2nd index is for configuration # range
   integer istrtu(0:pverp,0:nconfgmax+1)
! indices into iconu 
!   1st index is for interface # (including 
!     layer above top of model and at surface)
!   2nd index is for configuration # range
   integer j                 ! configuration index
   integer k1                ! level index
   integer k2                ! level index
   integer ksort(pverp)      ! level indices of cloud amounts to be sorted
   integer ktmp              ! temporary storage for sort when nxs = 2
   integer kx1(0:pverp)      ! level index for top of max-overlap region
   integer kx2(0:pverp)      ! level index for bottom of max-overlap region
   integer l                 ! index 
   integer l0                ! index
   integer mrgn              ! counter for nrgn
   integer mstr              ! counter for nstr
   integer n0                ! number of configurations with ccon(k,:)==0
   integer n1                ! number of configurations with ccon(k,:)==1
   integer nconfig           ! number of random-ovrlap configurations
   integer nconfigm          ! value of config before testing for areamin,
!    nconfgmax
   integer npasses           ! number of passes over the indexing loop
   integer nrgn              ! number of max overlap regions at current 
!    longitude
   integer nstr(pverp)       ! number of unique cloud configurations
!   ("streams") in a max-overlapped region
!   1st index is for region #
   integer nuniq             ! # of unique cloud configurations
   integer nuniqd(0:pverp)   ! # of unique cloud configurations: toa 
!   to level k
   integer nuniqu(0:pverp)   ! # of unique cloud configurations: surface
!   to level k 
   integer nxs               ! number of cloudy layers between k1 and k2 
   integer ptr0(nconfgmax)   ! indices of configurations with ccon(k,:)==0
   integer ptr1(nconfgmax)   ! indices of configurations with ccon(k,:)==1
   integer ptrc(nconfgmax)   ! pointer for configurations sorted by wgtv
!  integer findvalue         ! function for finding kth smallest element
!   in a vector
!  external findvalue

! 
! other
! 
   integer ns                ! spectral loop index
   integer i                 ! longitude loop index
   integer k                 ! level loop index
   integer km1               ! k - 1
   integer kp1               ! k + 1
   integer n                 ! loop index for daylight
   integer ndayc             ! number of daylight columns
   integer idayc(pcols)      ! daytime column indices
   integer indxsl            ! index for cloud particle properties
   integer ksz               ! dust size bin index
   integer krh               ! relative humidity bin index
   integer kaer              ! aerosol group index
   real(r8) wrh              ! weight for linear interpolation between lut points
   real(r8) :: rhtrunc       ! rh, truncated for the purposes of extrapolating
                             ! aerosol optical properties 
   real(r8) albdir(pcols,nspint) ! current spc intrvl srf alb to direct rad
   real(r8) albdif(pcols,nspint) ! current spc intrvl srf alb to diffuse rad
! 
   real(r8) wgtint           ! weight for specific spectral interval

! 
! diagnostic and accumulation arrays; note that sfltot, fswup, and
! fswdn are not used in the computation,but are retained for future use.
! 
   real(r8) solflx           ! solar flux in current interval
   real(r8) sfltot           ! spectrally summed total solar flux
   real(r8) totfld(0:pver)   ! spectrally summed flux divergence
   real(r8) fswup(0:pverp)   ! spectrally summed up flux
   real(r8) fswdn(0:pverp)   ! spectrally summed down flux
   real(r8) fswupc(0:pverp)  ! spectrally summed up clear sky flux
   real(r8) fswdnc(0:pverp)  ! spectrally summed down clear sky flux
! 
! cloud radiative property arrays
! 
!  real(r8) tauxcl(pcols,0:pver) ! water cloud extinction optical depth
!  real(r8) tauxci(pcols,0:pver) ! ice cloud extinction optical depth
   real(r8) wcl(pcols,0:pver) ! liquid cloud single scattering albedo
   real(r8) gcl(pcols,0:pver) ! liquid cloud asymmetry parameter
   real(r8) fcl(pcols,0:pver) ! liquid cloud forward scattered fraction
   real(r8) wci(pcols,0:pver) ! ice cloud single scattering albedo
   real(r8) gci(pcols,0:pver) ! ice cloud asymmetry parameter
   real(r8) fci(pcols,0:pver) ! ice cloud forward scattered fraction
!
! aerosol mass paths by species
!
  real(r8) usul(pcols,pver)   ! sulfate (so4)
  real(r8) ubg(pcols,pver)    ! background aerosol
  real(r8) usslt(pcols,pver)  ! sea-salt (sslt)
  real(r8) ucphil(pcols,pver) ! hydrophilic organic carbon (ocphi)
  real(r8) ucphob(pcols,pver) ! hydrophobic organic carbon (ocpho)
  real(r8) ucb(pcols,pver)    ! black carbon (bcphi + bcpho)
  real(r8) uvolc(pcols,pver) ! volcanic mass
  real(r8) udst(ndstsz,pcols,pver) ! dust

!
! local variables used for the external mixing of aerosol species
!
  real(r8) tau_sul             ! optical depth, sulfate
  real(r8) tau_bg              ! optical depth, background aerosol
  real(r8) tau_sslt            ! optical depth, sea-salt
  real(r8) tau_cphil           ! optical depth, hydrophilic carbon
  real(r8) tau_cphob           ! optical depth, hydrophobic carbon
  real(r8) tau_cb              ! optical depth, black carbon
  real(r8) tau_volc            ! optical depth, volcanic
  real(r8) tau_dst(ndstsz)     ! optical depth, dust, by size category
  real(r8) tau_dst_tot         ! optical depth, total dust
  real(r8) tau_tot             ! optical depth, total aerosol

  real(r8) tau_w_sul           ! optical depth * single scattering albedo, sulfate
  real(r8) tau_w_bg            ! optical depth * single scattering albedo, background aerosol
  real(r8) tau_w_sslt          ! optical depth * single scattering albedo, sea-salt
  real(r8) tau_w_cphil         ! optical depth * single scattering albedo, hydrophilic carbon
  real(r8) tau_w_cphob         ! optical depth * single scattering albedo, hydrophobic carbon
  real(r8) tau_w_cb            ! optical depth * single scattering albedo, black carbon
  real(r8) tau_w_volc          ! optical depth * single scattering albedo, volcanic
  real(r8) tau_w_dst(ndstsz)   ! optical depth * single scattering albedo, dust, by size
  real(r8) tau_w_dst_tot       ! optical depth * single scattering albedo, total dust
  real(r8) tau_w_tot           ! optical depth * single scattering albedo, total aerosol

  real(r8) tau_w_g_sul         ! optical depth * single scattering albedo * asymmetry parameter, sulfate
  real(r8) tau_w_g_bg          ! optical depth * single scattering albedo * asymmetry parameter, background aerosol
  real(r8) tau_w_g_sslt        ! optical depth * single scattering albedo * asymmetry parameter, sea-salt
  real(r8) tau_w_g_cphil       ! optical depth * single scattering albedo * asymmetry parameter, hydrophilic carbon
  real(r8) tau_w_g_cphob       ! optical depth * single scattering albedo * asymmetry parameter, hydrophobic carbon
  real(r8) tau_w_g_cb          ! optical depth * single scattering albedo * asymmetry parameter, black carbon
  real(r8) tau_w_g_volc        ! optical depth * single scattering albedo * asymmetry parameter, volcanic
  real(r8) tau_w_g_dst(ndstsz) ! optical depth * single scattering albedo * asymmetry parameter, dust, by size
  real(r8) tau_w_g_dst_tot     ! optical depth * single scattering albedo * asymmetry parameter, total dust
  real(r8) tau_w_g_tot         ! optical depth * single scattering albedo * asymmetry parameter, total aerosol

  real(r8) f_sul               ! forward scattering fraction, sulfate
  real(r8) f_bg                ! forward scattering fraction, background aerosol
  real(r8) f_sslt              ! forward scattering fraction, sea-salt
  real(r8) f_cphil             ! forward scattering fraction, hydrophilic carbon
  real(r8) f_cphob             ! forward scattering fraction, hydrophobic carbon
  real(r8) f_cb                ! forward scattering fraction, black carbon
  real(r8) f_volc              ! forward scattering fraction, volcanic
  real(r8) f_dst(ndstsz)       ! forward scattering fraction, dust, by size
  real(r8) f_dst_tot           ! forward scattering fraction, total dust
  real(r8) f_tot               ! forward scattering fraction, total aerosol

  real(r8) tau_w_f_sul         ! optical depth * forward scattering fraction * single scattering albedo, sulfate
  real(r8) tau_w_f_bg          ! optical depth * forward scattering fraction * single scattering albedo, background
  real(r8) tau_w_f_sslt        ! optical depth * forward scattering fraction * single scattering albedo, sea-salt
  real(r8) tau_w_f_cphil       ! optical depth * forward scattering fraction * single scattering albedo, hydrophilic c
  real(r8) tau_w_f_cphob       ! optical depth * forward scattering fraction * single scattering albedo, hydrophobic c
  real(r8) tau_w_f_cb          ! optical depth * forward scattering fraction * single scattering albedo, black c
  real(r8) tau_w_f_volc        ! optical depth * forward scattering fraction * single scattering albedo, volcanic
  real(r8) tau_w_f_dst(ndstsz) ! optical depth * forward scattering fraction * single scattering albedo, dust, by size
  real(r8) tau_w_f_dst_tot     ! optical depth * forward scattering fraction * single scattering albedo, total dust
  real(r8) tau_w_f_tot         ! optical depth * forward scattering fraction * single scattering albedo, total aerosol
  real(r8) w_dst_tot           ! single scattering albedo, total dust
  real(r8) w_tot               ! single scattering albedo, total aerosol
  real(r8) g_dst_tot           ! asymmetry parameter, total dust
  real(r8) g_tot               ! asymmetry parameter, total aerosol
  real(r8) ksuli               ! specific extinction interpolated between rh look-up-table points, sulfate
  real(r8) ksslti              ! specific extinction interpolated between rh look-up-table points, sea-salt
  real(r8) kcphili             ! specific extinction interpolated between rh look-up-table points, hydrophilic carbon
  real(r8) wsuli               ! single scattering albedo interpolated between rh look-up-table points, sulfate
  real(r8) wsslti              ! single scattering albedo interpolated between rh look-up-table points, sea-salt
  real(r8) wcphili             ! single scattering albedo interpolated between rh look-up-table points, hydrophilic carbon
  real(r8) gsuli               ! asymmetry parameter interpolated between rh look-up-table points, sulfate
  real(r8) gsslti              ! asymmetry parameter interpolated between rh look-up-table points, sea-salt
  real(r8) gcphili             ! asymmetry parameter interpolated between rh look-up-table points, hydrophilic carbon
! 
! aerosol radiative property arrays
! 
   real(r8) tauxar(pcols,0:pver) ! aerosol extinction optical depth
   real(r8) wa(pcols,0:pver) ! aerosol single scattering albedo
   real(r8) ga(pcols,0:pver) ! aerosol assymetry parameter
   real(r8) fa(pcols,0:pver) ! aerosol forward scattered fraction

! 
! various arrays and other constants:
! 
   real(r8) pflx(pcols,0:pverp) ! interface press, including extra layer
   real(r8) zenfac(pcols)    ! square root of cos solar zenith angle
   real(r8) sqrco2           ! square root of the co2 mass mixg ratio
   real(r8) tmp1             ! temporary constant array
   real(r8) tmp2             ! temporary constant array
   real(r8) pdel             ! pressure difference across layer
   real(r8) path             ! mass path of layer
   real(r8) ptop             ! lower interface pressure of extra layer
   real(r8) ptho2            ! used to compute mass path of o2
   real(r8) ptho3            ! used to compute mass path of o3
   real(r8) pthco2           ! used to compute mass path of co2
   real(r8) pthh2o           ! used to compute mass path of h2o
   real(r8) h2ostr           ! inverse sq. root h2o mass mixing ratio
   real(r8) wavmid(nspint)   ! spectral interval middle wavelength
   real(r8) trayoslp         ! rayleigh optical depth/standard pressure
   real(r8) tmp1l            ! temporary constant array
   real(r8) tmp2l            ! temporary constant array
   real(r8) tmp3l            ! temporary constant array
   real(r8) tmp1i            ! temporary constant array
   real(r8) tmp2i            ! temporary constant array
   real(r8) tmp3i            ! temporary constant array
   real(r8) rdenom           ! multiple scattering term
   real(r8) rdirexp          ! layer direct ref times exp transmission
   real(r8) tdnmexp          ! total transmission - exp transmission
   real(r8) psf(nspint)      ! frac of solar flux in spect interval
! 
! layer absorber amounts; note that 0 refers to the extra layer added
! above the top model layer
! 
   real(r8) uh2o(pcols,0:pver) ! layer absorber amount of h2o
   real(r8) uo3(pcols,0:pver) ! layer absorber amount of  o3
   real(r8) uco2(pcols,0:pver) ! layer absorber amount of co2
   real(r8) uo2(pcols,0:pver) ! layer absorber amount of  o2
   real(r8) uaer(pcols,0:pver) ! layer aerosol amount 
! 
! total column absorber amounts:
! 
   real(r8) uth2o(pcols)     ! total column  absorber amount of  h2o
   real(r8) uto3(pcols)      ! total column  absorber amount of  o3
   real(r8) utco2(pcols)     ! total column  absorber amount of  co2
   real(r8) uto2(pcols)      ! total column  absorber amount of  o2
! 
! these arrays are defined for pver model layers; 0 refers to the extra
! layer on top:
! 
   real(r8) rdir(nspint,pcols,0:pver) ! layer reflectivity to direct rad
   real(r8) rdif(nspint,pcols,0:pver) ! layer reflectivity to diffuse rad
   real(r8) tdir(nspint,pcols,0:pver) ! layer transmission to direct rad
   real(r8) tdif(nspint,pcols,0:pver) ! layer transmission to diffuse rad
   real(r8) explay(nspint,pcols,0:pver) ! solar beam exp trans. for layer

   real(r8) rdirc(nspint,pcols,0:pver) ! clear layer reflec. to direct rad
   real(r8) rdifc(nspint,pcols,0:pver) ! clear layer reflec. to diffuse rad
   real(r8) tdirc(nspint,pcols,0:pver) ! clear layer trans. to direct rad
   real(r8) tdifc(nspint,pcols,0:pver) ! clear layer trans. to diffuse rad
   real(r8) explayc(nspint,pcols,0:pver) ! solar beam exp trans. clear layer

   real(r8) flxdiv           ! flux divergence for layer
! 
! 
! radiative properties:
! 
! there are 1 classes of properties:
! (1. all-sky bulk properties
! (2. clear-sky properties
! 
! the first set of properties are generated during step 2 of the solution.
! 
! these arrays are defined at model interfaces; in 1st index (for level #),
! 0 is the top of the extra layer above the model top, and
! pverp is the earth surface.  2nd index is for cloud configuration
! defined over a whole column.
! 
   real(r8) exptdn(0:pverp,nconfgmax) ! sol. beam trans from layers above
   real(r8) rdndif(0:pverp,nconfgmax) ! ref to dif rad for layers above
   real(r8) rupdif(0:pverp,nconfgmax) ! ref to dif rad for layers below
   real(r8) rupdir(0:pverp,nconfgmax) ! ref to dir rad for layers below
   real(r8) tdntot(0:pverp,nconfgmax) ! total trans for layers above
! 
! bulk properties used during the clear-sky calculation.
! 
   real(r8) exptdnc(0:pverp) ! clr: sol. beam trans from layers above
   real(r8) rdndifc(0:pverp) ! clr: ref to dif rad for layers above
   real(r8) rupdifc(0:pverp) ! clr: ref to dif rad for layers below
   real(r8) rupdirc(0:pverp) ! clr: ref to dir rad for layers below
   real(r8) tdntotc(0:pverp) ! clr: total trans for layers above

   real(r8) fluxup(0:pverp)  ! up   flux at model interface
   real(r8) fluxdn(0:pverp)  ! down flux at model interface
   real(r8) wexptdn          ! direct solar beam trans. to surface

! moved to here from the module storage above, because these have to be thread-private.  jm 20100217
   real(r8) abarli           ! a coefficient for current spectral band
   real(r8) bbarli           ! b coefficient for current spectral band
   real(r8) cbarli           ! c coefficient for current spectral band
   real(r8) dbarli           ! d coefficient for current spectral band
   real(r8) ebarli           ! e coefficient for current spectral band
   real(r8) fbarli           ! f coefficient for current spectral band

   real(r8) abarii           ! a coefficient for current spectral band
   real(r8) bbarii           ! b coefficient for current spectral band
   real(r8) cbarii           ! c coefficient for current spectral band
   real(r8) dbarii           ! d coefficient for current spectral band
   real(r8) ebarii           ! e coefficient for current spectral band
   real(r8) fbarii           ! f coefficient for current spectral band
! jm 20100217

! 
!-----------------------------------------------------------------------
! start of calculation
!-----------------------------------------------------------------------
! 
!  write (6, '(a, x, i3)') 'radcswmx : chunk identifier', lchnk

   do i=1, ncol
! 
! initialize output fields:
! 
      fsds(i)     = 0.0_r8

      fsnirtoa(i) = 0.0_r8
      fsnrtoac(i) = 0.0_r8
      fsnrtoaq(i) = 0.0_r8

      fsns(i)     = 0.0_r8
      fsnsc(i)    = 0.0_r8
      fsdsc(i)    = 0.0_r8

      fsnt(i)     = 0.0_r8
      fsntc(i)    = 0.0_r8
      fsntoa(i)   = 0.0_r8
      fsntoac(i)  = 0.0_r8

      solin(i)    = 0.0_r8

      sols(i)     = 0.0_r8
      soll(i)     = 0.0_r8
      solsd(i)    = 0.0_r8
      solld(i)    = 0.0_r8

! initialize added downward/upward total and clear sky fluxes

         do k=1,pverp
            fsup(i,k)  = 0.0_r8
            fsupc(i,k) = 0.0_r8
            fsdn(i,k)  = 0.0_r8
            fsdnc(i,k) = 0.0_r8
            tauxcl(i,k-1) = 0.0_r8
            tauxci(i,k-1) = 0.0_r8
         end do

      do k=1, pver
         qrs(i,k) = 0.0_r8
      end do

      ! initialize aerosol diagnostic fields to 0.0 
      ! average can be obtained by dividing <aerod>/<frc_day>
      do kaer = 1, naer_groups
         do ns = 1, nspint
            frc_day(i) = 0.0_r8
            aertau(i,ns,kaer) = 0.0_r8
            aerssa(i,ns,kaer) = 0.0_r8
            aerasm(i,ns,kaer) = 0.0_r8
            aerfwd(i,ns,kaer) = 0.0_r8
         end do
      end do

   end do
! 
! compute starting, ending daytime loop indices:
!  *** note this logic assumes day and night points are contiguous so
!  *** will not work in general with chunked data structure.
! 
   ndayc = 0
   do i=1,ncol
      if (coszrs(i) > 0.0_r8) then
         ndayc = ndayc + 1
         idayc(ndayc) = i
      end if
   end do
! 
! if night everywhere, return:
! 
   if (ndayc == 0) return
! 
! perform other initializations
! 
   tmp1   = 0.5_r8/(gravit*sslp)
   tmp2   = delta/gravit
   sqrco2 = sqrt(co2mmr)

   do n=1,ndayc
      i=idayc(n)
! 
! define solar incident radiation and interface pressures:
! 
!        solin(i)  = scon*eccf*coszrs(i)
!wrf use solcon (mks) calculated outside
         solin(i)  = solcon*coszrs(i)*1000.
         pflx(i,0) = 0._r8
         do k=1,pverp
            pflx(i,k) = pint(i,k)
         end do
! 
! compute optical paths:
! 
         ptop      = pflx(i,1)
         ptho2     = o2mmr * ptop / gravit
         ptho3     = o3mmr(i,1) * ptop / gravit
         pthco2    = sqrco2 * (ptop / gravit)
         h2ostr    = sqrt( 1._r8 / h2ommr(i,1) )
         zenfac(i) = sqrt(coszrs(i))
         pthh2o    = ptop**2*tmp1 + (ptop*rga)* &
                    (h2ostr*zenfac(i)*delta)
         uh2o(i,0) = h2ommr(i,1)*pthh2o
         uco2(i,0) = zenfac(i)*pthco2
         uo2 (i,0) = zenfac(i)*ptho2
         uo3 (i,0) = ptho3
         uaer(i,0) = 0.0_r8
         do k=1,pver
            pdel      = pflx(i,k+1) - pflx(i,k)
            path      = pdel / gravit
            ptho2     = o2mmr * path
            ptho3     = o3mmr(i,k) * path
            pthco2    = sqrco2 * path
            h2ostr    = sqrt(1.0_r8/h2ommr(i,k))
            pthh2o    = (pflx(i,k+1)**2 - pflx(i,k)**2)*tmp1 + pdel*h2ostr*zenfac(i)*tmp2
            uh2o(i,k) = h2ommr(i,k)*pthh2o
            uco2(i,k) = zenfac(i)*pthco2
            uo2 (i,k) = zenfac(i)*ptho2
            uo3 (i,k) = ptho3
            usul(i,k) = aermmr(i,k,idxsul) * path 
            ubg(i,k) = aermmr(i,k,idxbg) * path 
            usslt(i,k) = aermmr(i,k,idxsslt) * path
            if (usslt(i,k) .lt. 0.0) then  ! usslt is sometimes small and negative, will be fixed
              usslt(i,k) = 0.0
            end if
            ucphil(i,k) = aermmr(i,k,idxocphi) * path
            ucphob(i,k) = aermmr(i,k,idxocpho) * path
            ucb(i,k) = ( aermmr(i,k,idxbcpho) + aermmr(i,k,idxbcphi) ) * path
            uvolc(i,k) =  aermmr(i,k,idxvolc)
            do ksz = 1, ndstsz
              udst(ksz,i,k) = aermmr(i,k,idxdustfirst-1+ksz) * path
            end do
         end do
! 
! compute column absorber amounts for the clear sky computation:
! 
         uth2o(i) = 0.0_r8
         uto3(i)  = 0.0_r8
         utco2(i) = 0.0_r8
         uto2(i)  = 0.0_r8

         do k=1,pver
            uth2o(i) = uth2o(i) + uh2o(i,k)
            uto3(i)  = uto3(i)  + uo3(i,k)
            utco2(i) = utco2(i) + uco2(i,k)
            uto2(i)  = uto2(i)  + uo2(i,k)
         end do
! 
! set cloud properties for top (0) layer; so long as tauxcl is zero,
! there is no cloud above top of model; the other cloud properties
! are arbitrary:
! 
         tauxcl(i,0)  = 0._r8
         wcl(i,0)     = 0.999999_r8
         gcl(i,0)     = 0.85_r8
         fcl(i,0)     = 0.725_r8
         tauxci(i,0)  = 0._r8
         wci(i,0)     = 0.999999_r8
         gci(i,0)     = 0.85_r8
         fci(i,0)     = 0.725_r8
! 
! aerosol 
! 
         tauxar(i,0)  = 0._r8
         wa(i,0)      = 0.925_r8
         ga(i,0)      = 0.850_r8
         fa(i,0)      = 0.7225_r8
! 
! end  do n=1,ndayc
! 
   end do
! 
! begin spectral loop
! 
   do ns=1,nspint
! 
! set index for cloud particle properties based on the wavelength,
! according to a. slingo (1989) equations 1-3:
! use index 1 (0.25 to 0.69 micrometers) for visible
! use index 2 (0.69 - 1.19 micrometers) for near-infrared
! use index 3 (1.19 to 2.38 micrometers) for near-infrared
! use index 4 (2.38 to 4.00 micrometers) for near-infrared
! 
! note that the minimum wavelength is encoded (with .001, .002, .003)
! in order to specify the index appropriate for the near-infrared
! cloud absorption properties
! 
      if(wavmax(ns) <= 0.7_r8) then
         indxsl = 1
      else if(wavmin(ns) == 0.700_r8) then
         indxsl = 2
      else if(wavmin(ns) == 0.701_r8) then
         indxsl = 3
      else if(wavmin(ns) == 0.702_r8 .or. wavmin(ns) > 2.38_r8) then
         indxsl = 4
      end if
! 
! set cloud extinction optical depth, single scatter albedo,
! asymmetry parameter, and forward scattered fraction:
! 
      abarli = abarl(indxsl)
      bbarli = bbarl(indxsl)
      cbarli = cbarl(indxsl)
      dbarli = dbarl(indxsl)
      ebarli = ebarl(indxsl)
      fbarli = fbarl(indxsl)
! 
      abarii = abari(indxsl)
      bbarii = bbari(indxsl)
      cbarii = cbari(indxsl)
      dbarii = dbari(indxsl)
      ebarii = ebari(indxsl)
      fbarii = fbari(indxsl)
! 
! adjustfraction within spectral interval to allow for the possibility of
! sub-divisions within a particular interval:
! 
      psf(ns) = 1.0_r8
      if(ph2o(ns)/=0._r8) psf(ns) = psf(ns)*ph2o(ns)
      if(pco2(ns)/=0._r8) psf(ns) = psf(ns)*pco2(ns)
      if(po2 (ns)/=0._r8) psf(ns) = psf(ns)*po2 (ns)

      do n=1,ndayc
         i=idayc(n)

         frc_day(i) = 1.0_r8
         do kaer = 1, naer_groups
            aertau(i,ns,kaer) = 0.0
            aerssa(i,ns,kaer) = 0.0
            aerasm(i,ns,kaer) = 0.0
            aerfwd(i,ns,kaer) = 0.0
         end do

            do k=1,pver
! 
! liquid
! 
               tmp1l = abarli + bbarli/rel(i,k)
               tmp2l = 1._r8 - cbarli - dbarli*rel(i,k)
               tmp3l = fbarli*rel(i,k)
! 
! ice
! 
               tmp1i = abarii + bbarii/rei(i,k)
               tmp2i = 1._r8 - cbarii - dbarii*rei(i,k)
               tmp3i = fbarii*rei(i,k)

               if (cld(i,k) >= cldmin .and. cld(i,k) >= cldeps) then
                  tauxcl(i,k) = cliqwp(i,k)*tmp1l
                  tauxci(i,k) = cicewp(i,k)*tmp1i
               else
                  tauxcl(i,k) = 0.0
                  tauxci(i,k) = 0.0
               endif
! 
! do not let single scatter albedo be 1.  delta-eddington solution
! for non-conservative case has different analytic form from solution
! for conservative case, and raddedmx is written for non-conservative case.
! 
               wcl(i,k) = min(tmp2l,.999999_r8)
               gcl(i,k) = ebarli + tmp3l
               fcl(i,k) = gcl(i,k)*gcl(i,k)
! 
               wci(i,k) = min(tmp2i,.999999_r8)
               gci(i,k) = ebarii + tmp3i
               fci(i,k) = gci(i,k)*gci(i,k)
! 
! set aerosol properties
! conversion factor to adjust aerosol extinction (m2/g)
! 
               rhtrunc = rh(i,k)
               rhtrunc = min(rh(i,k),1._r8)
!              if(rhtrunc.lt.0._r8) call endrun ('radcswmx')
               krh = min(floor( rhtrunc * nrh ) + 1, nrh - 1)
               wrh = rhtrunc * nrh - krh

               ! linear interpolation of optical properties between rh table points
               ksuli = ksul(krh + 1, ns) * (wrh + 1) - ksul(krh, ns) * wrh
               ksslti = ksslt(krh + 1, ns) * (wrh + 1) - ksslt(krh, ns) * wrh
               kcphili = kcphil(krh + 1, ns) * (wrh + 1) - kcphil(krh, ns) * wrh
               wsuli = wsul(krh + 1, ns) * (wrh + 1) - wsul(krh, ns) * wrh
               wsslti = wsslt(krh + 1, ns) * (wrh + 1) - wsslt(krh, ns) * wrh
               wcphili = wcphil(krh + 1, ns) * (wrh + 1) - wcphil(krh, ns) * wrh
               gsuli = gsul(krh + 1, ns) * (wrh + 1) - gsul(krh, ns) * wrh
               gsslti = gsslt(krh + 1, ns) * (wrh + 1) - gsslt(krh, ns) * wrh
               gcphili = gcphil(krh + 1, ns) * (wrh + 1) - gcphil(krh, ns) * wrh

               tau_sul = 1.e4 * ksuli * usul(i,k)
               tau_sslt = 1.e4 * ksslti * usslt(i,k)
               tau_cphil = 1.e4 * kcphili * ucphil(i,k)
               tau_cphob = 1.e4 * kcphob(ns) * ucphob(i,k)
               tau_cb = 1.e4 * kcb(ns) * ucb(i,k)
               tau_volc = 1.e3 * kvolc(ns) * uvolc(i,k)
               tau_dst(:) = 1.e4 * kdst(:,ns) * udst(:,i,k)
               tau_bg = 1.e4 * kbg(ns) * ubg(i,k)

               tau_w_sul = tau_sul * wsuli
               tau_w_sslt = tau_sslt * wsslti
               tau_w_cphil = tau_cphil * wcphili
               tau_w_cphob = tau_cphob * wcphob(ns)
               tau_w_cb = tau_cb * wcb(ns)
               tau_w_volc = tau_volc * wvolc(ns)
               tau_w_dst(:) = tau_dst(:) * wdst(:,ns)
               tau_w_bg = tau_bg * wbg(ns)

               tau_w_g_sul = tau_w_sul * gsuli
               tau_w_g_sslt = tau_w_sslt * gsslti
               tau_w_g_cphil = tau_w_cphil * gcphili
               tau_w_g_cphob = tau_w_cphob * gcphob(ns)
               tau_w_g_cb = tau_w_cb * gcb(ns)
               tau_w_g_volc = tau_w_volc * gvolc(ns)
               tau_w_g_dst(:) = tau_w_dst(:) * gdst(:,ns)
               tau_w_g_bg = tau_w_bg * gbg(ns)

               f_sul = gsuli * gsuli
               f_sslt = gsslti * gsslti
               f_cphil = gcphili * gcphili
               f_cphob = gcphob(ns) * gcphob(ns)
               f_cb = gcb(ns) * gcb(ns)
               f_volc = gvolc(ns) * gvolc(ns)
               f_dst(:) = gdst(:,ns) * gdst(:,ns)
               f_bg = gbg(ns) * gbg(ns)

               tau_w_f_sul = tau_w_sul * f_sul
               tau_w_f_bg = tau_w_bg * f_bg
               tau_w_f_sslt = tau_w_sslt * f_sslt
               tau_w_f_cphil = tau_w_cphil * f_cphil
               tau_w_f_cphob = tau_w_cphob * f_cphob
               tau_w_f_cb = tau_w_cb * f_cb
               tau_w_f_volc = tau_w_volc * f_volc
               tau_w_f_dst(:) = tau_w_dst(:) * f_dst(:)
!
! mix dust aerosol size bins
!   w_dst_tot, g_dst_tot, w_dst_tot are currently not used anywhere
!   but calculate them anyway for future use
!
               tau_dst_tot = sum(tau_dst)
               tau_w_dst_tot = sum(tau_w_dst)
               tau_w_g_dst_tot = sum(tau_w_g_dst)
               tau_w_f_dst_tot = sum(tau_w_f_dst)

               if (tau_dst_tot .gt. 0.0) then
                 w_dst_tot = tau_w_dst_tot / tau_dst_tot
               else
                 w_dst_tot = 0.0
               endif

               if (tau_w_dst_tot .gt. 0.0) then
                 g_dst_tot = tau_w_g_dst_tot / tau_w_dst_tot
                 f_dst_tot = tau_w_f_dst_tot / tau_w_dst_tot
               else
                 g_dst_tot = 0.0
                 f_dst_tot = 0.0
               endif
!
! mix aerosols
!
               tau_tot     = tau_sul + tau_sslt &
                           + tau_cphil + tau_cphob + tau_cb + tau_dst_tot
               tau_tot     = tau_tot + tau_bg + tau_volc

               tau_w_tot   = tau_w_sul + tau_w_sslt &
                           + tau_w_cphil + tau_w_cphob + tau_w_cb + tau_w_dst_tot
               tau_w_tot   = tau_w_tot + tau_w_bg + tau_w_volc

               tau_w_g_tot = tau_w_g_sul + tau_w_g_sslt &
                           + tau_w_g_cphil + tau_w_g_cphob + tau_w_g_cb + tau_w_g_dst_tot
               tau_w_g_tot = tau_w_g_tot + tau_w_g_bg + tau_w_g_volc

               tau_w_f_tot = tau_w_f_sul + tau_w_f_sslt &
                           + tau_w_f_cphil + tau_w_f_cphob + tau_w_f_cb + tau_w_f_dst_tot
               tau_w_f_tot = tau_w_f_tot + tau_w_f_bg  + tau_w_f_volc

               if (tau_tot .gt. 0.0) then
                 w_tot = tau_w_tot / tau_tot
               else
                 w_tot = 0.0
               endif

               if (tau_w_tot .gt. 0.0) then
                 g_tot = tau_w_g_tot / tau_w_tot
                 f_tot = tau_w_f_tot / tau_w_tot
               else
                 g_tot = 0.0
                 f_tot = 0.0
               endif

               tauxar(i,k) = tau_tot
               wa(i,k)     = min(w_tot, 0.999999_r8)
               if (g_tot.gt.1._r8) write(6,*) "g_tot > 1"
               if (g_tot.lt.-1._r8) write(6,*) "g_tot < -1"
!              if (g_tot.gt.1._r8) call endrun ('radcswmx')
!              if (g_tot.lt.-1._r8) call endrun ('radcswmx')
               ga(i,k)     = g_tot
               if (f_tot.gt.1._r8) write(6,*)"f_tot > 1"
               if (f_tot.lt.0._r8) write(6,*)"f_tot < 0"
!              if (f_tot.gt.1._r8) call endrun ('radcswmx')
!              if (f_tot.lt.0._r8) call endrun ('radcswmx')
               fa(i,k)     = f_tot

               aertau(i,ns,1) = aertau(i,ns,1) + tau_sul
               aertau(i,ns,2) = aertau(i,ns,2) + tau_sslt
               aertau(i,ns,3) = aertau(i,ns,3) + tau_cphil + tau_cphob + tau_cb
               aertau(i,ns,4) = aertau(i,ns,4) + tau_dst_tot
               aertau(i,ns,5) = aertau(i,ns,5) + tau_bg
               aertau(i,ns,6) = aertau(i,ns,6) + tau_volc
               aertau(i,ns,7) = aertau(i,ns,7) + tau_tot

               aerssa(i,ns,1) = aerssa(i,ns,1) + tau_w_sul
               aerssa(i,ns,2) = aerssa(i,ns,2) + tau_w_sslt
               aerssa(i,ns,3) = aerssa(i,ns,3) + tau_w_cphil + tau_w_cphob + tau_w_cb
               aerssa(i,ns,4) = aerssa(i,ns,4) + tau_w_dst_tot
               aerssa(i,ns,5) = aerssa(i,ns,5) + tau_w_bg
               aerssa(i,ns,6) = aerssa(i,ns,6) + tau_w_volc
               aerssa(i,ns,7) = aerssa(i,ns,7) + tau_w_tot

               aerasm(i,ns,1) = aerasm(i,ns,1) + tau_w_g_sul
               aerasm(i,ns,2) = aerasm(i,ns,2) + tau_w_g_sslt
               aerasm(i,ns,3) = aerasm(i,ns,3) + tau_w_g_cphil + tau_w_g_cphob + tau_w_g_cb
               aerasm(i,ns,4) = aerasm(i,ns,4) + tau_w_g_dst_tot
               aerasm(i,ns,5) = aerasm(i,ns,5) + tau_w_g_bg
               aerasm(i,ns,6) = aerasm(i,ns,6) + tau_w_g_volc
               aerasm(i,ns,7) = aerasm(i,ns,7) + tau_w_g_tot

               aerfwd(i,ns,1) = aerfwd(i,ns,1) + tau_w_f_sul
               aerfwd(i,ns,2) = aerfwd(i,ns,2) + tau_w_f_sslt
               aerfwd(i,ns,3) = aerfwd(i,ns,3) + tau_w_f_cphil + tau_w_f_cphob + tau_w_f_cb
               aerfwd(i,ns,4) = aerfwd(i,ns,4) + tau_w_f_dst_tot
               aerfwd(i,ns,5) = aerfwd(i,ns,5) + tau_w_f_bg
               aerfwd(i,ns,6) = aerfwd(i,ns,6) + tau_w_f_volc
               aerfwd(i,ns,7) = aerfwd(i,ns,7) + tau_w_f_tot

! 
! end do k=1,pver
! 
            end do

            ! normalize aerosol optical diagnostic fields
            do kaer = 1, naer_groups

               if (aerssa(i,ns,kaer) .gt. 0.0) then   ! aerssa currently holds product of tau and ssa
                  aerasm(i,ns,kaer) = aerasm(i,ns,kaer) / aerssa(i,ns,kaer)
                  aerfwd(i,ns,kaer) = aerfwd(i,ns,kaer) / aerssa(i,ns,kaer)
               else
                  aerasm(i,ns,kaer) = 0.0_r8
                  aerfwd(i,ns,kaer) = 0.0_r8
               end if

               if (aertau(i,ns,kaer) .gt. 0.0) then
                  aerssa(i,ns,kaer) = aerssa(i,ns,kaer) / aertau(i,ns,kaer)
               else
                  aerssa(i,ns,kaer) = 0.0_r8
               end if

            end do


! 
! end do n=1,ndayc
! 
      end do

! 
! set reflectivities for surface based on mid-point wavelength
! 
      wavmid(ns) = 0.5_r8*(wavmin(ns) + wavmax(ns))
! 
! wavelength less  than 0.7 micro-meter
! 
      if (wavmid(ns) < 0.7_r8 ) then
         do n=1,ndayc
            i=idayc(n)
               albdir(i,ns) = asdir(i)
               albdif(i,ns) = asdif(i)
         end do
! 
! wavelength greater than 0.7 micro-meter
! 
      else
         do n=1,ndayc
            i=idayc(n)
               albdir(i,ns) = aldir(i)
               albdif(i,ns) = aldif(i)
         end do
      end if
      trayoslp = raytau(ns)/sslp
! 
! layer input properties now completely specified; compute the
! delta-eddington solution reflectivities and transmissivities
! for each layer
! 
      call raddedmx(pver, pverp, pcols, coszrs   ,ndayc    ,idayc   , &
              abh2o(ns),abo3(ns) ,abco2(ns),abo2(ns) , &
              uh2o     ,uo3      ,uco2     ,uo2      , &
              trayoslp ,pflx     ,ns       , &
              tauxcl   ,wcl      ,gcl      ,fcl      , &
              tauxci   ,wci      ,gci      ,fci      , &
              tauxar   ,wa       ,ga       ,fa       , &
              rdir     ,rdif     ,tdir     ,tdif     ,explay  , &
              rdirc    ,rdifc    ,tdirc    ,tdifc    ,explayc )
! 
! end spectral loop
! 
   end do
! 
!----------------------------------------------------------------------
! 
! solution for max/random cloud overlap.  
! 
! steps:
! (1. delta-eddington solution for each layer (called above)
! 
! (2. the adding method is used to
! compute the reflectivity and transmissivity to direct and diffuse
! radiation from the top and bottom of the atmosphere for each
! cloud configuration.  this calculation is based upon the
! max-random overlap assumption.
! 
! (3. to solve for the fluxes, combine the
! bulk properties of the atmosphere above/below the region.
! 
! index calculations for steps 2-3 are performed outside spectral
! loop to avoid redundant calculations.  index calculations (with
! application of areamin & nconfgmax conditions) are performed 
! first to identify the minimum subset of terms for the configurations 
! satisfying the areamin & nconfgmax conditions. this minimum set is 
! used to identify the corresponding minimum subset of terms in 
! steps 2 and 3.
! 

   do n=1,ndayc
      i=idayc(n)

!----------------------------------------------------------------------
! index calculations for max overlap
! 
! the column is divided into sets of adjacent layers, called regions, 
! in which the clouds are maximally overlapped.  the clouds are
! randomly overlapped between different regions.  the number of
! regions in a column is set by nmxrgn, and the range of pressures
! included in each region is set by pmxrgn.  
! 
! the following calculations determine the number of unique cloud 
! configurations (assuming maximum overlap), called "streams",
! within each region. each stream consists of a vector of binary
! clouds (either 0 or 100% cloud cover).  over the depth of the region, 
! each stream requires a separate calculation of radiative properties. these
! properties are generated using the adding method from
! the radiative properties for each layer calculated by raddedmx.
! 
! the upward and downward-propagating streams are treated
! separately.
! 
! we will refer to a particular configuration of binary clouds
! within a single max-overlapped region as a "stream".  we will 
! refer to a particular arrangement of binary clouds over the entire column
! as a "configuration".
! 
! this section of the code generates the following information:
! (1. nrgn    : the true number of max-overlap regions (need not = nmxrgn)
! (2. nstr    : the number of streams in a region (>=1)
! (3. cstr    : flags for presence of clouds at each layer in each stream
! (4. wstr    : the fractional horizontal area of a grid box covered
! by each stream
! (5. kx1,2   : level indices for top/bottom of each region
! 
! the max-overlap calculation proceeds in 3 stages:
! (1. compute layer radiative properties in raddedmx.
! (2. combine these properties between layers 
! (3. combine properties to compute fluxes at each interface.  
! 
! most of the indexing information calculated here is used in steps 2-3
! after the call to raddedmx.
! 
! initialize indices for layers to be max-overlapped
! 
! loop to handle fix in totwgt=0. for original overlap config 
! from npasses = 0.
! 
         npasses = 0
         do
            do irgn = 0, nmxrgn(i)
               kx2(irgn) = 0
            end do
            mrgn = 0
! 
! outermost loop over regions (sets of adjacent layers) to be max overlapped
! 
            do irgn = 1, nmxrgn(i)
! 
! calculate min/max layer indices inside region.  
! 
               region_found = .false.
               if (kx2(irgn-1) < pver) then
                  k1 = kx2(irgn-1)+1
                  kx1(irgn) = k1
                  kx2(irgn) = k1-1
                  do k2 = pver, k1, -1
                     if (pmid(i,k2) <= pmxrgn(i,irgn)) then
                        kx2(irgn) = k2
                        mrgn = mrgn+1
                        region_found = .true.
                        exit
                     end if
                  end do
               else
                  exit
               endif

               if (region_found) then
! 
! sort cloud areas and corresponding level indices.  
! 
                  nxs = 0
                  if (cldeps > 0) then 
                     do k = k1,k2
                        if (cld(i,k) >= cldmin .and. cld(i,k) >= cldeps) then
                           nxs = nxs+1
                           ksort(nxs) = k
! 
! we need indices for clouds in order of largest to smallest, so
! sort 1-cld in ascending order
! 
                           asort(nxs) = 1.0_r8-(floor(cld(i,k)/cldeps)*cldeps)
                        end if
                     end do
                  else
                     do k = k1,k2
                        if (cld(i,k) >= cldmin) then
                           nxs = nxs+1
                           ksort(nxs) = k
! 
! we need indices for clouds in order of largest to smallest, so
! sort 1-cld in ascending order
! 
                           asort(nxs) = 1.0_r8-cld(i,k)
                        end if
                     end do
                  endif
! 
! if nxs eq 1, no need to sort. 
! if nxs eq 2, sort by swapping if necessary
! if nxs ge 3, sort using local sort routine
! 
                  if (nxs == 2) then
                     if (asort(2) < asort(1)) then
                        ktmp = ksort(1)
                        ksort(1) = ksort(2)
                        ksort(2) = ktmp

                        atmp = asort(1)
                        asort(1) = asort(2)
                        asort(2) = atmp
                     endif
                  else if (nxs >= 3) then
                     call sortarray(nxs,asort,ksort)
                  endif
! 
! construct wstr, cstr, nstr for this region
! 
                  cstr(k1:k2,1:nxs+1) = 0
                  mstr = 1
                  cld0 = 0.0_r8
                  do l = 1, nxs
                     if (asort(l) /= cld0) then
                        wstr(mstr,mrgn) = asort(l) - cld0
                        cld0 = asort(l)
                        mstr = mstr + 1
                     endif
                     cstr(ksort(l),mstr:nxs+1) = 1
                  end do
                  nstr(mrgn) = mstr
                  wstr(mstr,mrgn) = 1.0_r8 - cld0
! 
! end test of region_found = true
! 
               endif
! 
! end loop over regions irgn for max-overlap
! 
            end do
            nrgn = mrgn
! 
! finish construction of cstr for additional top layer
! 
            cstr(0,1:nstr(1)) = 0
! 
! index computations for step 2-3
! this section of the code generates the following information:
! (1. totwgt     step 3     total frac. area of configurations satisfying
! areamin & nconfgmax criteria
! (2. wgtv       step 3     frac. area of configurations 
! (3. ccon       step 2     binary flag for clouds in each configuration
! (4. nconfig    steps 2-3  number of configurations
! (5. nuniqu/d   step 2     number of unique cloud configurations for
! up/downwelling rad. between surface/toa
! and level k
! (6. istrtu/d   step 2     indices into iconu/d
! (7. iconu/d    step 2     cloud configurations which are identical
! for up/downwelling rad. between surface/toa
! and level k
! 
! number of configurations (all permutations of streams in each region)
! 
            nconfigm = product(nstr(1: nrgn))
! 
! construction of totwgt, wgtv, ccon, nconfig
! 
            istr(1: nrgn) = 1
            nconfig = 0
            totwgt = 0.0_r8
            new_term = .true.
            do iconfig = 1, nconfigm
               xwgt = 1.0_r8
               do mrgn = 1,  nrgn
                  xwgt = xwgt * wstr(istr(mrgn),mrgn)
               end do
               if (xwgt >= areamin) then
                  nconfig = nconfig + 1
                  if (nconfig <= nconfgmax) then
                     j = nconfig
                     ptrc(nconfig) = nconfig
                  else
                     nconfig = nconfgmax
                     if (new_term) then
                        j = findvalue(1,nconfig,wgtv,ptrc)
                     endif
                     if (wgtv(j) < xwgt) then
                        totwgt = totwgt - wgtv(j)
                        new_term = .true.
                     else
                        new_term = .false.
                     endif
                  endif
                  if (new_term) then
                     wgtv(j) = xwgt
                     totwgt = totwgt + xwgt
                     do mrgn = 1, nrgn
                        ccon(kx1(mrgn):kx2(mrgn),j) = cstr(kx1(mrgn):kx2(mrgn),istr(mrgn))
                     end do
                  endif
               endif

               mrgn =  nrgn
               istr(mrgn) = istr(mrgn) + 1
               do while (istr(mrgn) > nstr(mrgn) .and. mrgn > 1)
                  istr(mrgn) = 1
                  mrgn = mrgn - 1
                  istr(mrgn) = istr(mrgn) + 1
               end do
! 
! end do iconfig = 1, nconfigm
! 
            end do
! 
! if totwgt = 0 implement maximum overlap and make another pass
! if totwgt = 0 on this second pass then terminate.
! 
            if (totwgt > 0.) then
               exit
            else
               npasses = npasses + 1
               if (npasses >= 2 ) then
                  write(6,*)'radcswmx: maximum overlap of column ','failed'
                  call endrun
               endif
               nmxrgn(i)=1
               pmxrgn(i,1)=1.0e30
            end if
!
! end npasses = 0, do
!
         end do
! 
! 
! finish construction of ccon
! 
         ccon(0,:) = 0
         ccon(pverp,:) = 0
! 
! construction of nuniqu/d, istrtu/d, iconu/d using binary tree 
! 
         nuniqd(0) = 1
         nuniqu(pverp) = 1

         istrtd(0,1) = 1
         istrtu(pverp,1) = 1

         do j = 1, nconfig
            icond(0,j)=j
            iconu(pverp,j)=j
         end do

         istrtd(0,2) = nconfig+1
         istrtu(pverp,2) = nconfig+1

         do k = 1, pverp
            km1 = k-1
            nuniq = 0
            istrtd(k,1) = 1
            do l0 = 1, nuniqd(km1)
               is0 = istrtd(km1,l0)
               is1 = istrtd(km1,l0+1)-1
               n0 = 0
               n1 = 0
               do isn = is0, is1
                  j = icond(km1,isn)
                  if (ccon(k,j) == 0) then
                     n0 = n0 + 1
                     ptr0(n0) = j
                  endif
                  if (ccon(k,j) == 1) then
                     n1 = n1 + 1
                     ptr1(n1) = j
                  endif
               end do
               if (n0 > 0) then
                  nuniq = nuniq + 1
                  istrtd(k,nuniq+1) = istrtd(k,nuniq)+n0
                  icond(k,istrtd(k,nuniq):istrtd(k,nuniq+1)-1) =  ptr0(1:n0)
               endif
               if (n1 > 0) then
                  nuniq = nuniq + 1
                  istrtd(k,nuniq+1) = istrtd(k,nuniq)+n1
                  icond(k,istrtd(k,nuniq):istrtd(k,nuniq+1)-1) =  ptr1(1:n1)
               endif
            end do
            nuniqd(k) = nuniq
         end do

         do k = pver, 0, -1
            kp1 = k+1
            nuniq = 0
            istrtu(k,1) = 1
            do l0 = 1, nuniqu(kp1)
               is0 = istrtu(kp1,l0)
               is1 = istrtu(kp1,l0+1)-1
               n0 = 0
               n1 = 0
               do isn = is0, is1
                  j = iconu(kp1,isn)
                  if (ccon(k,j) == 0) then
                     n0 = n0 + 1
                     ptr0(n0) = j
                  endif
                  if (ccon(k,j) == 1) then
                     n1 = n1 + 1
                     ptr1(n1) = j
                  endif
               end do
               if (n0 > 0) then
                  nuniq = nuniq + 1
                  istrtu(k,nuniq+1) = istrtu(k,nuniq)+n0
                  iconu(k,istrtu(k,nuniq):istrtu(k,nuniq+1)-1) =  ptr0(1:n0)
               endif
               if (n1 > 0) then
                  nuniq = nuniq + 1
                  istrtu(k,nuniq+1) = istrtu(k,nuniq)+n1
                  iconu(k,istrtu(k,nuniq):istrtu(k,nuniq+1)-1) = ptr1(1:n1)
               endif
            end do
            nuniqu(k) = nuniq
         end do
! 
!----------------------------------------------------------------------
! end of index calculations
!----------------------------------------------------------------------


!----------------------------------------------------------------------
! start of flux calculations
!----------------------------------------------------------------------
! 
! initialize spectrally integrated totals:
! 
         do k=0,pver
            totfld(k) = 0.0_r8
            fswup (k) = 0.0_r8
            fswdn (k) = 0.0_r8
            fswupc (k) = 0.0_r8
            fswdnc (k) = 0.0_r8
         end do

         sfltot        = 0.0_r8
         fswup (pverp) = 0.0_r8
         fswdn (pverp) = 0.0_r8
         fswupc (pverp) = 0.0_r8
         fswdnc (pverp) = 0.0_r8
! 
! start spectral interval
! 
         do ns = 1,nspint
            wgtint = nirwgt(ns)
!----------------------------------------------------------------------
! step 2
! 
! 
! apply adding method to solve for radiative properties
! 
! first initialize the bulk properties at toa
! 
            rdndif(0,1:nconfig) = 0.0_r8
            exptdn(0,1:nconfig) = 1.0_r8
            tdntot(0,1:nconfig) = 1.0_r8
! 
! solve for properties involving downward propagation of radiation.
! the bulk properties are:
! 
! (1. exptdn   sol. beam dwn. trans from layers above
! (2. rdndif   ref to dif rad for layers above
! (3. tdntot   total trans for layers above
! 
            do k = 1, pverp
               km1 = k - 1
               do l0 = 1, nuniqd(km1)
                  is0 = istrtd(km1,l0)
                  is1 = istrtd(km1,l0+1)-1

                  j = icond(km1,is0)

                  xexpt   = exptdn(km1,j)
                  xrdnd   = rdndif(km1,j)
                  tdnmexp = tdntot(km1,j) - xexpt

                  if (ccon(km1,j) == 1) then
! 
! if cloud in layer, use cloudy layer radiative properties
! 
                     ytdnd = tdif(ns,i,km1)
                     yrdnd = rdif(ns,i,km1)

                     rdenom  = 1._r8/(1._r8-yrdnd*xrdnd)
                     rdirexp = rdir(ns,i,km1)*xexpt

                     zexpt = xexpt * explay(ns,i,km1)
                     zrdnd = yrdnd + xrdnd*(ytdnd**2)*rdenom
                     ztdnt = xexpt*tdir(ns,i,km1) + ytdnd*(tdnmexp + xrdnd*rdirexp)*rdenom
                  else
! 
! if clear layer, use clear-sky layer radiative properties
! 
                     ytdnd = tdifc(ns,i,km1)
                     yrdnd = rdifc(ns,i,km1)

                     rdenom  = 1._r8/(1._r8-yrdnd*xrdnd)
                     rdirexp = rdirc(ns,i,km1)*xexpt

                     zexpt = xexpt * explayc(ns,i,km1)
                     zrdnd = yrdnd + xrdnd*(ytdnd**2)*rdenom
                     ztdnt = xexpt*tdirc(ns,i,km1) + ytdnd* &
                                            (tdnmexp + xrdnd*rdirexp)*rdenom
                  endif

! 
! if 2 or more configurations share identical properties at a given level k,
! the properties (at level k) are computed once and copied to 
! all the configurations for efficiency.
! 
                  do isn = is0, is1
                     j = icond(km1,isn)
                     exptdn(k,j) = zexpt
                     rdndif(k,j) = zrdnd
                     tdntot(k,j) = ztdnt
                  end do
! 
! end do l0 = 1, nuniqd(k)
! 
               end do
! 
! end do k = 1, pverp
! 
            end do
! 
! solve for properties involving upward propagation of radiation.
! the bulk properties are:
! 
! (1. rupdif   ref to dif rad for layers below
! (2. rupdir   ref to dir rad for layers below
! 
! specify surface boundary conditions (surface albedos)
! 
            rupdir(pverp,1:nconfig) = albdir(i,ns)
            rupdif(pverp,1:nconfig) = albdif(i,ns)

            do k = pver, 0, -1
               do l0 = 1, nuniqu(k)
                  is0 = istrtu(k,l0)
                  is1 = istrtu(k,l0+1)-1

                  j = iconu(k,is0)

                  xrupd = rupdif(k+1,j)
                  xrups = rupdir(k+1,j)

                  if (ccon(k,j) == 1) then
! 
! if cloud in layer, use cloudy layer radiative properties
! 
                     yexpt = explay(ns,i,k)
                     yrupd = rdif(ns,i,k)
                     ytupd = tdif(ns,i,k)

                     rdenom  = 1._r8/( 1._r8 - yrupd*xrupd)
                     tdnmexp = (tdir(ns,i,k)-yexpt)
                     rdirexp = xrups*yexpt

                     zrupd = yrupd + xrupd*(ytupd**2)*rdenom
                     zrups = rdir(ns,i,k) + ytupd*(rdirexp + xrupd*tdnmexp)*rdenom
                  else
! 
! if clear layer, use clear-sky layer radiative properties
! 
                     yexpt = explayc(ns,i,k)
                     yrupd = rdifc(ns,i,k)
                     ytupd = tdifc(ns,i,k)

                     rdenom  = 1._r8/( 1._r8 - yrupd*xrupd)
                     tdnmexp = (tdirc(ns,i,k)-yexpt)
                     rdirexp = xrups*yexpt

                     zrupd = yrupd + xrupd*(ytupd**2)*rdenom
                     zrups = rdirc(ns,i,k) + ytupd*(rdirexp + xrupd*tdnmexp)*rdenom
                  endif

! 
! if 2 or more configurations share identical properties at a given level k,
! the properties (at level k) are computed once and copied to 
! all the configurations for efficiency.
! 
                  do isn = is0, is1
                     j = iconu(k,isn)
                     rupdif(k,j) = zrupd
                     rupdir(k,j) = zrups
                  end do
! 
! end do l0 = 1, nuniqu(k)
! 
               end do
! 
! end do k = pver,0,-1
! 
            end do
! 
!----------------------------------------------------------------------
! 
! step 3
! 
! compute up and down fluxes for each interface k.  this requires
! adding up the contributions from all possible permutations
! of streams in all max-overlap regions, weighted by the
! product of the fractional areas of the streams in each region
! (the random overlap assumption).  the adding principle has been
! used in step 2 to combine the bulk radiative properties 
! above and below the interface.
! 
            do k = 0,pverp
! 
! initialize the fluxes
! 
               fluxup(k)=0.0_r8
               fluxdn(k)=0.0_r8

               do iconfig = 1, nconfig
                  xwgt = wgtv(iconfig)
                  xexpt = exptdn(k,iconfig)
                  xtdnt = tdntot(k,iconfig)
                  xrdnd = rdndif(k,iconfig)
                  xrupd = rupdif(k,iconfig)
                  xrups = rupdir(k,iconfig)
! 
! flux computation
! 
                  rdenom = 1._r8/(1._r8 - xrdnd * xrupd)

                  fluxup(k) = fluxup(k) + xwgt *  &
                              ((xexpt * xrups + (xtdnt - xexpt) * xrupd) * rdenom)
                  fluxdn(k) = fluxdn(k) + xwgt *  &
                              (xexpt + (xtdnt - xexpt + xexpt * xrups * xrdnd) * rdenom)
! 
! end do iconfig = 1, nconfig
! 
               end do
! 
! normalize by total area covered by cloud configurations included
! in solution
! 
               fluxup(k)=fluxup(k) / totwgt
               fluxdn(k)=fluxdn(k) / totwgt                  
! 
! end do k = 0,pverp
! 
            end do
! 
! initialize the direct-beam flux at surface
! 
            wexptdn = 0.0_r8

            do iconfig = 1, nconfig
               wexptdn =  wexptdn + wgtv(iconfig) * exptdn(pverp,iconfig)
            end do

            wexptdn = wexptdn / totwgt
! 
! monochromatic computation completed; accumulate in totals
! 
            solflx   = solin(i)*frcsol(ns)*psf(ns)
            fsnt(i)  = fsnt(i) + solflx*(fluxdn(1) - fluxup(1))
            fsntoa(i)= fsntoa(i) + solflx*(fluxdn(0) - fluxup(0))
            fsns(i)  = fsns(i) + solflx*(fluxdn(pverp)-fluxup(pverp))
            sfltot   = sfltot + solflx
            fswup(0) = fswup(0) + solflx*fluxup(0)
            fswdn(0) = fswdn(0) + solflx*fluxdn(0)
! 
! down spectral fluxes need to be in mks; thus the .001 conversion factors
! 
            if (wavmid(ns) < 0.7_r8) then
               sols(i)  = sols(i) + wexptdn*solflx*0.001_r8
               solsd(i) = solsd(i)+(fluxdn(pverp)-wexptdn)*solflx*0.001_r8
            else
               soll(i)  = soll(i) + wexptdn*solflx*0.001_r8
               solld(i) = solld(i)+(fluxdn(pverp)-wexptdn)*solflx*0.001_r8
               fsnrtoaq(i) = fsnrtoaq(i) + solflx*(fluxdn(0) - fluxup(0))
            end if
            fsnirtoa(i) = fsnirtoa(i) + wgtint*solflx*(fluxdn(0) - fluxup(0))

            do k=0,pver
! 
! compute flux divergence in each layer using the interface up and down
! fluxes:
! 
               kp1 = k+1
               flxdiv = (fluxdn(k  ) - fluxdn(kp1)) + (fluxup(kp1) - fluxup(k  ))
               totfld(k)  = totfld(k)  + solflx*flxdiv
               fswdn(kp1) = fswdn(kp1) + solflx*fluxdn(kp1)
               fswup(kp1) = fswup(kp1) + solflx*fluxup(kp1)
            end do
! 
! perform clear-sky calculation
! 
            exptdnc(0) =   1.0_r8
            rdndifc(0) =   0.0_r8
            tdntotc(0) =   1.0_r8
            rupdirc(pverp) = albdir(i,ns)
            rupdifc(pverp) = albdif(i,ns)

            do k = 1, pverp
               km1 = k - 1
               xexpt = exptdnc(km1)
               xrdnd = rdndifc(km1)
               yrdnd = rdifc(ns,i,km1)
               ytdnd = tdifc(ns,i,km1)

               exptdnc(k) = xexpt*explayc(ns,i,km1)

               rdenom  = 1._r8/(1._r8 - yrdnd*xrdnd)
               rdirexp = rdirc(ns,i,km1)*xexpt
               tdnmexp = tdntotc(km1) - xexpt

               tdntotc(k) = xexpt*tdirc(ns,i,km1) + ytdnd*(tdnmexp + xrdnd*rdirexp)* &
                                rdenom
               rdndifc(k) = yrdnd + xrdnd*(ytdnd**2)*rdenom
            end do

            do k=pver,0,-1
               xrupd = rupdifc(k+1)
               yexpt = explayc(ns,i,k)
               yrupd = rdifc(ns,i,k)
               ytupd = tdifc(ns,i,k)

               rdenom = 1._r8/( 1._r8 - yrupd*xrupd)

               rupdirc(k) = rdirc(ns,i,k) + ytupd*(rupdirc(k+1)*yexpt + &
                            xrupd*(tdirc(ns,i,k)-yexpt))*rdenom
               rupdifc(k) = yrupd + xrupd*ytupd**2*rdenom
            end do

            do k=0,1
               rdenom    = 1._r8/(1._r8 - rdndifc(k)*rupdifc(k))
               fluxup(k) = (exptdnc(k)*rupdirc(k) + (tdntotc(k)-exptdnc(k))*rupdifc(k))* &
                           rdenom
               fluxdn(k) = exptdnc(k) + &
                           (tdntotc(k) - exptdnc(k) + exptdnc(k)*rupdirc(k)*rdndifc(k))* &
                           rdenom
               fswupc(k) = fswupc(k) + solflx*fluxup(k)
               fswdnc(k) = fswdnc(k) + solflx*fluxdn(k)
            end do
!           k = pverp
            do k=2,pverp
            rdenom      = 1._r8/(1._r8 - rdndifc(k)*rupdifc(k))
            fluxup(k)   = (exptdnc(k)*rupdirc(k) + (tdntotc(k)-exptdnc(k))*rupdifc(k))* &
                           rdenom
            fluxdn(k)   = exptdnc(k) + (tdntotc(k) - exptdnc(k) + &
                          exptdnc(k)*rupdirc(k)*rdndifc(k))*rdenom
            fswupc(k)   = fswupc(k) + solflx*fluxup(k)
            fswdnc(k)   = fswdnc(k) + solflx*fluxdn(k)
            end do

            fsntc(i)    = fsntc(i)+solflx*(fluxdn(1)-fluxup(1))
            fsntoac(i)  = fsntoac(i)+solflx*(fluxdn(0)-fluxup(0))
            fsnsc(i)    = fsnsc(i)+solflx*(fluxdn(pverp)-fluxup(pverp))
            fsdsc(i)    = fsdsc(i)+solflx*(fluxdn(pverp))
            fsnrtoac(i) = fsnrtoac(i)+wgtint*solflx*(fluxdn(0)-fluxup(0))
! 
! end of clear sky calculation
! 

! 
! end of spectral interval loop
! 
         end do
! 
! compute solar heating rate (j/kg/s)
! 
         do k=1,pver
            qrs(i,k) = -1.e-4*gravit*totfld(k)/(pint(i,k) - pint(i,k+1))
         end do

! added downward/upward total and clear sky fluxes

         do k=1,pverp
            fsup(i,k)  = fswup(k)
            fsupc(i,k) = fswupc(k)
            fsdn(i,k)  = fswdn(k)
            fsdnc(i,k) = fswdnc(k)
         end do
! 
! set the downwelling flux at the surface 
! 
         fsds(i) = fswdn(pverp)
! 
! end do n=1,ndayc
! 
   end do

!  write (6, '(a, x, i3)') 'radcswmx : exiting, chunk identifier', lchnk

   return
end subroutine radcswmx

subroutine raddedmx(pver, pverp, pcols, coszrs  ,ndayc   ,idayc   ,abh2o   , &
                    abo3    ,abco2   ,abo2    ,uh2o    ,uo3     , &
                    uco2    ,uo2     ,trayoslp,pflx    ,ns      , &
                    tauxcl  ,wcl     ,gcl     ,fcl     ,tauxci  , &
                    wci     ,gci     ,fci     ,tauxar  ,wa      , &
                    ga      ,fa      ,rdir    ,rdif    ,tdir    , &
                    tdif    ,explay  ,rdirc   ,rdifc   ,tdirc   , &
                    tdifc   ,explayc )
!----------------------------------------------------------------------- 
! 
! purpose: 
! computes layer reflectivities and transmissivities, from the top down
! to the surface using the delta-eddington solutions for each layer
! 
! method: 
! for more details , see briegleb, bruce p., 1992: delta-eddington
! approximation for solar radiation in the ncar community climate model,
! journal of geophysical research, vol 97, d7, pp7603-7612).
!
! modified for maximum/random cloud overlap by bill collins and john
!    truesdale
! 
! author: bill collins
! 
!-----------------------------------------------------------------------
!  use shr_kind_mod, only: r8 => shr_kind_r8
!  use ppgrid

   implicit none

   integer nspint           ! num of spctrl intervals across solar spectrum

   parameter ( nspint = 19 )
!
! minimum total transmission below which no layer computation are done:
!
   real(r8) trmin                ! minimum total transmission allowed
   real(r8) wray                 ! rayleigh single scatter albedo
   real(r8) gray                 ! rayleigh asymetry parameter
   real(r8) fray                 ! rayleigh forward scattered fraction

   parameter (trmin = 1.e-3)
   parameter (wray = 0.999999)
   parameter (gray = 0.0)
   parameter (fray = 0.1)
!
!------------------------------arguments--------------------------------
!
! input arguments
!
   integer, intent(in) :: pver, pverp, pcols
   real(r8), intent(in) :: coszrs(pcols)        ! cosine zenith angle
   real(r8), intent(in) :: trayoslp             ! tray/sslp
   real(r8), intent(in) :: pflx(pcols,0:pverp)  ! interface pressure
   real(r8), intent(in) :: abh2o                ! absorption coefficiant for h2o
   real(r8), intent(in) :: abo3                 ! absorption coefficiant for o3
   real(r8), intent(in) :: abco2                ! absorption coefficiant for co2
   real(r8), intent(in) :: abo2                 ! absorption coefficiant for o2
   real(r8), intent(in) :: uh2o(pcols,0:pver)   ! layer absorber amount of h2o
   real(r8), intent(in) :: uo3(pcols,0:pver)    ! layer absorber amount of  o3
   real(r8), intent(in) :: uco2(pcols,0:pver)   ! layer absorber amount of co2
   real(r8), intent(in) :: uo2(pcols,0:pver)    ! layer absorber amount of  o2
   real(r8), intent(in) :: tauxcl(pcols,0:pver) ! cloud extinction optical depth (liquid)
   real(r8), intent(in) :: wcl(pcols,0:pver)    ! cloud single scattering albedo (liquid)
   real(r8), intent(in) :: gcl(pcols,0:pver)    ! cloud asymmetry parameter (liquid)
   real(r8), intent(in) :: fcl(pcols,0:pver)    ! cloud forward scattered fraction (liquid)
   real(r8), intent(in) :: tauxci(pcols,0:pver) ! cloud extinction optical depth (ice)
   real(r8), intent(in) :: wci(pcols,0:pver)    ! cloud single scattering albedo (ice)
   real(r8), intent(in) :: gci(pcols,0:pver)    ! cloud asymmetry parameter (ice)
   real(r8), intent(in) :: fci(pcols,0:pver)    ! cloud forward scattered fraction (ice)
   real(r8), intent(in) :: tauxar(pcols,0:pver) ! aerosol extinction optical depth
   real(r8), intent(in) :: wa(pcols,0:pver)     ! aerosol single scattering albedo
   real(r8), intent(in) :: ga(pcols,0:pver)     ! aerosol asymmetry parameter
   real(r8), intent(in) :: fa(pcols,0:pver)     ! aerosol forward scattered fraction

   integer, intent(in) :: ndayc                 ! number of daylight columns
   integer, intent(in) :: idayc(pcols)          ! daylight column indices
   integer, intent(in) :: ns                    ! index of spectral interval
!
! input/output arguments
!
! following variables are defined for each layer; 0 refers to extra
! layer above top of model:
!
   real(r8), intent(inout) :: rdir(nspint,pcols,0:pver)   ! layer reflectivity to direct rad
   real(r8), intent(inout) :: rdif(nspint,pcols,0:pver)   ! layer reflectivity to diffuse rad
   real(r8), intent(inout) :: tdir(nspint,pcols,0:pver)   ! layer transmission to direct rad
   real(r8), intent(inout) :: tdif(nspint,pcols,0:pver)   ! layer transmission to diffuse rad
   real(r8), intent(inout) :: explay(nspint,pcols,0:pver) ! solar beam exp transm for layer
!
! corresponding quantities for clear-skies
!
   real(r8), intent(inout) :: rdirc(nspint,pcols,0:pver)  ! clear layer reflec. to direct rad
   real(r8), intent(inout) :: rdifc(nspint,pcols,0:pver)  ! clear layer reflec. to diffuse rad
   real(r8), intent(inout) :: tdirc(nspint,pcols,0:pver)  ! clear layer trans. to direct rad
   real(r8), intent(inout) :: tdifc(nspint,pcols,0:pver)  ! clear layer trans. to diffuse rad
   real(r8), intent(inout) :: explayc(nspint,pcols,0:pver)! solar beam exp transm clear layer
!
!---------------------------local variables-----------------------------
!
   integer i                 ! column indices
   integer k                 ! level index
   integer nn                ! index of column loops (max=ndayc)

   real(r8) taugab(pcols)        ! layer total gas absorption optical depth
   real(r8) tauray(pcols)        ! layer rayleigh optical depth
   real(r8) taucsc               ! layer cloud scattering optical depth
   real(r8) tautot               ! total layer optical depth
   real(r8) wtot                 ! total layer single scatter albedo
   real(r8) gtot                 ! total layer asymmetry parameter
   real(r8) ftot                 ! total layer forward scatter fraction
   real(r8) wtau                 !  rayleigh layer scattering optical depth
   real(r8) wt                   !  layer total single scattering albedo
   real(r8) ts                   !  layer scaled extinction optical depth
   real(r8) ws                   !  layer scaled single scattering albedo
   real(r8) gs                   !  layer scaled asymmetry parameter
!
!---------------------------statement functions-------------------------
!
! statement functions and other local variables
!
   real(r8) alpha                ! term in direct reflect and transmissivity
   real(r8) gamma                ! term in direct reflect and transmissivity
   real(r8) el                   ! term in alpha,gamma,n,u
   real(r8) taus                 ! scaled extinction optical depth
   real(r8) omgs                 ! scaled single particle scattering albedo
   real(r8) asys                 ! scaled asymmetry parameter
   real(r8) u                    ! term in diffuse reflect and
!    transmissivity
   real(r8) n                    ! term in diffuse reflect and
!    transmissivity
   real(r8) lm                   ! temporary for el
   real(r8) ne                   ! temporary for n
   real(r8) w                    ! dummy argument for statement function
   real(r8) uu                   ! dummy argument for statement function
   real(r8) g                    ! dummy argument for statement function
   real(r8) e                    ! dummy argument for statement function
   real(r8) f                    ! dummy argument for statement function
   real(r8) t                    ! dummy argument for statement function
   real(r8) et                   ! dummy argument for statement function
!
! intermediate terms for delta-eddington solution
!
   real(r8) alp                  ! temporary for alpha
   real(r8) gam                  ! temporary for gamma
   real(r8) ue                   ! temporary for u
   real(r8) arg                  ! exponential argument
   real(r8) extins               ! extinction
   real(r8) amg                  ! alp - gam
   real(r8) apg                  ! alp + gam
!
   alpha(w,uu,g,e) = .75_r8*w*uu*((1._r8 + g*(1._r8-w))/(1._r8 - e*e*uu*uu))
   gamma(w,uu,g,e) = .50_r8*w*((3._r8*g*(1._r8-w)*uu*uu + 1._r8)/(1._r8-e*e*uu*uu))
   el(w,g)         = sqrt(3._r8*(1._r8-w)*(1._r8 - w*g))
   taus(w,f,t)     = (1._r8 - w*f)*t
   omgs(w,f)       = (1._r8 - f)*w/(1._r8 - w*f)
   asys(g,f)       = (g - f)/(1._r8 - f)
   u(w,g,e)        = 1.5_r8*(1._r8 - w*g)/e
   n(uu,et)        = ((uu+1._r8)*(uu+1._r8)/et ) - ((uu-1._r8)*(uu-1._r8)*et)
!
!-----------------------------------------------------------------------
!
! compute layer radiative properties
!
! compute radiative properties (reflectivity and transmissivity for
!    direct and diffuse radiation incident from above, under clear
!    and cloudy conditions) and transmission of direct radiation
!    (under clear and cloudy conditions) for each layer.
!
   do k=0,pver
      do nn=1,ndayc
         i=idayc(nn)
            tauray(i) = trayoslp*(pflx(i,k+1)-pflx(i,k))
            taugab(i) = abh2o*uh2o(i,k) + abo3*uo3(i,k) + abco2*uco2(i,k) + abo2*uo2(i,k)
            tautot = tauxcl(i,k) + tauxci(i,k) + tauray(i) + taugab(i) + tauxar(i,k)
            taucsc = tauxcl(i,k)*wcl(i,k) + tauxci(i,k)*wci(i,k) + tauxar(i,k)*wa(i,k)
            wtau   = wray*tauray(i)
            wt     = wtau + taucsc
            wtot   = wt/tautot
            gtot   = (wtau*gray + gcl(i,k)*wcl(i,k)*tauxcl(i,k) &
                     + gci(i,k)*wci(i,k)*tauxci(i,k) + ga(i,k) *wa(i,k) *tauxar(i,k))/wt
            ftot   = (wtau*fray + fcl(i,k)*wcl(i,k)*tauxcl(i,k) &
                     + fci(i,k)*wci(i,k)*tauxci(i,k) + fa(i,k) *wa(i,k) *tauxar(i,k))/wt
            ts   = taus(wtot,ftot,tautot)
            ws   = omgs(wtot,ftot)
            gs   = asys(gtot,ftot)
            lm   = el(ws,gs)
            alp  = alpha(ws,coszrs(i),gs,lm)
            gam  = gamma(ws,coszrs(i),gs,lm)
            ue   = u(ws,gs,lm)
!
!     limit argument of exponential to 25, in case lm very large:
!
            arg  = min(lm*ts,25._r8)
            extins = exp(-arg)
            ne = n(ue,extins)
            rdif(ns,i,k) = (ue+1._r8)*(ue-1._r8)*(1._r8/extins - extins)/ne
            tdif(ns,i,k)   =   4._r8*ue/ne
!
!     limit argument of exponential to 25, in case coszrs is very small:
!
            arg       = min(ts/coszrs(i),25._r8)
            explay(ns,i,k) = exp(-arg)
            apg = alp + gam
            amg = alp - gam
            rdir(ns,i,k) = amg*(tdif(ns,i,k)*explay(ns,i,k)-1._r8) + apg*rdif(ns,i,k)
            tdir(ns,i,k) = apg*tdif(ns,i,k) + (amg*rdif(ns,i,k)-(apg-1._r8))*explay(ns,i,k)
!
!     under rare conditions, reflectivies and transmissivities can be
!     negative; zero out any negative values
!
            rdir(ns,i,k) = max(rdir(ns,i,k),0.0_r8)
            tdir(ns,i,k) = max(tdir(ns,i,k),0.0_r8)
            rdif(ns,i,k) = max(rdif(ns,i,k),0.0_r8)
            tdif(ns,i,k) = max(tdif(ns,i,k),0.0_r8)
!
!     clear-sky calculation
!
            if (tauxcl(i,k) == 0.0_r8 .and. tauxci(i,k) == 0.0_r8) then

               rdirc(ns,i,k) = rdir(ns,i,k)
               tdirc(ns,i,k) = tdir(ns,i,k)
               rdifc(ns,i,k) = rdif(ns,i,k)
               tdifc(ns,i,k) = tdif(ns,i,k)
               explayc(ns,i,k) = explay(ns,i,k)
            else
               tautot = tauray(i) + taugab(i) + tauxar(i,k)
               taucsc = tauxar(i,k)*wa(i,k)
!
! wtau already computed for all-sky
!
               wt     = wtau + taucsc
               wtot   = wt/tautot
               gtot   = (wtau*gray + ga(i,k)*wa(i,k)*tauxar(i,k))/wt
               ftot   = (wtau*fray + fa(i,k)*wa(i,k)*tauxar(i,k))/wt
               ts   = taus(wtot,ftot,tautot)
               ws   = omgs(wtot,ftot)
               gs   = asys(gtot,ftot)
               lm   = el(ws,gs)
               alp  = alpha(ws,coszrs(i),gs,lm)
               gam  = gamma(ws,coszrs(i),gs,lm)
               ue   = u(ws,gs,lm)
!
!     limit argument of exponential to 25, in case lm very large:
!
               arg  = min(lm*ts,25._r8)
               extins = exp(-arg)
               ne = n(ue,extins)
               rdifc(ns,i,k) = (ue+1._r8)*(ue-1._r8)*(1._r8/extins - extins)/ne
               tdifc(ns,i,k)   =   4._r8*ue/ne
!
!     limit argument of exponential to 25, in case coszrs is very small:
!
               arg       = min(ts/coszrs(i),25._r8)
               explayc(ns,i,k) = exp(-arg)
               apg = alp + gam
               amg = alp - gam
               rdirc(ns,i,k) = amg*(tdifc(ns,i,k)*explayc(ns,i,k)-1._r8)+ &
                               apg*rdifc(ns,i,k)
               tdirc(ns,i,k) = apg*tdifc(ns,i,k) + (amg*rdifc(ns,i,k) - (apg-1._r8))* &
                               explayc(ns,i,k)
!
!     under rare conditions, reflectivies and transmissivities can be
!     negative; zero out any negative values
!
               rdirc(ns,i,k) = max(rdirc(ns,i,k),0.0_r8)
               tdirc(ns,i,k) = max(tdirc(ns,i,k),0.0_r8)
               rdifc(ns,i,k) = max(rdifc(ns,i,k),0.0_r8)
               tdifc(ns,i,k) = max(tdifc(ns,i,k),0.0_r8)
            end if
         end do
   end do

   return
end subroutine raddedmx

subroutine radinp(lchnk   ,ncol    , pcols, pver, pverp,     &
                  pmid    ,pint    ,o3vmr   , pmidrd  ,&
                  pintrd  ,eccf    ,o3mmr   )
!----------------------------------------------------------------------- 
! 
! purpose: 
! set latitude and time dependent arrays for input to solar
! and longwave radiation.
! convert model pressures to cgs, and compute ozone mixing ratio, needed for
! the solar radiation.
! 
! method: 
! <describe the algorithm(s) used in the routine.> 
! <also include any applicable external references.> 
! 
! author: ccm1, cms contact j. kiehl
! 
!-----------------------------------------------------------------------
!  use shr_kind_mod, only: r8 => shr_kind_r8
!  use ppgrid
!  use time_manager, only: get_curr_calday

   implicit none

!------------------------------arguments--------------------------------
!
! input arguments
!
   integer, intent(in) :: lchnk                ! chunk identifier
   integer, intent(in) :: pcols, pver, pverp
   integer, intent(in) :: ncol                 ! number of atmospheric columns

   real(r8), intent(in) :: pmid(pcols,pver)    ! pressure at model mid-levels (pascals)
   real(r8), intent(in) :: pint(pcols,pverp)   ! pressure at model interfaces (pascals)
   real(r8), intent(in) :: o3vmr(pcols,pver)   ! ozone volume mixing ratio
!
! output arguments
!
   real(r8), intent(out) :: pmidrd(pcols,pver)  ! pressure at mid-levels (dynes/cm*2)
   real(r8), intent(out) :: pintrd(pcols,pverp) ! pressure at interfaces (dynes/cm*2)
   real(r8), intent(out) :: eccf                ! earth-sun distance factor
   real(r8), intent(out) :: o3mmr(pcols,pver)   ! ozone mass mixing ratio

!
!---------------------------local variables-----------------------------
!
   integer i                ! longitude loop index
   integer k                ! vertical loop index

   real(r8) :: calday           ! current calendar day
   real(r8) vmmr                ! ozone volume mixing ratio
   real(r8) delta               ! solar declination angle

!
!-----------------------------------------------------------------------
!
!  calday = get_curr_calday()
   eccf = 1. ! declared intent(out) so fill a value (not used in wrf)
!  call shr_orb_decl (calday  ,eccen     ,mvelpp  ,lambm0  ,obliqr  , &
!                     delta   ,eccf)

!
! convert pressure from pascals to dynes/cm2
!
   do k=1,pver
      do i=1,ncol
         pmidrd(i,k) = pmid(i,k)*10.0
         pintrd(i,k) = pint(i,k)*10.0
      end do
   end do
   do i=1,ncol
      pintrd(i,pverp) = pint(i,pverp)*10.0
   end do
!
! convert ozone volume mixing ratio to mass mixing ratio:
!
   vmmr = amo/amd
   do k=1,pver
      do i=1,ncol
         o3mmr(i,k) = vmmr*o3vmr(i,k)
      end do
   end do
!
   return
end subroutine radinp
subroutine radoz2(lchnk   ,ncol    ,pcols, pver, pverp, o3vmr   ,pint    ,plol    ,plos, ntoplw    )
!----------------------------------------------------------------------- 
! 
! purpose: 
! computes the path length integrals to the model interfaces given the
! ozone volume mixing ratio
! 
! method: 
! <describe the algorithm(s) used in the routine.> 
! <also include any applicable external references.> 
! 
! author: ccm1, cms contact j. kiehl
! 
!-----------------------------------------------------------------------
!  use shr_kind_mod, only: r8 => shr_kind_r8
!  use ppgrid
!  use comozp

   implicit none
!------------------------------input arguments--------------------------
!
   integer, intent(in) :: lchnk                ! chunk identifier
   integer, intent(in) :: ncol                 ! number of atmospheric columns
   integer, intent(in) :: pcols, pver, pverp

   real(r8), intent(in) :: o3vmr(pcols,pver)   ! ozone volume mixing ratio
   real(r8), intent(in) :: pint(pcols,pverp)   ! model interface pressures

   integer, intent(in) :: ntoplw               ! topmost level/layer longwave is solved for

!
!----------------------------output arguments---------------------------
!
   real(r8), intent(out) :: plol(pcols,pverp)   ! ozone prs weighted path length (cm)
   real(r8), intent(out) :: plos(pcols,pverp)   ! ozone path length (cm)

!
!---------------------------local workspace-----------------------------
!
   integer i                ! longitude index
   integer k                ! level index
!
!-----------------------------------------------------------------------
!
! evaluate the ozone path length integrals to interfaces;
! factors of .1 and .01 to convert pressures from cgs to mks:
!
   do i=1,ncol
      plos(i,ntoplw) = 0.1 *cplos*o3vmr(i,ntoplw)*pint(i,ntoplw)
      plol(i,ntoplw) = 0.01*cplol*o3vmr(i,ntoplw)*pint(i,ntoplw)*pint(i,ntoplw)
   end do
   do k=ntoplw+1,pverp
      do i=1,ncol
         plos(i,k) = plos(i,k-1) + 0.1*cplos*o3vmr(i,k-1)*(pint(i,k) - pint(i,k-1))
         plol(i,k) = plol(i,k-1) + 0.01*cplol*o3vmr(i,k-1)* &
                    (pint(i,k)*pint(i,k) - pint(i,k-1)*pint(i,k-1))
      end do
   end do
!
   return
end subroutine radoz2


subroutine radozn (lchnk, ncol, pcols, pver,pmid, pin, levsiz, ozmix, o3vmr)
!----------------------------------------------------------------------- 
! 
! purpose: interpolate ozone from current time-interpolated values to model levels
! 
! method: use pressure values to determine interpolation levels
! 
! author: bruce briegleb
! 
!--------------------------------------------------------------------------
!  use shr_kind_mod, only: r8 => shr_kind_r8
!  use ppgrid
!  use phys_grid,     only: get_lat_all_p, get_lon_all_p
!  use comozp
!  use abortutils, only: endrun
!--------------------------------------------------------------------------
   implicit none
!--------------------------------------------------------------------------
!
! arguments
!
   integer, intent(in) :: lchnk               ! chunk identifier
   integer, intent(in) :: pcols, pver
   integer, intent(in) :: ncol                ! number of atmospheric columns
   integer, intent(in) :: levsiz              ! number of ozone layers

   real(r8), intent(in) :: pmid(pcols,pver)   ! level pressures (mks)
   real(r8), intent(in) :: pin(levsiz)        ! ozone data level pressures (mks)
   real(r8), intent(in) :: ozmix(pcols,levsiz) ! ozone mixing ratio

   real(r8), intent(out) :: o3vmr(pcols,pver) ! ozone volume mixing ratio
!
! local storage
!
   integer i                   ! longitude index
   integer k, kk, kkstart      ! level indices
   integer kupper(pcols)       ! level indices for interpolation
   integer kount               ! counter
   integer lats(pcols)         ! latitude indices
   integer lons(pcols)         ! latitude indices

   real(r8) dpu                ! upper level pressure difference
   real(r8) dpl                ! lower level pressure difference
!
! initialize latitude indices
!
!  call get_lat_all_p(lchnk, ncol, lats)
!  call get_lon_all_p(lchnk, ncol, lons)
!
! initialize index array
!
   do i=1,ncol
      kupper(i) = 1
   end do

   do k=1,pver
!
! top level we need to start looking is the top level for the previous k
! for all longitude points
!
      kkstart = levsiz
      do i=1,ncol
         kkstart = min0(kkstart,kupper(i))
      end do
      kount = 0
!
! store level indices for interpolation
!
      do kk=kkstart,levsiz-1
         do i=1,ncol
            if (pin(kk).lt.pmid(i,k) .and. pmid(i,k).le.pin(kk+1)) then
               kupper(i) = kk
               kount = kount + 1
            end if
         end do
!
! if all indices for this level have been found, do the interpolation and
! go to the next level
!
         if (kount.eq.ncol) then
            do i=1,ncol
               dpu = pmid(i,k) - pin(kupper(i))
               dpl = pin(kupper(i)+1) - pmid(i,k)
               o3vmr(i,k) = (ozmix(i,kupper(i))*dpl + &
                             ozmix(i,kupper(i)+1)*dpu)/(dpl + dpu)
            end do
            goto 35
         end if
      end do
!
! if we've fallen through the kk=1,levsiz-1 loop, we cannot interpolate and
! must extrapolate from the bottom or top ozone data level for at least some
! of the longitude points.
!
      do i=1,ncol
         if (pmid(i,k) .lt. pin(1)) then
            o3vmr(i,k) = ozmix(i,1)*pmid(i,k)/pin(1)
         else if (pmid(i,k) .gt. pin(levsiz)) then
            o3vmr(i,k) = ozmix(i,levsiz)
         else
            dpu = pmid(i,k) - pin(kupper(i))
            dpl = pin(kupper(i)+1) - pmid(i,k)
            o3vmr(i,k) = (ozmix(i,kupper(i))*dpl + &
                          ozmix(i,kupper(i)+1)*dpu)/(dpl + dpu)
         end if
      end do

      if (kount.gt.ncol) then
         call endrun ('radozn: bad ozone data: non-monotonicity suspected')
      end if
35    continue
   end do

   return
end subroutine radozn


#endif

end module module_ra_cam
